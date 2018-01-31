program main 

  use data_global 
  use mod_grid 

  implicit none 

  real ( dp ) :: diff = 1.0 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: ae , aw , an , as , at , ab , ap 
  integer :: nicv , njcv , nkcv 
  integer :: level 
  integer :: i , j , k 
  integer :: i2 , j2 , k2 
  real ( dp ) :: aed , awd , asd , and , atd , abd 
  real ( dp ) :: time1 , time2 
  real ( dp ) :: tmpr = 0 
  integer :: tmpi 
  integer , allocatable , dimension ( : ) :: niter 
  integer :: vcycle_count 
  integer :: tmpxxx = 0 
  real ( 8 ) :: frhs , brhs , truesol

  call grid_read ( "grid_binary.out") 

  call cpu_time ( time1 ) 

  MG_LEVELS = 1 

#include "MACRO.F"

  allocate ( niter ( MG_LEVELS ) ) 

  niter ( : ) = 10 
  niter ( 1 ) = 6 
  niter ( MG_LEVELS ) = 100 
  if ( MG_LEVELS == 1 ) niter ( MG_LEVELS ) = 1500

  level = 1 
  nicv = grid ( level ) %nicv 
  njcv = grid ( level ) %njcv 
  nkcv = grid ( level ) %nkcv 

  allocate ( ap ( nicv , njcv , nkcv ) ) 
  allocate ( ae ( nicv , njcv , nkcv ) ) 
  allocate ( aw ( nicv , njcv , nkcv ) ) 
  allocate ( an ( nicv , njcv , nkcv ) ) 
  allocate ( as ( nicv , njcv , nkcv ) ) 
  allocate ( at ( nicv , njcv , nkcv ) ) 
  allocate ( ab ( nicv , njcv , nkcv ) ) 

  do level = 1 , MG_LEVELS 
    grid ( level ) %f = 0 
    grid ( level ) %v = 0 
    grid ( level ) %r = 0 
  enddo 

  !-------------------------------------------------
  vcycle_count = 0 
  tmpr = 1 
  BIGMG : DO WHILE ( abs(tmpr) > 1.e-10 ) 

    vcycle_count = vcycle_count + 1 
    print * , "V CYCLE COUNT:", vcycle_count 

    grid ( 1 ) %f = 0 
    do level = 2 , MG_LEVELS 
      grid ( level ) %f = 0 
      grid ( level ) %v = 0 
    enddo 

    DOWN : do level = 1 , MG_LEVELS 

      print * , ""
      print * , "DOWN LEVEL:", LEVEL 

      nicv = grid ( level ) %nicv 
      njcv = grid ( level ) %njcv 
      nkcv = grid ( level ) %nkcv 

      do k = 1 , nkcv 
        do j = 1 , njcv 
          do i = 1 , nicv 

            ! central difference for diffusion:
            aed =  diff * areax_ ( i , j , k ) / ( xc_ ( i + 1 , j , k ) - xc_ ( i , j , k ) ) ! TODO: this doesnt eq to Lpe !!!
            awd =  diff * areax_ ( i - 1 , j , k ) / ( xc_ ( i , j , k ) - xc_ ( i - 1 , j , k ) ) 
            and =  diff * areay_ ( i , j , k ) / ( yc_ ( i , j + 1 , k ) - yc_ ( i , j , k ) ) 
            asd =  diff * areay_ ( i , j - 1 , k ) / ( yc_ ( i , j , k ) - yc_ ( i , j - 1 , k ) ) 
            atd =  diff * areaz_ ( i , j , k ) / ( zc_ ( i , j , k + 1 ) - zc_ ( i , j , k ) ) 
            abd =  diff * areaz_ ( i , j , k - 1 ) / ( zc_ ( i , j , k ) - zc_ ( i , j , k - 1 ) ) 

            ae ( i , j , k ) = aed 
            aw ( i , j , k ) = awd 
            an ( i , j , k ) = and 
            as ( i , j , k ) = asd 
            at ( i , j , k ) = atd 
            ab ( i , j , k ) = abd 

            ap ( i , j , k ) = - ( aed + awd + and + asd + atd + abd ) 

            if ( level == 1 ) then ! -------------
              grid ( level ) %f ( i , j , k ) = frhs ( xc_ ( i , j , k ) , yc_ ( i , j , k ) , zc_ ( i , j , k ) ) * vol_ ( i , j , k ) 
            endif 

          enddo 
        enddo 
      enddo 

#include "BC.F"

      ! TODO: call smoother
      call solver_gs ( grid ( 1 ) %nicv , grid ( 1 ) %njcv , grid ( 1 ) %nkcv , nicv , njcv , nkcv , ae , aw , an , as , at , ab , ap , grid ( level ) %f , grid ( level ) %v , niter ( level ) ) 
!      call solver_sip ( grid ( 1 ) %nicv , grid ( 1 ) %njcv , grid ( 1 ) %nkcv , nicv , njcv , nkcv , ae , aw , an , as , at , ab , ap , grid ( level ) %f , grid ( level ) %v , niter ( level ) ) 

      ! Restriction
      if ( level /= MG_LEVELS ) then 

        do k = 1 , nkcv 
          do j = 1 , njcv 
            do i = 1 , nicv 
              grid ( level ) %r ( i , j , k ) = grid ( level ) %f ( i , j , k ) - &
              ap ( i , j , k ) * grid ( level ) %v ( i , j , k ) - &
              ae ( i , j , k ) * grid ( level ) %v ( i + 1 , j , k ) - &
              aw ( i , j , k ) * grid ( level ) %v ( i - 1 , j , k ) - &
              an ( i , j , k ) * grid ( level ) %v ( i , j + 1 , k ) - &
              as ( i , j , k ) * grid ( level ) %v ( i , j - 1 , k ) - &
              at ( i , j , k ) * grid ( level ) %v ( i , j , k + 1 ) - &
              ab ( i , j , k ) * grid ( level ) %v ( i , j , k - 1 ) 
            enddo 
          enddo 
        enddo 

        k2 = 1 
        do k = 1 , nkcv , 2 
          j2 = 1 
          do j = 1 , njcv , 2 
            i2 = 1 
            do i = 1 , nicv , 2 
              grid ( level + 1 ) %f ( i2 , j2 , k2 ) = ( grid ( level ) %r ( i , j , k ) + grid ( level ) %r ( i + 1 , j , k ) + &
              grid ( level ) %r ( i , j + 1 , k ) + grid ( level ) %r ( i + 1 , j + 1 , k ) + &
              grid ( level ) %r ( i , j , k + 1 ) + grid ( level ) %r ( i + 1 , j , k + 1 ) + &
              grid ( level ) %r ( i , j + 1 , k + 1 ) + grid ( level ) %r ( i + 1 , j + 1 , k + 1 ) ) / 1.0 
              i2 = i2 + 1 
            enddo 
            j2 = j2 + 1 
          enddo 
          k2 = k2 + 1 
        enddo 
      endif 

    enddo DOWN 

    UP : do level = MG_LEVELS - 1 , 1 , - 1  ! fixme

      print * , ""
      print * , "UP LEVEL:", level 

      nicv = grid ( level ) %nicv 
      njcv = grid ( level ) %njcv 
      nkcv = grid ( level ) %nkcv 

      do k = 1 , nkcv 
        do j = 1 , njcv 
          do i = 1 , nicv 

            ! central difference for diffusion:
            aed = - diff * areax_ ( i , j , k ) / ( xc_ ( i + 1 , j , k ) - xc_ ( i , j , k ) ) ! TODO: this doesnt eq to Lpe !!!
            awd = - diff * areax_ ( i - 1 , j , k ) / ( xc_ ( i , j , k ) - xc_ ( i - 1 , j , k ) ) 
            and = - diff * areay_ ( i , j , k ) / ( yc_ ( i , j + 1 , k ) - yc_ ( i , j , k ) ) 
            asd = - diff * areay_ ( i , j - 1 , k ) / ( yc_ ( i , j , k ) - yc_ ( i , j - 1 , k ) ) 
            atd = - diff * areaz_ ( i , j , k ) / ( zc_ ( i , j , k + 1 ) - zc_ ( i , j , k ) ) 
            abd = - diff * areaz_ ( i , j , k - 1 ) / ( zc_ ( i , j , k ) - zc_ ( i , j , k - 1 ) ) 

            ae ( i , j , k ) = aed 
            aw ( i , j , k ) = awd 
            an ( i , j , k ) = and 
            as ( i , j , k ) = asd 
            at ( i , j , k ) = atd 
            ab ( i , j , k ) = abd 

            ap ( i , j , k ) = - ( aed + awd + and + asd + atd + abd ) 

          enddo 
        enddo 
      enddo 

      ! BC modifications
      aw ( 1 , : , : ) = 0 
      ae ( nicv , : , : ) = 0 
      as ( : , 1 , : ) = 0 
      an ( : , njcv , : ) = 0 
      ab ( : , : , 1 ) = 0 
      at ( : , : , nkcv ) = 0 

      ! TODO: call smoother
      if ( level /= MG_LEVELS ) then 
        call solver_gs ( grid ( 1 ) %nicv , grid ( 1 ) %njcv , grid ( 1 ) %nkcv , nicv , njcv , nkcv , ae , aw , an , as , at , ab , ap , grid ( level ) %f , grid ( level ) %v , niter ( level ) ) 
        !call solver_sip ( grid ( 1 ) %nicv , grid ( 1 ) %njcv , grid ( 1 ) %nkcv , nicv , njcv , nkcv , ae , aw , an , as , at , ab , ap , grid ( level ) %f , grid ( level ) %v , niter ( level ) ) 
      endif 

      ! Prolongation
      if ( level /= 1 ) then 
        k2 = 1 
        do k = 1 , nkcv 
          j2 = 1 
          do j = 1 , njcv 
            i2 = 1 
            do i = 1 , nicv 
              ! lower part
              grid ( level - 1 ) %v ( i2 , j2 , k2 ) = grid ( level - 1 ) %v ( i2 , j2 , k2 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i , j - 1 , k ) + &
              9 * grid ( level ) %v ( i - 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j , k - 1 ) + &
              3 * grid ( level ) %v ( i , j - 1 , k - 1 ) + &
              3 * grid ( level ) %v ( i - 1 , j - 1 , k ) + &
              3 * grid ( level ) %v ( i - 1 , j , k - 1 ) + &
              grid ( level ) %v ( i - 1 , j - 1 , k - 1 ) &
              ) 

              grid ( level - 1 ) %v ( i2 + 1 , j2 , k2 ) = grid ( level - 1 ) %v ( i2 + 1 , j2 , k2 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i , j - 1 , k ) + &
              9 * grid ( level ) %v ( i + 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j , k - 1 ) + &
              3 * grid ( level ) %v ( i , j - 1 , k - 1 ) + &
              3 * grid ( level ) %v ( i + 1 , j , k - 1 ) + &
              3 * grid ( level ) %v ( i + 1 , j - 1 , k ) + &
              grid ( level ) %v ( i + 1 , j - 1 , k - 1 ) &
              ) 

              grid ( level - 1 ) %v ( i2 + 1 , j2 + 1 , k2 ) = grid ( level - 1 ) %v ( i2 + 1 , j2 + 1 , k2 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i + 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j + 1 , k ) + &
              9 * grid ( level ) %v ( i , j , k - 1 ) + &
              3 * grid ( level ) %v ( i + 1 , j + 1 , k ) + &
              3 * grid ( level ) %v ( i + 1 , j , k - 1 ) + &
              3 * grid ( level ) %v ( i , j + 1 , k - 1 ) + &
              grid ( level ) %v ( i + 1 , j + 1 , k - 1 ) &
              ) 

              grid ( level - 1 ) %v ( i2 , j2 + 1 , k2 ) = grid ( level - 1 ) %v ( i2 , j2 + 1 , k2 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i , j , k - 1 ) + &
              9 * grid ( level ) %v ( i , j + 1 , k ) + &
              9 * grid ( level ) %v ( i - 1 , j , k ) + &
              3 * grid ( level ) %v ( i - 1 , j + 1 , k ) + &
              3 * grid ( level ) %v ( i - 1 , j , k - 1 ) + &
              3 * grid ( level ) %v ( i , j + 1 , k - 1 ) + &
              grid ( level ) %v ( i - 1 , j + 1 , k - 1 ) &
              ) 

              ! upper part
              grid ( level - 1 ) %v ( i2 , j2 , k2 + 1 ) = grid ( level - 1 ) %v ( i2 , j2 , k2 + 1 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i , j - 1 , k ) + &
              9 * grid ( level ) %v ( i - 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j , k + 1 ) + &
              3 * grid ( level ) %v ( i , j - 1 , k + 1 ) + &
              3 * grid ( level ) %v ( i - 1 , j , k + 1 ) + &
              3 * grid ( level ) %v ( i - 1 , j - 1 , k ) + &
              grid ( level ) %v ( i - 1 , j - 1 , k + 1 ) &
              ) 

              grid ( level - 1 ) %v ( i2 + 1 , j2 , k2 + 1 ) = grid ( level - 1 ) %v ( i2 + 1 , j2 , k2 + 1 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i , j , k + 1 ) + &
              9 * grid ( level ) %v ( i + 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j - 1 , k ) + &
              3 * grid ( level ) %v ( i + 1 , j - 1 , k ) + &
              3 * grid ( level ) %v ( i + 1 , j , k + 1 ) + &
              3 * grid ( level ) %v ( i , j - 1 , k + 1 ) + &
              grid ( level ) %v ( i + 1 , j - 1 , k + 1 ) &
              ) 

              grid ( level - 1 ) %v ( i2 + 1 , j2 + 1 , k2 + 1 ) = grid ( level - 1 ) %v ( i2 + 1 , j2 + 1 , k2 + 1 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i + 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j + 1 , k ) + &
              9 * grid ( level ) %v ( i , j , k + 1 ) + &
              3 * grid ( level ) %v ( i + 1 , j , k + 1 ) + &
              3 * grid ( level ) %v ( i , j + 1 , k + 1 ) + &
              3 * grid ( level ) %v ( i + 1 , j + 1 , k ) + &
              grid ( level ) %v ( i + 1 , j + 1 , k + 1 ) &
              ) 

              grid ( level - 1 ) %v ( i2 , j2 + 1 , k2 + 1 ) = grid ( level - 1 ) %v ( i2 , j2 + 1 , k2 + 1 ) + &
              ( 1.0 / 64.0 ) * ( &
              27 * grid ( level ) %v ( i , j , k ) + &
              9 * grid ( level ) %v ( i , j , k + 1 ) + &
              9 * grid ( level ) %v ( i - 1 , j , k ) + &
              9 * grid ( level ) %v ( i , j + 1 , k ) + &
              3 * grid ( level ) %v ( i , j + 1 , k + 1 ) + &
              3 * grid ( level ) %v ( i - 1 , j + 1 , k ) + &
              3 * grid ( level ) %v ( i - 1 , j , k + 1 ) + &
              grid ( level ) %v ( i - 1 , j + 1 , k + 1 ) &
              ) 

              i2 = i2 + 2 
            enddo 
            j2 = j2 + 2 
          enddo 
          k2 = k2 + 2 
        enddo 

      endif 

    enddo UP 
    !-------------------------------------------------
    ! Compute R.M.S.
    level = 1 
    do k = 1 , nkcv 
      do j = 1 , njcv 
        do i = 1 , nicv 
          ! central difference for diffusion:
          aed =  diff * areax_ ( i , j , k ) / ( xc_ ( i + 1 , j , k ) - xc_ ( i , j , k ) ) ! TODO: this doesnt eq to Lpe !!!
          awd =  diff * areax_ ( i - 1 , j , k ) / ( xc_ ( i , j , k ) - xc_ ( i - 1 , j , k ) ) 
          and =  diff * areay_ ( i , j , k ) / ( yc_ ( i , j + 1 , k ) - yc_ ( i , j , k ) ) 
          asd =  diff * areay_ ( i , j - 1 , k ) / ( yc_ ( i , j , k ) - yc_ ( i , j - 1 , k ) ) 
          atd =  diff * areaz_ ( i , j , k ) / ( zc_ ( i , j , k + 1 ) - zc_ ( i , j , k ) ) 
          abd =  diff * areaz_ ( i , j , k - 1 ) / ( zc_ ( i , j , k ) - zc_ ( i , j , k - 1 ) ) 
          ae ( i , j , k ) = aed 
          aw ( i , j , k ) = awd 
          an ( i , j , k ) = and 
          as ( i , j , k ) = asd 
          at ( i , j , k ) = atd 
          ab ( i , j , k ) = abd 
          ap ( i , j , k ) = - ( aed + awd + and + asd + atd + abd ) 
            if ( level == 1 ) then ! -------------
              grid ( level ) %f ( i , j , k ) = frhs ( xc_ ( i , j , k ) , yc_ ( i , j , k ) , zc_ ( i , j , k ) ) * vol_ ( i , j , k ) 
            endif 
        enddo 
      enddo 
    enddo 

#include "BC.F"

level=1
    tmpr = 0 
    do k = 1 , nkcv 
      do j = 1 , njcv 
        do i = 1 , nicv 
          grid ( level ) %r ( i , j , k ) = grid ( level ) %f ( i , j , k ) - &
          ap ( i , j , k ) * grid ( level ) %v ( i , j , k ) - &
          ae ( i , j , k ) * grid ( level ) %v ( i + 1 , j , k ) - &
          aw ( i , j , k ) * grid ( level ) %v ( i - 1 , j , k ) - &
          an ( i , j , k ) * grid ( level ) %v ( i , j + 1 , k ) - &
          as ( i , j , k ) * grid ( level ) %v ( i , j - 1 , k ) - &
          at ( i , j , k ) * grid ( level ) %v ( i , j , k + 1 ) - &
          ab ( i , j , k ) * grid ( level ) %v ( i , j , k - 1 ) 
          if ( abs(tmpr) < abs(grid ( level ) %r ( i , j , k )) ) tmpr = ( grid ( level ) %r ( i , j , k ) )
        enddo 
      enddo 
    enddo 
    print * , "MAX RES:", tmpr 
    if ( vCycle_count == 1 ) then 
      exit 
    endif 
  enddo BIGMG 
  !-------------------------------------------------
  ! Compute R.M.S.
  112 level = 1 
  tmpr = 0 
  do k = 1 , nkcv 
    do j = 1 , njcv 
      do i = 1 , nicv 
        grid ( level ) %r ( i , j , k ) = grid ( level ) %f ( i , j , k ) - &
        ap ( i , j , k ) * grid ( level ) %v ( i , j , k ) - &
        ae ( i , j , k ) * grid ( level ) %v ( i + 1 , j , k ) - &
        aw ( i , j , k ) * grid ( level ) %v ( i - 1 , j , k ) - &
        an ( i , j , k ) * grid ( level ) %v ( i , j + 1 , k ) - &
        as ( i , j , k ) * grid ( level ) %v ( i , j - 1 , k ) - &
        at ( i , j , k ) * grid ( level ) %v ( i , j , k + 1 ) - &
        ab ( i , j , k ) * grid ( level ) %v ( i , j , k - 1 ) 
        tmpr = tmpr + ( grid ( level ) %r ( i , j , k ) ) ** 2 
      enddo 
    enddo 
  enddo 
  tmpr = sqrt ( tmpr / ( nicv * njcv * nkcv ) ) 
  print * , "RMS=", tmpr 
  call cpu_time ( time2 ) 
  print * , "MG TIME:", time2 - time1 
  print * , "TOTAL NITER:", total_niter 
  call tecwrite ( 1 ) 

endprogram main 

subroutine tecwrite ( level ) 

  use mod_grid 

  implicit none 
  integer :: i , j , k 
  integer , intent ( in ) :: level 
reAL(8) :: truesol

  write ( * , "(A41)") "WRITING TECPLOT FILE: A_OUTPUT_TECPLOT.plt"
  open ( unit = 102 , status = 'replace ', file = "A_OUTPUT_TEC.plt") 
  write ( 102 , * ) 'TITLE = "Example: Simple 3D-Volume Data"'
  write ( 102 , * ) 'VARIABLES = "X", "Y", "Z", "V", "f" "Error"'
  write ( 102 , * ) 'ZONE I = ', grid(level)%nicv + 0 , 'J = ', grid(level)%njcv + 0 , 'K = ', grid(level)%nkcv + 0 , 'F = POINT '
  do k = 1 , grid ( level ) %nkcv 
    do j = 1 , grid ( level ) %njcv 
      do i = 1 , grid ( level ) %nicv 
        write ( 102 , * ) grid ( level ) %xc ( i , j , k ) , grid ( level ) %yc ( i , j , k ) , grid ( level ) %zc ( i , j , k ) , &
        grid ( level ) %v ( i , j , k ) , grid ( level ) %f ( i , j , k ) , grid ( level ) %v ( i , j , k ) - &
truesol(xc_(i,j,k) , yc_(i,j,k),zc_(i,j,k))
        !if ( grid ( level ) %v ( i , j , k ) /= 3.0 ) then 
        !  print * , i , j , k , grid ( level ) %v ( i , j , k ) , grid ( level ) %nkcv , grid ( level ) %njcv , grid ( level ) %nicv 
        !stop 115
        !endif 
      enddo 
    enddo 
  enddo 
  flush ( 102 ) 
  close ( 102 ) 
  print * , "DONE WRITING..."
  print * , ""
  continue 

endsubroutine tecwrite 

!===================================================================================================
!
! Implements Stone's "Strongly Implicit Procedure - SIP"
!
!===================================================================================================
subroutine solver_sip ( nicv0 , njcv0 , nkcv0 , nicv , njcv , nkcv , ae , aw , an , as , at , ab , ap , q , fi , niter ) 

  use data_global 

  implicit none 

  integer :: nicv0 , njcv0 , nkcv0 
  integer :: nicv , njcv , nkcv 
  real ( dp ) , dimension ( nicv0 , njcv0 , nkcv0 ) , intent ( in ) :: ae , aw , an , as , at , ab , ap 
  real ( dp ) , dimension ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) , intent ( out ) :: fi , q 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: BB , BT , BW , BE , BS , BN , BP 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: RES 
  real ( dp ) :: RSM , RES1 , RESOR 
  real ( dp ) :: P1 , P2 , P3 
  integer :: i , j , k , n 
  real ( dp ) , parameter :: alfa = 0.92 
  integer :: sip_count 
  real ( dp ) :: time1 , time2 
  integer , intent ( in ) :: niter 

  !print *,"SIP:", "0:" , nicv0, njcv0,nkcv0
  !print *,"SIP:", ":" , nicv, njcv,nkcv
  if ( nicv == 128 ) total_niter = total_niter + niter 

  allocate ( BB ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( BT ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( BN ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( BS ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( BE ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( BW ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( BP ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 
  allocate ( RES ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) ) 

  BB = 0 
  BT = 0 
  BN = 0 
  BS = 0 
  BE = 0 
  BW = 0 
  RES = 0 
  RSM = 0 

  call cpu_time ( time1 ) 

  ! calculate coefficients of L and U matrices
  do k = 1 , nkcv 
    do i = 1 , nicv 
      do j = 1 , njcv 

        BB ( i , j , k ) = AB ( i , j , k ) / ( 1. + ALFA * ( BN ( i , j , k - 1 ) + BE ( i , j , k - 1 ) ) ) 
        BW ( i , j , k ) = AW ( i , j , k ) / ( 1. + ALFA * ( BN ( i - 1 , j , k ) + BT ( i - 1 , j , k ) ) ) 
        BS ( i , j , k ) = AS ( i , j , k ) / ( 1. + ALFA * ( BE ( i , j - 1 , k ) + BT ( i , j - 1 , k ) ) ) 

        P1 = ALFA * ( BB ( i , j , k ) * BN ( i , j , k - 1 ) + BW ( i , j , k ) * BN ( i - 1 , j , k ) ) 
        P2 = ALFA * ( BB ( i , j , k ) * BE ( i , j , k - 1 ) + BS ( i , j , k ) * BE ( i , j - 1 , k ) ) 
        P3 = ALFA * ( BW ( i , j , k ) * BT ( i - 1 , j , k ) + BS ( i , j , k ) * BT ( i , j - 1 , k ) ) 

        BP ( i , j , k ) = 1. / ( AP ( i , j , k ) + P1 + P2 + P3 - &
        BB ( i , j , k ) * BT ( i , j , k - 1 ) - &
        BW ( i , j , k ) * BE ( i - 1 , j , k ) - &
        BS ( i , j , k ) * BN ( i , j - 1 , k ) + SMALL ) 

        BN ( i , j , k ) = ( AN ( i , j , k ) - P1 ) * BP ( i , j , k ) 
        BE ( i , j , k ) = ( AE ( i , j , k ) - P2 ) * BP ( i , j , k ) 
        BT ( i , j , k ) = ( AT ( i , j , k ) - P3 ) * BP ( i , j , k ) 

      enddo 
    enddo 
  enddo 

  sip_count = 0 

  SIP_ITERATIONS : do n = 1 , niter 

    sip_count = sip_count + 1 

    res1 = 0.0 

    do k = 1 , nkcv 
      do i = 1 , nicv 
        do j = 1 , njcv 

          RES ( i , j , k ) = - ( AE ( i , j , k ) * FI ( i + 1 , j , k ) + AW ( i , j , k ) * FI ( i - 1 , j , k ) + &
          AN ( i , j , k ) * FI ( i , j + 1 , k ) + AS ( i , j , k ) * FI ( i , j - 1 , k ) + &
          AT ( i , j , k ) * FI ( i , j , k + 1 ) + AB ( i , j , k ) * FI ( i , j , k - 1 ) ) + &
          Q ( i , j , k ) - AP ( i , j , k ) * FI ( i , j , k ) 

          RES1 = RES1 + DABS ( RES ( i , j , k ) ) 

          RES ( i , j , k ) = ( RES ( i , j , k ) - BB ( i , j , k ) * RES ( i , j , k - 1 ) - &
          BW ( i , j , k ) * RES ( i - 1 , j , k ) - &
          BS ( i , j , k ) * RES ( i , j - 1 , k ) ) * BP ( i , j , k ) 

        enddo 
      enddo 
    enddo 

    IF ( N == 1 ) then 
      RESOR = RES1 
    endif 

    RSM = RES1 / ( RESOR + SMALL ) 

    do k = nkcv , 1 , - 1 
      do i = nicv , 1 , - 1 
        do j = njcv , 1 , - 1 

          RES ( i , j , k ) = RES ( i , j , k ) - BN ( i , j , k ) * RES ( i , j + 1 , k ) - &
          BE ( i , j , k ) * RES ( i + 1 , j , k ) - &
          BT ( i , j , k ) * RES ( i , j , k + 1 ) 

          FI ( i , j , k ) = FI ( i , j , k ) + RES ( i , j , k ) 

        enddo 
      enddo 
    enddo 

    if ( RSM < 0.01 ) then 
      call cpu_time ( time2 ) 
      write ( * , "(A28,I4,2X,A11,F10.7,2X,A15,F7.3)") "NUMBER OF SIP  ITERATIONS = ", sip_count , &
      "RESIDUAL = ", rsm, "ELAPSED TIME = ", time2 - time1 
      goto 1000 
    endif 

  enddo SIP_ITERATIONS 

  call cpu_time ( time2 ) 
  write ( * , * ) "***ERROR: NOT CONVERGED"
  write ( * , * ) "***ERROR: NOT CONVERGED"
  write ( * , "(A28,I4,2X,A11,F10.7,2X,A15,F7.3)") "NUMBER OF SIP  ITERATIONS = ", sip_count , &
  "RESIDUAL = ", rsm, "ELAPSED TIME = ", time2 - time1 

  1000 continue 

  deallocate ( Bt ) 
  deallocate ( Bw ) 
  deallocate ( Be ) 
  deallocate ( Bs ) 
  deallocate ( Bn ) 
  deallocate ( Bp ) 
  deallocate ( res ) 

endsubroutine solver_sip 
!===================================================================================================
!
! Implements Stone's "Strongly Implicit Procedure - SIP"
!
!===================================================================================================
subroutine solver_gs ( nicv0 , njcv0 , nkcv0 , nicv , njcv , nkcv , ae , aw , an , as , at , ab , ap , q , fi , niter ) 

  use data_global 

  implicit none 

  integer , intent ( in ) :: nicv0 , njcv0 , nkcv0 
  integer , intent ( in ) :: nicv , njcv , nkcv 
  integer , intent ( in ) :: niter 
  real ( dp ) , dimension ( nicv0 , njcv0 , nkcv0 ) , intent ( in ) :: ae , aw , an , as , at , ab , ap 
  real ( dp ) , dimension ( 0 : nicv + 1 , 0 : njcv + 1 , 0 : nkcv + 1 ) , intent ( inout ) :: fi , q 
  integer :: i , j , k , n 

  if ( nicv == 128 ) total_niter = total_niter + niter 

  SIP_ITERATIONS : do n = 1 , niter 

    do k = 1 , nkcv 
      do i = 1 , nicv 
        do j = 1 , njcv 

          FI ( i , j , k ) = ( 1. / AP ( i , j , k ) ) * ( - ( AE ( i , j , k ) * FI ( i + 1 , j , k ) + AW ( i , j , k ) * FI ( i - 1 , j , k ) + &
          AN ( i , j , k ) * FI ( i , j + 1 , k ) + AS ( i , j , k ) * FI ( i , j - 1 , k ) + &
          AT ( i , j , k ) * FI ( i , j , k + 1 ) + AB ( i , j , k ) * FI ( i , j , k - 1 ) ) + &
          Q ( i , j , k ) ) 

        enddo 
      enddo 
    enddo 
  enddo sip_iterations 
  1000 continue 

endsubroutine solver_gs 

real ( 8 ) function frhs ( x , y , z ) 
  use data_global 
  implicit none 
  real ( dp ) , intent ( in ) :: x , y , z 
real(8) :: pi
  pi=3.14
  frhs = sin(2*pi*x/0.125)*sin(2*pi*y/0.25)*sin(2*pi*z/0.25) 
endfunction frhs 

REAL ( 8 ) FUNCTION BRHS ( ISIDE , X , Y , Z ) 
  implicit none 
  INTEGER ISIDE 
  REAL ( 8 ) X , Y , Z 
  !
  REAL COS , EXP , SIN 
  INTRINSIC COS , EXP , SIN 
  ! Boundary conditions

  IF ( ISIDE .EQ. 5 ) THEN 
    BRHS = - 2.0 * SIN ( 3.0 * X + Y - 2.0 * Z ) - EXP ( X - Z ) 
  ELSE 
    BRHS = COS ( 3.0 * X + Y - 2.0 * Z ) + EXP ( X - Z ) + 1.0 
  END IF 
  RETURN 
ENDfunction brhs 

REAL ( 8 ) FUNCTION truesol ( X , Y , Z ) 
  implicit none 
  REAL ( 8 ) X , Y , Z 
  !
  REAL COS , EXP , SIN 
  INTRINSIC COS , EXP , SIN 

 truesol = COS ( 3.0 * X + Y - 2.0 * Z ) + EXP ( X - Z ) + 1.0 

  RETURN 
ENDfunction truesol
