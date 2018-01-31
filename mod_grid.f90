!===================================================================================================
! Multigrid adopted grid
!===================================================================================================
module mod_grid 

  use data_global 

  implicit none 

  save 

  integer :: MG_LEVELS 

  type grid_t 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: v , f , r 
  integer :: nicv , njcv , nkcv 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: x , y , z 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: xc , yc , zc 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: vol 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: areax , areay , areaz 
  type ( Vector3D ) , allocatable , dimension ( : , : , : ) :: s1 , s2 , s3 
  real ( dp ) , allocatable , dimension ( : , : , : ) :: fx , fy , fz 
  endtype grid_t 

  type ( grid_t ) , dimension ( : ) , allocatable :: grid 

  contains 

  !===================================================================================================
  ! Allocate
  !===================================================================================================
  subroutine grid_read ( filename ) 

    implicit none 

    integer :: nicv , njcv , nkcv 
    integer :: i , j , k 
    integer :: level 
    real ( dp ) :: l 
    character ( * ) , intent ( in ) :: filename 

    open ( unit = 100 , status = 'old ', file = filename , form = 'unformatted ') 

    MG_LEVELS = 3 

    allocate ( grid ( MG_LEVELS ) ) 

    ! Allocate grid

    read ( 100 ) grid ( 1 ) %NICV , grid ( 1 ) %NJCV , grid ( 1 ) %NKCV 

    do level = 2 , MG_LEVELS 
      grid ( level ) %nicv = grid ( level - 1 ) %nicv / 2 
      grid ( level ) %njcv = grid ( level - 1 ) %njcv / 2 
      grid ( level ) %nkcv = grid ( level - 1 ) %nkcv / 2 
    enddo 

    do level = 1 , MG_LEVELS 

      nicv = grid ( level ) %nicv 
      njcv = grid ( level ) %njcv 
      nkcv = grid ( level ) %nkcv 

      allocate ( grid ( level ) %x ( - 1 : NICV + 1 , - 1 : NJCV + 1 , - 1 : NKCV + 1 ) ) 
      allocate ( grid ( level ) %y ( - 1 : NICV + 1 , - 1 : NJCV + 1 , - 1 : NKCV + 1 ) ) 
      allocate ( grid ( level ) %z ( - 1 : NICV + 1 , - 1 : NJCV + 1 , - 1 : NKCV + 1 ) ) 

      allocate ( grid ( level ) %xc ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 
      allocate ( grid ( level ) %yc ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 
      allocate ( grid ( level ) %zc ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 

      allocate ( grid ( level ) %vol ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 

      allocate ( grid ( level ) %v ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 
      allocate ( grid ( level ) %f ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 
      allocate ( grid ( level ) %r ( 0 : NICV + 1 , 0 : NJCV + 1 , 0 : NKCV + 1 ) ) 

      allocate ( grid ( level ) %areax ( 0 : NICV , NJCV , NKCV ) ) 
      allocate ( grid ( level ) %areay ( NICV , 0 : NJCV , NKCV ) ) 
      allocate ( grid ( level ) %areaz ( NICV , NJCV , 0 : NKCV ) ) 

      allocate ( grid ( level ) %s1 ( 0 : NICV , NJCV , NKCV ) ) 
      allocate ( grid ( level ) %s2 ( NICV , 0 : NJCV , NKCV ) ) 
      allocate ( grid ( level ) %s3 ( NICV , NJCV , 0 : NKCV ) ) 

      allocate ( grid ( level ) %fx ( 0 : NICV , NJCV , NKCV ) ) 
      allocate ( grid ( level ) %fy ( NICV , 0 : NJCV , NKCV ) ) 
      allocate ( grid ( level ) %fz ( NICV , NJCV , 0 : NKCV ) ) 

    enddo 

    ! Let's reopen gridfile
    close ( 100 ) 
    open ( unit = 100 , status = 'old ', file = filename , form = 'unformatted ') 
    !-------------------------------------------------
    MG_GRID : do level = 1 , MG_LEVELS 

      write ( * , * ) "LEVEL:", level 
      read ( 100 ) grid ( level ) %NICV , grid ( level ) %NJCV , grid ( level ) %NKCV 
      write ( * , * ) "GRID:", grid ( level ) %NICV , grid ( level ) %NJCV , grid ( level ) %NKCV 

      nicv = grid ( level ) %nicv 
      njcv = grid ( level ) %njcv 
      nkcv = grid ( level ) %nkcv 

#define x_ grid(level)%x
#define y_ grid(level)%y
#define z_ grid(level)%z

#define xc_ grid(level)%xc
#define yc_ grid(level)%yc
#define zc_ grid(level)%zc

#define vol_ grid(level)%vol 

#define areax_ grid(level)%areax
#define areay_ grid(level)%areay
#define areaz_ grid(level)%areaz

#define s1_ grid(level)%s1 
#define s2_ grid(level)%s2 
#define s3_ grid(level)%s3 

#define fx_ grid(level)%fx 
#define fy_ grid(level)%fy 
#define fz_ grid(level)%fz 

      do i = - 1 , nicv + 1 
        do j = - 1 , njcv + 1 
          do k = - 1 , nkcv + 1 
            read ( 100 ) x_ ( i , j , k ) , y_ ( i , j , k ) , z_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      do i = 0 , nicv + 1 
        do j = 0 , njcv + 1 
          do k = 0 , nkcv + 1 
            read ( 100 ) xc_ ( i , j , k ) , yc_ ( i , j , k ) , zc_ ( i , j , k ) , vol_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      !-------------------------------------------------
      ! Read cell face area
      !-------------------------------------------------

      do i = 0 , nicv 
        do j = 1 , njcv 
          do k = 1 , nkcv 
            read ( 100 ) areax_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      do i = 1 , nicv 
        do j = 0 , njcv 
          do k = 1 , nkcv 
            read ( 100 ) areay_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      do i = 1 , nicv 
        do j = 1 , njcv 
          do k = 0 , nkcv 
            read ( 100 ) areaz_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      !-------------------------------------------------
      ! Read cell face normal vector
      !-------------------------------------------------

      do i = 0 , nicv 
        do j = 1 , njcv 
          do k = 1 , nkcv 
            read ( 100 ) s1_ ( i , j , k ) %x , s1_ ( i , j , k ) %y , s1_ ( i , j , k ) %z 
          enddo 
        enddo 
      enddo 

      do i = 1 , nicv 
        do j = 0 , njcv 
          do k = 1 , nkcv 
            read ( 100 ) s2_ ( i , j , k ) %x , s2_ ( i , j , k ) %y , s2_ ( i , j , k ) %z 
          enddo 
        enddo 
      enddo 

      do i = 1 , nicv 
        do j = 1 , njcv 
          do k = 0 , nkcv 
            read ( 100 ) s3_ ( i , j , k ) %x , s3_ ( i , j , k ) %y , s3_ ( i , j , k ) %z 
          enddo 
        enddo 
      enddo 

      do i = 0 , nicv 
        do j = 1 , njcv 
          do k = 1 , nkcv 
            read ( 100 ) fx_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      do i = 1 , nicv 
        do j = 0 , njcv 
          do k = 1 , nkcv 
            read ( 100 ) fy_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

      do i = 1 , nicv 
        do j = 1 , njcv 
          do k = 0 , nkcv 
            read ( 100 ) fz_ ( i , j , k ) 
          enddo 
        enddo 
      enddo 

    enddo MG_GRID 

    print * , ""
    print * , "DONE READING GRID"
    print * , ""

  endsubroutine grid_read 

endmodule mod_grid 

