!===================================================================================================
!
! Global data: precission, small & great number and etc.
!
!===================================================================================================
module data_global 

  implicit none 

  save 

  integer :: total_niter = 0 


  integer , parameter :: dp = 8 
  real ( dp ) :: small 
  real ( dp ) :: great 
  real ( dp ) , parameter :: zero = 0.0 

  integer , dimension ( 10 ) :: ns 
  real ( dp ) , dimension ( 10 ) :: urf 
  real ( dp ) , dimension ( 3 ) :: grad_deff_corr = (/ 0.0 , 0.0 , 0.0 /) 
  real ( dp ) , dimension ( 10 ) :: resmax 

  integer , parameter :: iu = 1 
  integer , parameter :: iv = 2 
  integer , parameter :: iw = 3 
  integer , parameter :: ip = 4 
  integer , parameter :: ipp = 5 
  integer , parameter :: ien = 6 
  integer , parameter :: ite = 7 
  integer , parameter :: ied = 8 

  logical :: LCDS , LSMART 
  real ( dp ) , dimension ( 10 ) :: gama 

  type Vector3D 
  real ( dp ) :: x , y , z 
  endtype Vector3D 

endmodule data_global 
