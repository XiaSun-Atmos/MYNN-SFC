! ###########################################################################################
! CI test driver for MYNN SFC scheme
! ###########################################################################################
program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use module_ccpp_mynn_sfc_tests,  only: ccpp_test  
  implicit none
  call ccpp_test()
end program driver
