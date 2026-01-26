! ###########################################################################################
! CI test driver for MYNN SFC scheme
! ###########################################################################################
program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use module_sf_mynnsfc_driver,  only: psi_init
  use module_sf_mynnsfc_ccpp_tests
  use module_sf_mynnsfc_wrf_tests
  implicit none
   !Global configuration options, to be moved to namelist variables:
 integer, parameter :: psi_opt         = 0       !0: mynn, 1: gfs
   real, dimension(0:1000 ),save :: psim_stab,psim_unstab, &
                                             psih_stab,psih_unstab
                                              !--- output arguments:
 character(len=512) :: errmsg
 integer :: errflg

  errmsg = ' '
  errflg = 0
  call psi_init(psi_opt,errmsg, errflg )
  write(*,*) 'psim_stab ',psim_stab 
  call ccpp_test()
  call wrf_test()

end program driver
