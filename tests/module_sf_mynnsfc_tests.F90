! Module for MYNN SFC scheme tests
module module_ccpp_mynn_sfc_tests
  implicit none
  contains
    !=================================================================================================================    

    subroutine init_mynn_sfc_flags_for_test_all_true()
      write(*,*) '--- calling  init_mynn_sfc_flags_for_test_all_true()'
    end subroutine init_mynn_sfc_flags_for_test_all_true
    !=================================================================================================================    
    subroutine init_ccpp_data_for_test(nlev)
      integer, intent(in) :: nlev
      write(*,*) '--- calling  opening data file(nlev)'
      ! opening data file

    end subroutine init_ccpp_data_for_test
    !=================================================================================================================
    subroutine ccpp_test()
      use module_sf_mynnsfc_driver

      write(*,*) '--- entring ccpp_test subroutine'
      ! Initialize input data for tests
      call init_ccpp_data_for_test(nlev=59)

      
      ! Initialize MYNN SFC flags for tests
      call init_mynn_sfc_flags_for_test_all_true()

      ! Initialize MYNN SFC
      write(*,*) '--- calling  module_sf_mynnsfc_driver_init()'
    end subroutine ccpp_test
    !=================================================================================================================
    
end module module_ccpp_mynn_sfc_tests           
