! ###########################################################################################
! Main program
! ###########################################################################################
program driver
  use machine, only: wp => kind_phys, sp => kind_sngl_prec, dp => kind_dbl_prec
  use scheme,  only: scheme_run
  
  implicit none

  character(len=256) :: file_datain
  character(len=256) :: file_dataout
  
  ! Read in sample data
  file_datain = 'input_data.nc'
  call read_sample_data(file_datain)
  
  ! Call scheme
  call scheme_run()

  ! Write output
  file_dataout = 'driver_data.nc'
  call write_output(file_dataout)

contains
  ! #########################################################################################
  ! Subroutine to read sample data for offline driver.
  ! #########################################################################################
  subroutine read_sample_data(file)
    use netcdf
    character(len=256), intent(in) :: file
    !
    integer :: ncid, status
    
    print*,'reading in sample data from ',trim(file)

    ! Open file
    status = nf90_open(path=trim(file), mode=nf90_nowrite, ncid=ncid)
    if (status /= nf90_noerr) call handle_err(status)

  end subroutine read_sample_data

  ! #########################################################################################
  ! Subroutine to write output to file (NetCDF)
  ! #########################################################################################
  subroutine write_output(file)
    use netcdf
    character(len=256), intent(in) :: file
    !
    integer :: ncid, status

    print*,'writing driver output to file ',trim(file)

    ! Create file
    status = nf90_create(path=trim(file), cmode=nf90_clobber, ncid=ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! Define dimension

    ! Define variables

    ! Exit define mode
    status = nf90_enddef(ncid)
    if (status /= nf90_noerr) call handle_err(status)

    ! Populate file

    ! Close file
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    
  end subroutine write_output

  ! #########################################################################################
  ! Error handling for NetCDF
  ! #########################################################################################
  subroutine handle_err(status)
    use netcdf
    integer, intent(in) :: status
    if (status /= nf90_noerr) then
       print*,trim(nf90_strerror(status))

    end if
  end subroutine handle_err
  
end program driver
