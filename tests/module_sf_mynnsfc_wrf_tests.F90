! Module for MYNN SFC scheme tests
module module_sf_mynnsfc_wrf_tests
  implicit none

  contains
    !=================================================================================================================    

    subroutine init_mynn_sfc_flags_for_test_all_true()
      write(*,*) '--- calling  init_mynn_sfc_flags_for_test_all_true()'
      ! for future use
    end subroutine init_mynn_sfc_flags_for_test_all_true
    !=================================================================================================================    
    subroutine init_input_data_for_test()
      integer :: iostat, line_num
      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20
      
      write(*,*) '--- opening data files ---'
      ! Open input file
      close(unit=input_unit)
      open(unit=input_unit, file='./data/input_lnd.txt', status='old', action='read', iostat=iostat)

      if (iostat /= 0) then
          print *, 'Error opening input file'
          stop
      end if

      ! Open output file
      open(unit=output_unit, file='./data/wrf_output_lnd.txt', status='replace', action='write', iostat=iostat)
      write(output_unit,'(A5, A5, A10, A10, A10, A10, A10, A10, A10, A10, A10)')                &
            'itimestep', 'iter', 'T2', 'Q2', 'TH2', 'U10', 'V10', 'HFX', 'LH', 'UST_lnd','PBLH'
      if (iostat /= 0) then
          print *, 'Error opening output file'
          close(output_unit)
          stop
      end if

    end subroutine init_input_data_for_test

    !===============================================================================
    ! Subroutine to process each line
    !===============================================================================
    subroutine process_line(line, out_unit, line_number, flag_iter,  &                         
           U1D,V1D,T1D,QV1D,P1D,dz8w1d,                              &
           U1D2,V1D2,dz2w1d,                                         &
           PSFCPA,PBLH,MAVAIL,XLAND,DX,                              &
           ISFFLX,isftcflx,iz0tlnd,psi_opt,                          &
           compute_flux,compute_diag,                                &
           sigmaf,vegtype,shdmax,ivegsrc,                            & 
           z0pert,ztpert,                                            &
           redrag,sfc_z0_type,                                       &
           itimestep,iter,flag_restart,lsm,lsm_ruc,                  &
           wet,          dry,          icy,                          &
           tskin_wat,    tskin_lnd,    tskin_ice,                    &
           tsurf_wat,    tsurf_lnd,    tsurf_ice,                    &
           qsfc_wat,     qsfc_lnd,     qsfc_ice,                     &
           snowh_wat,    snowh_lnd,    snowh_ice,                    &
           ZNT_wat,      ZNT_lnd,      ZNT_ice,                      &
           UST_wat,      UST_lnd,      UST_ice,                      &
           cm_wat,       cm_lnd,       cm_ice,                       &
           ch_wat,       ch_lnd,       ch_ice,                       &
           rb_wat,       rb_lnd,       rb_ice,                       &
           stress_wat,   stress_lnd,   stress_ice,                   &
           psix_wat,     psix_lnd,     psix_ice,                     &
           psit_wat,     psit_lnd,     psit_ice,                     &
           psix10_wat,   psix10_lnd,   psix10_ice,                   &
           psit2_wat,    psit2_lnd,    psit2_ice,                    &
           HFLX_wat,     HFLX_lnd,     HFLX_ice,                     &
           QFLX_wat,     QFLX_lnd,     QFLX_ice,                     &
           ch,CHS,CHS2,CQS2,CPM,                                     &
           ZNT,USTM,ZOL,MOL,RMOL,                                    &
           PSIM,PSIH,                                                &
           HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                           &
           QGH,QSFC,                                                 &
           U10,V10,TH2,T2,Q2,                                        &
           GZ1OZ0,WSPD,wstar,qstar,                                  &
           spp_sfc,  rstoch1D )         
        
        implicit none
        
        character(len=2000), intent(in) :: line
        integer, intent(in) :: out_unit
        integer, intent(in) :: line_number
        
        ! Variables to store parsed values
        logical, intent(out) :: flag_iter, compute_flux, compute_diag, redrag, flag_restart
        logical, intent(out) :: wet, dry, icy
        real, intent(out) :: U1D, V1D, T1D, QV1D, P1D, dz8w1d, U1D2, V1D2, dz2w1d, &
                PSFCPA, PBLH, MAVAIL, XLAND, DX, sigmaf, shdmax, z0pert,           &
                ztpert
        real, intent(out) :: tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice, &
                qsfc_wat, qsfc_lnd, qsfc_ice, snowh_wat, snowh_lnd, snowh_ice,                 &
                ZNT_wat, ZNT_lnd, ZNT_ice, UST_wat, UST_lnd, UST_ice,                          &
                cm_wat, cm_lnd, cm_ice, ch_wat, ch_lnd, ch_ice,                                &
                rb_wat, rb_lnd, rb_ice, stress_wat, stress_lnd, stress_ice,                    &                   
                psix_wat,     psix_lnd,     psix_ice,                                          &
                psit_wat,     psit_lnd,     psit_ice,                                          &
                psix10_wat,   psix10_lnd,   psix10_ice,                                        &
                psit2_wat,    psit2_lnd,    psit2_ice,                                         &
                HFLX_wat,     HFLX_lnd,     HFLX_ice,                                          &
                QFLX_wat,     QFLX_lnd,     QFLX_ice,                                          &
                ch,CHS,CHS2,CQS2,CPM,                                                          &
                ZNT,USTM,ZOL,MOL,RMOL,                                                         &
                PSIM,PSIH,                                                                     &
                HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                                                &
                QGH,QSFC,                                                                      &
                U10,V10,TH2,T2,Q2,                                                             &
                GZ1OZ0,WSPD,wstar,qstar,                                                       &
                rstoch1D  
        integer :: read_stat           
        integer, intent(out) :: ISFFLX, isftcflx, iz0tlnd, psi_opt, vegtype,                   &
                   ivegsrc, sfc_z0_type, itimestep, iter, lsm, lsm_ruc, spp_sfc
        
        ! Read values from the line with specified format
        read (line, '(L5,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'          // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,' // &
                     'I5,I5,I5,I5,'                                     // &
                     'L5,L5,'                                           // &
                     'E15.3,I5,E15.3,I5,'                               // &
                     'E15.3,E15.3,'                                     // &
                     'L5,I5,'                                           // &
                     'I5,I5,L5,I5,I5,'                                  // &
                     'L5,L5,L5,'                                        // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // & !stress_* values
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,'                               // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,'                                     // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,E15.3,'       // &
                     'E15.3,E15.3,'                                     // &
                     'E15.3,E15.3,E15.3,E15.3,E15.3,'                   // &
                     'E15.3,E15.3,E15.3,E15.3,'                         // &
                     'I5,E15.3)', iostat=read_stat)                        &
                      flag_iter,                                           &
                      U1D,   V1D,   T1D,   QV1D,   P1D,   dz8w1d,          &
                      U1D2,   V1D2,   dz2w1d,                              &
                      PSFCPA,   PBLH,   MAVAIL,   XLAND,   DX,             &
                      ISFFLX,isftcflx,iz0tlnd,psi_opt,                     &
                      compute_flux,compute_diag,                           &
                      sigmaf,   vegtype,   shdmax,   ivegsrc,              &
                      z0pert,   ztpert,                                    &
                      redrag,sfc_z0_type,                                  &
                      itimestep,iter,flag_restart,lsm,lsm_ruc,             &
                            wet,             dry,             icy,         &
                      tskin_wat,       tskin_lnd,       tskin_ice,         &
                      tsurf_wat,       tsurf_lnd,       tsurf_ice,         &
                       qsfc_wat,        qsfc_lnd,        qsfc_ice,         &
                      snowh_wat,       snowh_lnd,       snowh_ice,         &
                        ZNT_wat,         ZNT_lnd,         ZNT_ice,         &
                        UST_wat,         UST_lnd,         UST_ice,         &
                         cm_wat,          cm_lnd,          cm_ice,         &
                         ch_wat,          ch_lnd,          ch_ice,         &
                         rb_wat,          rb_lnd,          rb_ice,         &
                     stress_wat,      stress_lnd,      stress_ice,         &
                     psix_wat,     psix_lnd,     psix_ice,                 & 
                     psit_wat,     psit_lnd,     psit_ice,                 & 
                     psix10_wat,   psix10_lnd,   psix10_ice,               &
                     psit2_wat,    psit2_lnd,    psit2_ice,                &
                     HFLX_wat,     HFLX_lnd,     HFLX_ice,                 &
                     QFLX_wat,     QFLX_lnd,     QFLX_ice,                 &
                     ch,CHS,CHS2,CQS2,CPM,                                 &
                     ZNT,USTM,ZOL,MOL,RMOL,                                &
                     PSIM,PSIH,                                            &
                     HFLX,HFX,QFLX,QFX,LH,FLHC,FLQC,                       &
                     QGH,QSFC,                                             &
                     U10,V10,TH2,T2,Q2,                                    &
                     GZ1OZ0,WSPD,wstar,qstar,                              &
                     spp_sfc,rstoch1D                                              
              
      ! Debug: print the line being read
      write(0,*) "Line length:", len_trim(line)
      write(0,*) "U1D=",U1D, "flag_iter=",flag_iter, "V1D=", V1D," ZNT", ZNT_lnd, "WSPD=",WSPD, "UST=",UST_lnd
        
    end subroutine process_line

    !=================================================================================================================
 
    subroutine wrf_test()

      use module_sf_mynnsfc_land, only : mynnsfc_land

      integer :: iostat, line_num
      integer, parameter :: i=1, j=1

      character(len=2000) :: input_line
      integer, parameter :: input_unit = 10
      integer, parameter :: output_unit = 20

      logical :: flag_iter 
      logical :: compute_flux, compute_diag, redrag, flag_restart
      real ::  u_1, v_1, t_1, qv_1, p_1, dz8w_1, u_2, v_2, dz8w_2,                  &
              psfcpa, pblh, mavial, xland, dx, sigmaf, shdmax, z0pert,              &
              ztpert
      integer :: spp_sfc
      real :: tskin, tsurf, qsfc, snowh, znt, ust, cm, ch, rb, stress, psix, psit,  &
              psix10, psit2, chs, chs2, cqs2, cpm, ustm, zol, mol, rmol,            &
              psim, psih, hfx, qfx, lh, flhc, flqc, u10, v10, th2, t2, q2,          &
              gz1oz0, wspd, wstar, qstar, rstoch_1
      integer :: vegtype, read_stat, isfflx, psi_opt, ivegsrc, sfc_z0_type,         &
              itimestep, iter, lsm, lsm_ruc, isftcflx, iz0tlnd


      ! Additional variables to read process_line
      logical :: wet, dry, icy
      real :: ch_ice, ch_lnd, ch_wat, cm_ice, cm_lnd, cm_wat, dz2w1d, dz8w1d, hflx, hflx_ice, &
              hflx_lnd, hflx_wat, mavail, p1d, psit2_ice, psit2_lnd,                          &
              psit2_wat, psit_ice, psit_lnd, psit_wat, psix_ice, psix_lnd, psix_wat, qflx,    &
              qflx_ice, qflx_lnd, qflx_wat, qgh, rb_wat, rb_lnd, rb_ice,                      &
              stress_wat, stress_lnd, stress_ice, psix10_wat, psix10_lnd, psix10_ice,         &
              ZNT_wat, ZNT_lnd, ZNT_ice, UST_wat, UST_lnd, UST_ice,                           &
              tskin_wat, tskin_lnd, tskin_ice, tsurf_wat, tsurf_lnd, tsurf_ice,               &
              qsfc_wat, qsfc_lnd, qsfc_ice, snowh_wat, snowh_lnd, snowh_ice,                  &
              U1D, V1D, T1D, QV1D,U1D2, V1D2,rstoch1D

      !-----------------------------
      ! input stability function tables
      !-----------------------------
      real,dimension(0:1000) :: psim_stab,psim_unstab, &
                                                      psih_stab,psih_unstab
      ! Variables for error handling
      character(len=512) :: errmsg
      integer :: errflg
      ! Initialize error variables
      errmsg = ''
      errflg = 0

      write(*,*) '--- entering wrf_test subroutine ---'    
      ! Initialize input data for tests
      call init_input_data_for_test()

      ! Read header
      read(input_unit, '(A)', iostat=iostat) input_line
      ! read(input_unit, '(A)', iostat=iostat) input_line
      ! Process each line
      line_num = 0
      do
          read(input_unit, '(A)', iostat=iostat) input_line
          write(0,*) input_line
          
          ! Check for end of file or error
          if (iostat < 0) exit  ! End of file
          if (iostat > 0) then
              print *, 'Error reading line', line_num + 1
              exit
          end if
          
          line_num = line_num + 1
          
          ! Call subroutine to process the line
          call process_line(input_line, output_unit, line_num, flag_iter,     &
           U1D,  V1D,  T1D,  QV1D,  P1D,  dz8w1d,                             &
           U1D2,  V1D2,  dz2w1d,                                              &
           PSFCPA,  PBLH,  MAVAIL,  XLAND,  DX,                               &
           ISFFLX,  isftcflx,  iz0tlnd,  psi_opt,                             &
           compute_flux,  compute_diag,                                       &
           sigmaf,  vegtype,  shdmax,  ivegsrc,                               &
           z0pert,  ztpert,                                                   &
           redrag,  sfc_z0_type,                                              &
           itimestep,  iter,  flag_restart,  lsm,  lsm_ruc,                   &
                  wet,            dry,            icy,                        &
            tskin_wat,      tskin_lnd,      tskin_ice,                        &
            tsurf_wat,      tsurf_lnd,      tsurf_ice,                        &
             qsfc_wat,       qsfc_lnd,       qsfc_ice,                        &
            snowh_wat,      snowh_lnd,      snowh_ice,                        &
              ZNT_wat,        ZNT_lnd,        ZNT_ice,                        &
              UST_wat,        UST_lnd,        UST_ice,                        &
               cm_wat,         cm_lnd,         cm_ice,                        &
               ch_wat,         ch_lnd,         ch_ice,                        &
               rb_wat,         rb_lnd,         rb_ice,                        &
           stress_wat,     stress_lnd,     stress_ice,                        &
             psix_wat,       psix_lnd,       psix_ice,                        &
             psit_wat,       psit_lnd,       psit_ice,                        &
             psix10_wat,     psix10_lnd,     psix10_ice,                      &
             psit2_wat,      psit2_lnd,      psit2_ice,                       &
             HFLX_wat,       HFLX_lnd,       HFLX_ice,                        &
             QFLX_wat,       QFLX_lnd,       QFLX_ice,                        &
             ch,  CHS,  CHS2,  CQS2,  CPM,                                    &
             ZNT,  USTM,  ZOL,  MOL,  RMOL,                                   &
             PSIM,  PSIH,                                                     &
             HFLX,  HFX,  QFLX,  QFX,  LH,  FLHC,  FLQC,                      &
             QGH,  QSFC,                                                      &
             U10,  V10,  TH2,  T2,  Q2,                                       &
             GZ1OZ0,  WSPD,  wstar,  qstar,                                   &
             spp_sfc,  rstoch1D)

         ! Initialize MYNN SFC
          write(*,*) '--- calling  mynnsfc_land() ---'


          call mynnsfc_land (flag_iter=flag_iter, itimestep=itimestep, i=i, j=j,      &
               dx=DX, xland=XLAND,                                                    &
              !3d input - transformed to single point
               u_1=U1D, v_1=V1D, t_1=T1D, qv_1=QV1D,                                  &
               p_1=P1D, dz8w_1=dz8w1d, rho_1=CQS2, u_2=U1D2, & !rho_1
               v_2=V1D2, dz8w_2=dz2w1d,                                               &
               !GFS-related input
               sigmaf=sigmaf, vegtype= vegtype, shdmax=shdmax, ivegsrc=ivegsrc,       &  !intent(in)
               z0pert=z0pert, ztpert=ztpert, redrag=redrag, sfc_z0_type=sfc_z0_type , &  !intent(in)
               !2d variables - transformed to single point
               pblh=pblh, znt=ZNT_lnd , psfcpa=psfcpa, mavail=mavail,                 &  !intent(in)
               tskin=tskin_lnd, tsurf=tsurf_lnd, snowh=snowh_lnd,                     &  !intent(in)
               chs=ch, chs2=CHS2, cqs2=CQS2, cqs=CQS2,                                &   !cqs not being used
               ust=UST_lnd, ustm=USTM, stress=stress_lnd, qsfc=qsfc_lnd,              &  !intent(inout) 
               rmol=RMOL, zol=ZOL, mol=MOL ,                                          &
               psim=PSIM , psih=PSIH, hfx=HFX, qfx=QFX,                               &
               u10=U10 , v10=V10, th2=TH2,                                            &
               t2=T2, q2=Q2, flhc=FLHC, flqc=FLQC ,                                   &
               lh=LH, gz1oz0=GZ1OZ0, wspd=WSPD , rb=rb_lnd,                           &
               cpm=CPM , ch=ch, cm=cm_lnd, rstoch_1=rstoch1D,                         &
               wstar=wstar,qstar=qstar,                                               &
               ck=CQS2, cka=CQS2, cd=CQS2, cda=CQS2 ,                                 & ! ck, cka, cd, cda not being used
               psix=psix_lnd, psit=psit_lnd, psix10=psix10_lnd, psit2=psit2_lnd,      & !fm,fh,fm10,fh2: intent(inout)
               !namelist configuration options
               spp_sfc=spp_sfc, sf_mynn_sfcflux_land=iz0tlnd, ISFFLX=ISFFLX,          &
               flag_restart=flag_restart,flag_cycle=.false., psi_opt=psi_opt,         &
               compute_flux=compute_flux,compute_diag=compute_diag,                   &
               iter=iter, lsm=lsm, lsm_ruc=lsm_ruc,                                   &
               !stability function tables
                psim_stab= psim_stab ,psim_unstab=psim_unstab,psih_stab=psih_stab ,psih_unstab=psih_unstab, &
               errmsg= errmsg, errflg=errflg                                    )
 
         write(0,*) "T2=",t2,'chs=',ch,'ust=',UST_lnd,'hfx=',hfx,'wstar=',wstar
         write(0,*) "Read status:", read_stat
         open(20, file = './data/wrf_output_lnd.txt')
         write(20,'(I5, I5, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)')  &
              itimestep,iter,t2,q2,th2,u10,v10,hfx,lh,UST_lnd,pblh

      end do
    end subroutine wrf_test

end module module_sf_mynnsfc_wrf_tests           
