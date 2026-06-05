!
MODULE cpl_mod
!
! This module excahanges and transforms all fields that should be recieved or sent to 
! coupler. The transformation of the fields are done from the grid 1D-array in phylmd 
! to the regular 2D grid accepted by the coupler. Cumulation of the fields for each 
! timestep is done in here. 
!
! Each type of surface that recevie fields from the coupler have a subroutine named 
! cpl_receive_XXX_fields and each surface that have fields to be sent to the coupler 
! have a subroutine named cpl_send_XXX_fields.
!
!*************************************************************************************

! Use statements
!*************************************************************************************
  USE dimphy, ONLY : klon
  USE mod_phys_lmdz_para
  USE ioipsl
  USE iophy

! The module oasis is always used. Without the cpp key CPP_COUPLE only the parameters 
! in the module are compiled and not the subroutines.
  USE oasis
  USE write_field_phy
  USE time_phylmdz_mod, ONLY: day_step_phy
  
! Global attributes
!*************************************************************************************
  IMPLICIT NONE
  PRIVATE

  ! All subroutine are public except cpl_send_all
  PUBLIC :: cpl_init, cpl_receive_frac, cpl_receive_ocean_fields, cpl_receive_seaice_fields, &
       cpl_send_ocean_fields, cpl_send_seaice_fields, cpl_send_land_fields, &
       cpl_send_landice_fields, gath2cpl, cpl_inca
  

! Declaration of module variables
!*************************************************************************************
! variable for coupling period
  INTEGER, SAVE :: nexca
  !$OMP THREADPRIVATE(nexca)

! variables for cumulating fields during a coupling periode :
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_sols, cpl_nsol, cpl_rain
  !$OMP THREADPRIVATE(cpl_sols,cpl_nsol,cpl_rain)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_snow, cpl_evap, cpl_tsol
  !$OMP THREADPRIVATE(cpl_snow,cpl_evap,cpl_tsol)

  REAL, ALLOCATABLE, SAVE:: cpl_delta_sst(:), cpl_delta_sal(:), cpl_dter(:), cpl_dser(:), cpl_dt_ds(:)
  !$OMP THREADPRIVATE(cpl_delta_sst, cpl_delta_sal, cpl_dter, cpl_dser)
  !$OMP THREADPRIVATE(cpl_dt_ds)
  
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_fder, cpl_albe, cpl_taux, cpl_tauy
  !$OMP THREADPRIVATE(cpl_fder,cpl_albe,cpl_taux,cpl_tauy)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_windsp
  !$OMP THREADPRIVATE(cpl_windsp)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_sens_rain, cpl_sens_snow
  !$OMP THREADPRIVATE(cpl_sens_rain, cpl_sens_snow)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_taumod
  !$OMP THREADPRIVATE(cpl_taumod)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_atm_co2
  !$OMP THREADPRIVATE(cpl_atm_co2)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_rriv2D, cpl_rcoa2D, cpl_rlic2D
  !$OMP THREADPRIVATE(cpl_rriv2D,cpl_rcoa2D,cpl_rlic2D)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: rlic_in_frac2D  ! fraction for continental ice
  !$OMP THREADPRIVATE(rlic_in_frac2D)

! variables read from coupler :
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_sst     ! sea surface temperature
  !$OMP THREADPRIVATE(read_sst)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_sit     ! sea ice temperature
  !$OMP THREADPRIVATE(read_sit)

  REAL, ALLOCATABLE, SAVE:: read_sss(:, :)
  ! bulk salinity of the surface layer of the ocean, in ppt
  !$OMP THREADPRIVATE(read_sss)

  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_sic     ! sea ice fraction
  !$OMP THREADPRIVATE(read_sic)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_alb_sic ! albedo at sea ice
  !$OMP THREADPRIVATE(read_alb_sic)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_u0, read_v0 ! ocean surface current
  !$OMP THREADPRIVATE(read_u0,read_v0)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: read_co2     ! ocean co2 flux 
  !$OMP THREADPRIVATE(read_co2)
  INTEGER, ALLOCATABLE, DIMENSION(:), SAVE  :: unity
  !$OMP THREADPRIVATE(unity)
  INTEGER, SAVE                             :: nidct, nidcs
  !$OMP THREADPRIVATE(nidct,nidcs)

! variables to be sent to the coupler
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_sols2D, cpl_nsol2D, cpl_rain2D
  !$OMP THREADPRIVATE(cpl_sols2D, cpl_nsol2D, cpl_rain2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_snow2D, cpl_evap2D, cpl_tsol2D
  !$OMP THREADPRIVATE(cpl_snow2D, cpl_evap2D, cpl_tsol2D)

  REAL, ALLOCATABLE, SAVE:: cpl_delta_sst_2D(:, :), cpl_delta_sal_2D(:, :) 
  REAL, ALLOCATABLE, SAVE:: cpl_dter_2D(:, :), cpl_dser_2D(:, :), cpl_dt_ds_2D(:, :)
  !$OMP THREADPRIVATE(cpl_delta_sst_2D, cpl_delta_sal_2D)
  !$OMP THREADPRIVATE(cpl_dter_2D, cpl_dser_2D, cpl_dt_ds_2D)

  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_fder2D, cpl_albe2D
  !$OMP THREADPRIVATE(cpl_fder2D, cpl_albe2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_taux2D, cpl_tauy2D
  !$OMP THREADPRIVATE(cpl_taux2D, cpl_tauy2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE :: cpl_taumod2D
  !$OMP THREADPRIVATE(cpl_taumod2D)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_windsp2D
  !$OMP THREADPRIVATE(cpl_windsp2D)
  REAL, ALLOCATABLE, DIMENSION(:,:,:), SAVE   :: cpl_sens_rain2D, cpl_sens_snow2D
  !$OMP THREADPRIVATE(cpl_sens_rain2D, cpl_sens_snow2D)
  REAL, ALLOCATABLE, DIMENSION(:,:), SAVE   :: cpl_atm_co22D
  !$OMP THREADPRIVATE(cpl_atm_co22D)

!!!!!!!!!! variable for calving
  INTEGER, PARAMETER :: nb_zone_calving = 3
  REAL,ALLOCATABLE, DIMENSION(:,:,:),SAVE :: area_calving 
  !$OMP THREADPRIVATE(area_calving)
  REAL,ALLOCATABLE, DIMENSION(:,:),SAVE :: cell_area2D 
  !$OMP THREADPRIVATE(cell_area2D)
  INTEGER, SAVE :: ind_calving(nb_zone_calving) 
  !$OMP THREADPRIVATE(ind_calving)

  LOGICAL,SAVE :: cpl_old_calving
  !$OMP THREADPRIVATE(cpl_old_calving)

  ! Id for fields sent to ocean
  INTEGER, PARAMETER :: ids_tauxxu = 1
  INTEGER, PARAMETER :: ids_tauyyu = 2
  INTEGER, PARAMETER :: ids_tauzzu = 3
  INTEGER, PARAMETER :: ids_tauxxv = 4
  INTEGER, PARAMETER :: ids_tauyyv = 5
  INTEGER, PARAMETER :: ids_tauzzv = 6
  INTEGER, PARAMETER :: ids_windsp = 7
  INTEGER, PARAMETER :: ids_shfice = 8
  INTEGER, PARAMETER :: ids_shfoce = 9
  INTEGER, PARAMETER :: ids_shftot = 10
  INTEGER, PARAMETER :: ids_nsfice = 11
  INTEGER, PARAMETER :: ids_nsfoce = 12
  INTEGER, PARAMETER :: ids_nsftot = 13
  INTEGER, PARAMETER :: ids_dflxdt = 14
  INTEGER, PARAMETER :: ids_totrai = 15
  INTEGER, PARAMETER :: ids_totsno = 16
  INTEGER, PARAMETER :: ids_toteva = 17
  INTEGER, PARAMETER :: ids_icevap = 18
  INTEGER, PARAMETER :: ids_ocevap = 19
  INTEGER, PARAMETER :: ids_calvin = 20
  INTEGER, PARAMETER :: ids_liqrun = 21
  INTEGER, PARAMETER :: ids_runcoa = 22
  INTEGER, PARAMETER :: ids_rivflu = 23
  INTEGER, PARAMETER :: ids_atmco2 = 24
  INTEGER, PARAMETER :: ids_taumod = 25
  INTEGER, PARAMETER :: ids_qraioc = 26
  INTEGER, PARAMETER :: ids_qsnooc = 27
  INTEGER, PARAMETER :: ids_qraiic = 28
  INTEGER, PARAMETER :: ids_qsnoic = 29
  INTEGER, PARAMETER :: ids_delta_sst = 30, ids_delta_sal = 31, ids_dter = 32, &
       ids_dser = 33, ids_dt_ds = 34
  INTEGER, PARAMETER :: ids_atmn2o = 35
  INTEGER, PARAMETER :: ids_atmndp = 36
  INTEGER, PARAMETER :: ids_atmnh3 = 37

  INTEGER, PARAMETER :: maxsend    = 37  ! Maximum number of fields to send
  INTEGER, PARAMETER :: maxsend_phys = 34 ! Maximum number of fields to send in LMDZ phys - the last one will be send by Inca

  ! Id for fields received from ocean

  INTEGER, PARAMETER :: idr_sisutw = 1
  INTEGER, PARAMETER :: idr_icecov = 2
  INTEGER, PARAMETER :: idr_icealw = 3
  INTEGER, PARAMETER :: idr_icetem = 4
  INTEGER, PARAMETER :: idr_curenx = 5
  INTEGER, PARAMETER :: idr_cureny = 6
  INTEGER, PARAMETER :: idr_curenz = 7
  INTEGER, PARAMETER :: idr_oceco2 = 8
  ! bulk salinity of the surface layer of the ocean, in ppt
  INTEGER, PARAMETER :: idr_sss = 9
  INTEGER, PARAMETER :: idr_ocedms = 10
  INTEGER, PARAMETER :: idr_ocen2o = 11
  INTEGER, PARAMETER :: idr_ocenh3 = 12

  INTEGER, PARAMETER :: maxrecv      = 12     ! Maximum number of fields to receive
  INTEGER, PARAMETER :: maxrecv_phys = 9      ! Maximum number of fields to receive in physiq (without fields received in INCA model )
                                              ! will be changed in next version - INCA fields will be received in LMDZ (like for ORCHIDEE fields)
                                              ! and then send by routine in INCA model

  LOGICAL, SAVE :: cpl_current
!$OMP THREADPRIVATE(cpl_current)

CONTAINS
!
!************************************************************************************
!
  SUBROUTINE cpl_init(dtime, rlon, rlat)
    USE IOIPSL
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, fco2_ocn_day
    USE surface_data
    USE indice_sol_mod
    USE chemistry_cycle_mod, ONLY : dms_cycle_cpl, n2o_cycle_cpl, ndp_cycle_cpl, nh3_cycle_cpl
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat, grid1dTo2d_glo, klon_glo, grid_type, unstructured, regular_lonlat
    USE time_phylmdz_mod, ONLY: annee_ref, day_ini, itau_phy, itaufin_phy
    USE print_control_mod, ONLY: lunout
    USE geometry_mod, ONLY : longitude_deg, latitude_deg, ind_cell_glo, cell_area
    USE ioipsl_getin_p_mod, ONLY: getin_p
    use config_ocean_skin_m, only: activate_ocean_skin

! Input arguments
!*************************************************************************************
    REAL, INTENT(IN)                  :: dtime
    REAL, DIMENSION(klon), INTENT(IN) :: rlon, rlat

! Local variables
!*************************************************************************************
    INTEGER                           :: error, sum_error, ig, i
    INTEGER                           :: jf, nhoridct
    INTEGER                           :: nhoridcs
    INTEGER                           :: idtime
    INTEGER                           :: idayref
    INTEGER                           :: npas ! only for OASIS2
    REAL                              :: zjulian
    REAL, DIMENSION(nbp_lon,nbp_lat)  :: zx_lon, zx_lat
    CHARACTER(len = 20)               :: modname = 'cpl_init'
    CHARACTER(len = 80)               :: abort_message
    CHARACTER(len=80)                 :: clintocplnam, clfromcplnam
    REAL, DIMENSION(klon_mpi)         :: rlon_mpi, rlat_mpi, cell_area_mpi
    INTEGER, DIMENSION(klon_mpi)           :: ind_cell_glo_mpi
    REAL, DIMENSION(nbp_lon,jj_nb)         :: lon2D, lat2D
    INTEGER :: mask_calving(nbp_lon,jj_nb,nb_zone_calving)
    REAL :: pos
    !! WARNING: cpl_current_omp should NOT be put in a THREADPRIVATE statement, it is shared between tasks
    LOGICAL, SAVE                      :: cpl_current_omp

!***************************************
! Use old calving or not (default new calving method)
! New calving method should be used with DYNAMICO and when using new coupling
! weights.
    cpl_old_calving=.FALSE.
    CALL getin_p("cpl_old_calving",cpl_old_calving)
    WRITE(lunout,*)' cpl_old_calving = ', cpl_old_calving


!*************************************************************************************
! Calculate coupling period
!
!*************************************************************************************
     
    npas = itaufin_phy
!    nexca = 86400 / dtime
    nexca = t_coupl / dtime
    WRITE(lunout,*)' ##### Ocean couple #####'
    WRITE(lunout,*)' Valeurs des pas de temps'
    WRITE(lunout,*)' npas = ', npas
    WRITE(lunout,*)' nexca = ', nexca
    
!*************************************************************************************
! Allocate variables
!
!*************************************************************************************
    error = 0
    sum_error = 0

    ALLOCATE(unity(klon), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_sols(klon,2), stat = error) 
    sum_error = sum_error + error
    ALLOCATE(cpl_nsol(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rain(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_snow(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_evap(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_tsol(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_fder(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_albe(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_taux(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_tauy(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_windsp(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_taumod(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_sens_rain(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_sens_snow(klon,2), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rriv2D(nbp_lon,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rcoa2D(nbp_lon,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(cpl_rlic2D(nbp_lon,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(rlic_in_frac2D(nbp_lon,jj_nb), stat=error)
    sum_error = sum_error + error
    ALLOCATE(read_sst(nbp_lon, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_sic(nbp_lon, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_sit(nbp_lon, jj_nb), stat = error)
    sum_error = sum_error + error

    if (activate_ocean_skin >= 1) then
       ALLOCATE(read_sss(nbp_lon, jj_nb), stat = error)
       sum_error = sum_error + error
    
       if (activate_ocean_skin == 2) then
          ALLOCATE(cpl_delta_sst(klon), cpl_delta_sal(klon), cpl_dter(klon), &
               cpl_dser(klon), cpl_dt_ds(klon), stat = error)
          sum_error = sum_error + error
       end if
    end if

    ALLOCATE(read_alb_sic(nbp_lon, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_u0(nbp_lon, jj_nb), stat = error)
    sum_error = sum_error + error
    ALLOCATE(read_v0(nbp_lon, jj_nb), stat = error)
    sum_error = sum_error + error

    IF (carbon_cycle_cpl) THEN
       ALLOCATE(read_co2(nbp_lon, jj_nb), stat = error)
       sum_error = sum_error + error
       ALLOCATE(cpl_atm_co2(klon,2), stat = error)
       sum_error = sum_error + error

! Allocate variable in carbon_cycle_mod
       IF (.NOT.ALLOCATED(fco2_ocn_day)) ALLOCATE(fco2_ocn_day(klon), stat = error)
       sum_error = sum_error + error
    ENDIF

! calving initialization
    ALLOCATE(area_calving(nbp_lon, jj_nb, nb_zone_calving), stat = error)
    sum_error = sum_error + error
    ALLOCATE(cell_area2D(nbp_lon, jj_nb), stat = error)    
    sum_error = sum_error + error

    CALL gather_omp(longitude_deg,rlon_mpi)
    CALL gather_omp(latitude_deg,rlat_mpi)
    CALL gather_omp(ind_cell_glo,ind_cell_glo_mpi)
    CALL gather_omp(cell_area,cell_area_mpi)
      
    IF (is_omp_master) THEN
      CALL Grid1DTo2D_mpi(rlon_mpi,lon2D)
      CALL Grid1DTo2D_mpi(rlat_mpi,lat2D)
      CALL Grid1DTo2D_mpi(cell_area_mpi,cell_area2D)
      !--the next line is required for lat-lon grid and should have no impact
      !--for an unstructured grid for which nbp_lon=1
      !--if north pole in process mpi then divide cell area of pole cell by number of replicates
      IF (is_north_pole_dyn) cell_area2D(:,1)=cell_area2D(:,1)/FLOAT(nbp_lon)
      !--if south pole in process mpi then divide cell area of pole cell by number of replicates
      IF (is_south_pole_dyn) cell_area2D(:,jj_nb)=cell_area2D(:,jj_nb)/FLOAT(nbp_lon)
      mask_calving(:,:,:) = 0 
      WHERE ( lat2D >= 40) mask_calving(:,:,1) = 1
      WHERE ( lat2D < 40 .AND. lat2D > -50) mask_calving(:,:,2) = 1
      WHERE ( lat2D <= -50) mask_calving(:,:,3) = 1
    
    
      DO i=1,nb_zone_calving
        area_calving(:,:,i)=mask_calving(:,:,i)*cell_area2D(:,:)
        pos=1
        IF (i>1) pos = 1 + ((nbp_lon*nbp_lat-1)*(i-1))/(nb_zone_calving-1)
      
        ind_calving(i)=0
        IF (grid_type==unstructured) THEN

          DO ig=1,klon_mpi
            IF (ind_cell_glo_mpi(ig)==pos) ind_calving(i)=ig
          ENDDO

        ELSE IF (grid_type==regular_lonlat) THEN
          IF ((ij_begin<=pos .AND. ij_end>=pos) .OR. (ij_begin<=pos .AND. is_south_pole_dyn )) THEN
            ind_calving(i)=pos-(jj_begin-1)*nbp_lon
          ENDIF
        ENDIF
     
      ENDDO
    ENDIF
    
    IF (sum_error /= 0) THEN
       abort_message='Pb allocation variables couplees'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
!*************************************************************************************
! Initialize the allocated varaibles
!
!*************************************************************************************
    DO ig = 1, klon
       unity(ig) = ig
    ENDDO

!*************************************************************************************
! Initialize coupling
!
!*************************************************************************************
    idtime = INT(dtime)
#ifdef CPP_COUPLE
    WRITE(lunout,*) ' '
    WRITE(lunout,*) ' '
    WRITE(lunout,*) ' INIT OCEAN COUPLING'
    WRITE(lunout,*) ' *******************'
    WRITE(lunout,*) ' '
    WRITE(lunout,*) ' '

    ! Define if coupling ocean currents or not
!$OMP MASTER
    cpl_current_omp = .FALSE.
    CALL getin('cpl_current', cpl_current_omp)
!$OMP END MASTER
!$OMP BARRIER
    cpl_current = cpl_current_omp
    WRITE(lunout,*) 'Couple ocean currents, cpl_current = ',cpl_current

    ! Define coupling variables
    ! Atmospheric variables to send

!$OMP MASTER
    ALLOCATE(infosend(midcpl)%fld(maxsend), inforecv(midcpl)%fld(maxrecv) )

    infosend(midcpl)%fld(:)%action = .FALSE.
    infosend(midcpl)%fld(:)%nlvl = 1

    infosend(midcpl)%fld(ids_tauxxu)%action = .TRUE. ; infosend(midcpl)%fld(ids_tauxxu)%name = 'COTAUXXU'
    infosend(midcpl)%fld(ids_tauyyu)%action = .TRUE. ; infosend(midcpl)%fld(ids_tauyyu)%name = 'COTAUYYU'
    infosend(midcpl)%fld(ids_tauzzu)%action = .TRUE. ; infosend(midcpl)%fld(ids_tauzzu)%name = 'COTAUZZU'
    infosend(midcpl)%fld(ids_tauxxv)%action = .TRUE. ; infosend(midcpl)%fld(ids_tauxxv)%name = 'COTAUXXV'
    infosend(midcpl)%fld(ids_tauyyv)%action = .TRUE. ; infosend(midcpl)%fld(ids_tauyyv)%name = 'COTAUYYV'
    infosend(midcpl)%fld(ids_tauzzv)%action = .TRUE. ; infosend(midcpl)%fld(ids_tauzzv)%name = 'COTAUZZV'
    infosend(midcpl)%fld(ids_windsp)%action = .TRUE. ; infosend(midcpl)%fld(ids_windsp)%name = 'COWINDSP'
    infosend(midcpl)%fld(ids_shfice)%action = .TRUE. ; infosend(midcpl)%fld(ids_shfice)%name = 'COSHFICE'
    infosend(midcpl)%fld(ids_nsfice)%action = .TRUE. ; infosend(midcpl)%fld(ids_nsfice)%name = 'CONSFICE'
    infosend(midcpl)%fld(ids_dflxdt)%action = .TRUE. ; infosend(midcpl)%fld(ids_dflxdt)%name = 'CODFLXDT'
    infosend(midcpl)%fld(ids_calvin)%action = .TRUE. ; infosend(midcpl)%fld(ids_calvin)%name = 'COCALVIN'

    if (activate_ocean_skin == 2) then
       infosend(midcpl)%fld(ids_delta_sst)%action = .TRUE.
       infosend(midcpl)%fld(ids_delta_sst)%name = 'CODELSST'
       infosend(midcpl)%fld(ids_delta_sal)%action = .TRUE.
       infosend(midcpl)%fld(ids_delta_sal)%name = 'CODELSSS'
       infosend(midcpl)%fld(ids_dter)%action = .TRUE.
       infosend(midcpl)%fld(ids_dter)%name = 'CODELTER'
       infosend(midcpl)%fld(ids_dser)%action = .TRUE.
       infosend(midcpl)%fld(ids_dser)%name = 'CODELSER'
       infosend(midcpl)%fld(ids_dt_ds)%action = .TRUE.
       infosend(midcpl)%fld(ids_dt_ds)%name = 'CODTDS'
    end if

    IF (version_ocean=='nemo') THEN
        infosend(midcpl)%fld(ids_shftot)%action = .TRUE. ; infosend(midcpl)%fld(ids_shftot)%name = 'COQSRMIX'
        infosend(midcpl)%fld(ids_nsftot)%action = .TRUE. ; infosend(midcpl)%fld(ids_nsftot)%name = 'COQNSMIX'
        infosend(midcpl)%fld(ids_totrai)%action = .TRUE. ; infosend(midcpl)%fld(ids_totrai)%name = 'COTOTRAI'
        infosend(midcpl)%fld(ids_totsno)%action = .TRUE. ; infosend(midcpl)%fld(ids_totsno)%name = 'COTOTSNO'
        infosend(midcpl)%fld(ids_toteva)%action = .TRUE. ; infosend(midcpl)%fld(ids_toteva)%name = 'COTOTEVA'
        infosend(midcpl)%fld(ids_icevap)%action = .TRUE. ; infosend(midcpl)%fld(ids_icevap)%name = 'COICEVAP'
        infosend(midcpl)%fld(ids_liqrun)%action = .TRUE. ; infosend(midcpl)%fld(ids_liqrun)%name = 'COLIQRUN'
        infosend(midcpl)%fld(ids_taumod)%action = .TRUE. ; infosend(midcpl)%fld(ids_taumod)%name = 'COTAUMOD'
        IF (carbon_cycle_cpl) THEN
            infosend(midcpl)%fld(ids_atmco2)%action = .TRUE. ; infosend(midcpl)%fld(ids_atmco2)%name = 'COATMCO2'
        ENDIF
        IF (n2o_cycle_cpl) THEN
            infosend(midcpl)%fld(ids_atmn2o)%action = .TRUE. ; infosend(midcpl)%fld(ids_atmn2o)%name = 'COATMN2O'
        ENDIF
        IF (ndp_cycle_cpl) THEN
            infosend(midcpl)%fld(ids_atmndp)%action = .TRUE. ; infosend(midcpl)%fld(ids_atmndp)%name = 'COATMNDP'
        ENDIF
        IF (nh3_cycle_cpl) THEN
            infosend(midcpl)%fld(ids_atmnh3)%action = .TRUE. ; infosend(midcpl)%fld(ids_atmnh3)%name = 'COATMNH3'
        ENDIF
        infosend(midcpl)%fld(ids_qraioc)%action = .TRUE. ; infosend(midcpl)%fld(ids_qraioc)%name = 'COQRAIOC'
        infosend(midcpl)%fld(ids_qsnooc)%action = .TRUE. ; infosend(midcpl)%fld(ids_qsnooc)%name = 'COQSNOOC'
        infosend(midcpl)%fld(ids_qraiic)%action = .TRUE. ; infosend(midcpl)%fld(ids_qraiic)%name = 'COQRAIIC'
        infosend(midcpl)%fld(ids_qsnoic)%action = .TRUE. ; infosend(midcpl)%fld(ids_qsnoic)%name = 'COQSNOIC'

    ELSE IF (version_ocean=='opa8') THEN
        infosend(midcpl)%fld(ids_shfoce)%action = .TRUE. ; infosend(midcpl)%fld(ids_shfoce)%name = 'COSHFOCE'
        infosend(midcpl)%fld(ids_nsfoce)%action = .TRUE. ; infosend(midcpl)%fld(ids_nsfoce)%name = 'CONSFOCE'
        infosend(midcpl)%fld(ids_icevap)%action = .TRUE. ; infosend(midcpl)%fld(ids_icevap)%name = 'COTFSICE'
        infosend(midcpl)%fld(ids_ocevap)%action = .TRUE. ; infosend(midcpl)%fld(ids_ocevap)%name = 'COTFSOCE'
        infosend(midcpl)%fld(ids_totrai)%action = .TRUE. ; infosend(midcpl)%fld(ids_totrai)%name = 'COTOLPSU'
        infosend(midcpl)%fld(ids_totsno)%action = .TRUE. ; infosend(midcpl)%fld(ids_totsno)%name = 'COTOSPSU'
        infosend(midcpl)%fld(ids_runcoa)%action = .TRUE. ; infosend(midcpl)%fld(ids_runcoa)%name = 'CORUNCOA'
        infosend(midcpl)%fld(ids_rivflu)%action = .TRUE. ; infosend(midcpl)%fld(ids_rivflu)%name = 'CORIVFLU'
    ENDIF

    ! Oceanic variables to receive
    inforecv(midcpl)%fld(:)%action = .FALSE.
    inforecv(midcpl)%fld(:)%nlvl = 1

    inforecv(midcpl)%fld(idr_sisutw)%action = .TRUE. ; inforecv(midcpl)%fld(idr_sisutw)%name = 'SISUTESW'
    inforecv(midcpl)%fld(idr_icecov)%action = .TRUE. ; inforecv(midcpl)%fld(idr_icecov)%name = 'SIICECOV'
    inforecv(midcpl)%fld(idr_icealw)%action = .TRUE. ; inforecv(midcpl)%fld(idr_icealw)%name = 'SIICEALW'
    inforecv(midcpl)%fld(idr_icetem)%action = .TRUE. ; inforecv(midcpl)%fld(idr_icetem)%name = 'SIICTEMW'

    if (activate_ocean_skin >= 1) then
       inforecv(midcpl)%fld(idr_sss)%action = .TRUE.
       inforecv(midcpl)%fld(idr_sss)%name = 'SISUSALW'
    end if

    IF (cpl_current ) THEN
       inforecv(midcpl)%fld(idr_curenx)%action = .TRUE. ; inforecv(midcpl)%fld(idr_curenx)%name = 'CURRENTX'
       inforecv(midcpl)%fld(idr_cureny)%action = .TRUE. ; inforecv(midcpl)%fld(idr_cureny)%name = 'CURRENTY'
       inforecv(midcpl)%fld(idr_curenz)%action = .TRUE. ; inforecv(midcpl)%fld(idr_curenz)%name = 'CURRENTZ'
    ENDIF

    IF (carbon_cycle_cpl ) THEN
       inforecv(midcpl)%fld(idr_oceco2)%action = .TRUE. ; inforecv(midcpl)%fld(idr_oceco2)%name = 'SICO2FLX'
    ENDIF
    IF (dms_cycle_cpl) THEN
       inforecv(midcpl)%fld(idr_ocedms)%action = .TRUE. ; inforecv(midcpl)%fld(idr_ocedms)%name = 'SIDMSFLX'
    ENDIF
    IF (n2o_cycle_cpl) THEN
       inforecv(midcpl)%fld(idr_ocen2o)%action = .TRUE. ; inforecv(midcpl)%fld(idr_ocen2o)%name = 'SIN2OFLX'
    ENDIF
    IF (nh3_cycle_cpl) THEN
       inforecv(midcpl)%fld(idr_ocenh3)%action = .TRUE. ; inforecv(midcpl)%fld(idr_ocenh3)%name = 'SINH3FLX'
    ENDIF

    ! Configure coupler
    CALL cpl_vardef(midcpl)

!$OMP END MASTER
#endif

!*************************************************************************************
! initialize NetCDF output
!
!*************************************************************************************
    IF (is_sequential) THEN
       idayref = day_ini
       CALL ymds2ju(annee_ref, 1, idayref, 0.0, zjulian)
       CALL grid1dTo2d_glo(rlon,zx_lon)
       DO i = 1, nbp_lon
          zx_lon(i,1) = rlon(i+1)
          zx_lon(i,nbp_lat) = rlon(i+1)
       ENDDO
       CALL grid1dTo2d_glo(rlat,zx_lat)
       clintocplnam="cpl_atm_tauflx"
       CALL histbeg(clintocplnam,nbp_lon,zx_lon(:,1),nbp_lat,zx_lat(1,:),&
            1,nbp_lon,1,nbp_lat, itau_phy,zjulian,dtime,nhoridct,nidct) 
! no vertical axis
       CALL histdef(nidct, 'tauxe','tauxe', &
            "-",nbp_lon,nbp_lat, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tauyn','tauyn', &
            "-",nbp_lon,nbp_lat, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tmp_lon','tmp_lon', &
            "-",nbp_lon,nbp_lat, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       CALL histdef(nidct, 'tmp_lat','tmp_lat', &
            "-",nbp_lon,nbp_lat, nhoridct, 1, 1, 1, -99, 32, "inst", dtime,dtime)
       DO jf=1,maxsend_phys
         IF (infosend(midcpl)%fld(i)%action) THEN
             CALL histdef(nidct, infosend(midcpl)%fld(i)%name ,infosend(midcpl)%fld(i)%name , &
                "-",nbp_lon,nbp_lat,nhoridct,1,1,1,-99,32,"inst",dtime,dtime)
         ENDIF
       ENDDO
       CALL histend(nidct)
       CALL histsync(nidct)
       
       clfromcplnam="cpl_atm_sst"
       CALL histbeg(clfromcplnam,nbp_lon,zx_lon(:,1),nbp_lat,zx_lat(1,:),1,nbp_lon,1,nbp_lat, &
            0,zjulian,dtime,nhoridcs,nidcs) 
! no vertical axis
       DO jf=1,maxrecv_phys
         IF (inforecv(midcpl)%fld(i)%action) THEN
             CALL histdef(nidcs,inforecv(midcpl)%fld(i)%name ,inforecv(midcpl)%fld(i)%name , &
                "-",nbp_lon,nbp_lat,nhoridcs,1,1,1,-99,32,"inst",dtime,dtime)
         ENDIF
       ENDDO
       CALL histend(nidcs)
       CALL histsync(nidcs)

    ENDIF    ! is_sequential
    
!*************************************************************************************
! compatibility test
!
!*************************************************************************************
    IF (carbon_cycle_cpl .AND. version_ocean=='opa8') THEN
       abort_message='carbon_cycle_cpl does not work with opa8'
       CALL abort_physic(modname,abort_message,1)
    ENDIF

  END SUBROUTINE cpl_init

!
!*************************************************************************************
!

  SUBROUTINE cpl_inca

     USE chemistry_cycle_mod, ONLY : dms_cycle_cpl, n2o_cycle_cpl, ndp_cycle_cpl, nh3_cycle_cpl
     USE lmdz_cppkeys_wrapper, ONLY: CPPKEY_INCA

!$OMP MASTER
     IF (CPPKEY_INCA) THEN
        IF (dms_cycle_cpl .OR. n2o_cycle_cpl .OR. ndp_cycle_cpl .OR. nh3_cycle_cpl) THEN
           CALL init_inca_oasis(inforecv(midcpl)%fld(idr_ocedms:idr_ocenh3),infosend(midcpl)%fld(ids_atmn2o:ids_atmnh3))
        ENDIF
    END IF
!$OMP END MASTER
  END SUBROUTINE cpl_inca

!
!*************************************************************************************
!
 
  SUBROUTINE cpl_receive_frac(itime, dtime, pctsrf, is_modified)
! This subroutine receives from coupler for both ocean and seaice
! 4 fields : read_sst, read_sic, read_sit and read_alb_sic. 
! The new sea-ice-land-landice fraction is returned. The others fields 
! are stored in this module.
    USE yomcst_mod_h
USE surface_data
    USE geometry_mod, ONLY : longitude_deg, latitude_deg
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl
    USE indice_sol_mod
    USE print_control_mod, ONLY: lunout
    USE time_phylmdz_mod, ONLY: start_time, itau_phy
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
    use config_ocean_skin_m, only: activate_ocean_skin



! Arguments
!************************************************************************************
    INTEGER, INTENT(IN)                        :: itime
    REAL, INTENT(IN)                           :: dtime
    REAL, DIMENSION(klon,nbsrf), INTENT(INOUT) :: pctsrf
    LOGICAL, INTENT(OUT)                       :: is_modified

! Local variables
!************************************************************************************
    INTEGER                                 :: j, i, time_sec
    INTEGER                                 :: istart,iend
    INTEGER                                 :: itau_w
    INTEGER, DIMENSION(nbp_lon*nbp_lat)     :: ndexcs
    CHARACTER(len = 20)                     :: modname = 'cpl_receive_frac'
    CHARACTER(len = 80)                     :: abort_message
    REAL, DIMENSION(klon)                   :: read_sic1D
    REAL, DIMENSION(nbp_lon,jj_nb,maxrecv_phys)      :: tab_read_flds
    REAL, DIMENSION(klon,nbsrf)             :: pctsrf_old
    REAL, DIMENSION(klon_mpi)               :: rlon_mpi, rlat_mpi
    REAL, DIMENSION(nbp_lon,jj_nb)             :: tmp_lon, tmp_lat
    REAL, DIMENSION(nbp_lon,jj_nb)             :: tmp_r0
    REAL, DIMENSION(nbp_lon*jj_nb,1)             :: field

!*************************************************************************************
! Start calculation
! Get fields from coupler
!
!*************************************************************************************

    is_modified=.FALSE.

! Check if right moment to receive from coupler
    IF (MOD(itime, nexca) == 1) THEN
       is_modified=.TRUE.
 
       time_sec=(itime-1)*dtime
#ifdef CPP_COUPLE
!$OMP MASTER
! ======================================================================
! L. Fairhead (09/2003) adapted From L.Z.X Li: this section reads the SST
! and Sea-Ice provided by the coupler. Adaptation to psmile library
!======================================================================
    WRITE (lunout,*) ' '
    WRITE (lunout,*) 'Fromcpl: Reading fields from CPL, time_sec=',time_sec
    WRITE (lunout,*) ' '

    istart=ii_begin
    IF (is_south_pole_dyn) THEN
       iend=(jj_end-jj_begin)*nbp_lon+nbp_lon
    ELSE
       iend=(jj_end-jj_begin)*nbp_lon+ii_end
    ENDIF

    DO i = 1, maxrecv_phys
      IF (inforecv(midcpl)%fld(i)%action .AND. inforecv(midcpl)%fld(i)%nid .NE. -1) THEN
          field(:,:) = -99999.
          CALL cpl_rcv( midcpl, i, time_sec, field(istart:iend,:) )
          tab_read_flds(:,:,i) = RESHAPE(field(:,1),(/nbp_lon,jj_nb/))
      ENDIF
    END DO
!$OMP END MASTER
#endif
    
! NetCDF output of received fields
       IF (is_sequential) THEN
          ndexcs(:) = 0
          itau_w = itau_phy + itime + start_time * day_step_phy
          DO i = 1, maxrecv_phys
            IF (inforecv(midcpl)%fld(i)%action) THEN
                CALL histwrite(nidcs,inforecv(midcpl)%fld(i)%name,itau_w,tab_read_flds(:,:,i),nbp_lon*(nbp_lat),ndexcs)
            ENDIF
          ENDDO
       ENDIF


! Save each field in a 2D array. 
!$OMP MASTER
       read_sst(:,:)     = tab_read_flds(:,:,idr_sisutw)  ! Sea surface temperature
       read_sic(:,:)     = tab_read_flds(:,:,idr_icecov)  ! Sea ice concentration
       read_alb_sic(:,:) = tab_read_flds(:,:,idr_icealw)  ! Albedo at sea ice
       read_sit(:,:)     = tab_read_flds(:,:,idr_icetem)  ! Sea ice temperature
       if (activate_ocean_skin >= 1) read_sss(:,:) = tab_read_flds(:,:,idr_sss)
!$OMP END MASTER

       IF (cpl_current) THEN

! Transform the longitudes and latitudes on 2D arrays
          CALL gather_omp(longitude_deg,rlon_mpi)
          CALL gather_omp(latitude_deg,rlat_mpi)
!$OMP MASTER
          CALL Grid1DTo2D_mpi(rlon_mpi,tmp_lon)
          CALL Grid1DTo2D_mpi(rlat_mpi,tmp_lat)

! Transform the currents from cartesian to spheric coordinates
! tmp_r0 should be zero
          CALL geo2atm(nbp_lon, jj_nb, tab_read_flds(:,:,idr_curenx), &
             tab_read_flds(:,:,idr_cureny), tab_read_flds(:,:,idr_curenz), &
               tmp_lon, tmp_lat, &
               read_u0(:,:), read_v0(:,:), tmp_r0(:,:))
!$OMP END MASTER

      ELSE
          read_u0(:,:) = 0.
          read_v0(:,:) = 0.
      ENDIF

       IF (carbon_cycle_cpl) THEN
!$OMP MASTER
           read_co2(:,:) = tab_read_flds(:,:,idr_oceco2) ! CO2 flux
!$OMP END MASTER
       ENDIF

!*************************************************************************************
!  Transform seaice fraction (read_sic : ocean-seaice mask) into global 
!  fraction (pctsrf : ocean-seaice-land-landice mask)
!
!*************************************************************************************
       CALL cpl2gath(read_sic, read_sic1D, klon, unity)

       pctsrf_old(:,:) = pctsrf(:,:)
       DO i = 1, klon
          ! treatment only of points with ocean and/or seaice
          ! old land-ocean mask can not be changed
          IF (pctsrf_old(i,is_oce) + pctsrf_old(i,is_sic) > 0.) THEN
             pctsrf(i,is_sic) = (pctsrf_old(i,is_oce) + pctsrf_old(i,is_sic)) &
                  * read_sic1D(i)
             pctsrf(i,is_oce) = (pctsrf_old(i,is_oce) + pctsrf_old(i,is_sic)) &
                  - pctsrf(i,is_sic)
          ENDIF
       ENDDO

    ENDIF ! if time to receive

  END SUBROUTINE cpl_receive_frac

!
!*************************************************************************************
!

  SUBROUTINE cpl_receive_ocean_fields(knon, knindex, tsurf_new, u0_new, &
       v0_new, sss)
!
! This routine returns the field for the ocean that has been read from the coupler
! (done earlier with cpl_receive_frac). The field is the temperature.
! The temperature is transformed into 1D array with valid points from index 1 to knon.
!
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, fco2_ocn_day
    USE indice_sol_mod
    use config_ocean_skin_m, only: activate_ocean_skin

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex

! Output arguments
!*************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)      :: tsurf_new

    REAL, INTENT(OUT):: sss(:) ! (klon)
    ! bulk salinity of the surface layer of the ocean, in ppt

    REAL, DIMENSION(klon), INTENT(OUT)      :: u0_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: v0_new

! Local variables
!*************************************************************************************
    INTEGER                  :: i
    INTEGER, DIMENSION(klon) :: index
    REAL, DIMENSION(klon)    :: sic_new

!*************************************************************************************
! Transform read_sst into compressed 1D variable tsurf_new
!
!*************************************************************************************
    CALL cpl2gath(read_sst, tsurf_new, knon, knindex)
    if (activate_ocean_skin >= 1) CALL cpl2gath(read_sss, sss, knon, knindex)
    CALL cpl2gath(read_sic, sic_new, knon, knindex)
    CALL cpl2gath(read_u0, u0_new, knon, knindex)
    CALL cpl2gath(read_v0, v0_new, knon, knindex)

!*************************************************************************************
! Transform read_co2 into uncompressed 1D variable fco2_ocn_day added directly in 
! the module carbon_cycle_mod
!
!*************************************************************************************
    IF (carbon_cycle_cpl) THEN
       DO i=1,klon
          index(i)=i
       ENDDO
       CALL cpl2gath(read_co2, fco2_ocn_day, klon, index)
    ENDIF

!*************************************************************************************
! The fields received from the coupler have to be weighted with the fraction of ocean 
! in relation to the total sea-ice+ocean
!
!*************************************************************************************
    DO i=1, knon
       tsurf_new(i) = tsurf_new(i)/(1. - sic_new(i))
    ENDDO

  END SUBROUTINE cpl_receive_ocean_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_receive_seaice_fields(knon, knindex, &
       tsurf_new, alb_new, u0_new, v0_new)
!
! This routine returns the fields for the seaice that have been read from the coupler
! (done earlier with cpl_receive_frac). These fields are the temperature and 
! albedo at sea ice surface and fraction of sea ice.
! The fields are transformed into 1D arrays with valid points from index 1 to knon. 
!

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex

! Output arguments
!*************************************************************************************
    REAL, DIMENSION(klon), INTENT(OUT)      :: tsurf_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: alb_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: u0_new
    REAL, DIMENSION(klon), INTENT(OUT)      :: v0_new

! Local variables
!*************************************************************************************
    INTEGER               :: i
    REAL, DIMENSION(klon) :: sic_new

!*************************************************************************************
! Transform fields read from coupler from 2D into compressed 1D variables
!
!*************************************************************************************
    CALL cpl2gath(read_sit, tsurf_new, knon, knindex)
    CALL cpl2gath(read_alb_sic, alb_new, knon, knindex)
    CALL cpl2gath(read_sic, sic_new, knon, knindex)
    CALL cpl2gath(read_u0, u0_new, knon, knindex)
    CALL cpl2gath(read_v0, v0_new, knon, knindex)

!*************************************************************************************
! The fields received from the coupler have to be weighted with the sea-ice 
! concentration (in relation to the total sea-ice + ocean).
!
!*************************************************************************************
    DO i= 1, knon
       tsurf_new(i) = tsurf_new(i) / sic_new(i)
       alb_new(i)   = alb_new(i)   / sic_new(i)
    ENDDO

  END SUBROUTINE cpl_receive_seaice_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_ocean_fields(itime, knon, knindex, &
       swdown, lwdown, fluxlat, fluxsens, &
       precip_rain, precip_snow, evap, tsurf, fder, albsol, taux, tauy, windsp,&
       sens_prec_liq, sens_prec_sol, lat_prec_liq, lat_prec_sol, delta_sst, &
       delta_sal, dTer, dSer, dt_ds)

    ! This subroutine cumulates some fields for each time-step during
    ! a coupling period. At last time-step in a coupling period the
    ! fields are transformed to the grid accepted by the coupler. No
    ! sending to the coupler will be done from here (it is done in
    ! cpl_send_seaice_fields). Crucial hypothesis is that the surface
    ! fractions do not change between coupling time-steps.

    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl, co2_send
    USE indice_sol_mod
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
    use config_ocean_skin_m, only: activate_ocean_skin

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: itime
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, DIMENSION(klon), INTENT(IN)       :: swdown, lwdown 
    REAL, DIMENSION(klon), INTENT(IN)       :: fluxlat, fluxsens
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: evap, tsurf, fder, albsol
    REAL, DIMENSION(klon), INTENT(IN)       :: taux, tauy, windsp
    REAL, INTENT(IN):: sens_prec_liq(:), sens_prec_sol(:) ! (knon)
    REAL, DIMENSION(klon), INTENT(IN)       :: lat_prec_liq, lat_prec_sol
    
    REAL, intent(in):: delta_sst(:) ! (knon)
    ! Ocean-air interface temperature minus bulk SST, in
    ! K. Defined only if activate_ocean_skin >= 1.

    real, intent(in):: delta_sal(:) ! (knon)
    ! Ocean-air interface salinity minus bulk salinity, in ppt.

    REAL, intent(in):: dter(:) ! (knon)
    ! Temperature variation in the diffusive microlayer, that is
    ! ocean-air interface temperature minus subskin temperature. In
    ! K.

    REAL, intent(in):: dser(:) ! (knon)
    ! Salinity variation in the diffusive microlayer, that is
    ! ocean-air interface salinity minus subskin salinity. In ppt.

    real, intent(in):: dt_ds(:) ! (knon)
    ! (tks / tkt) * dTer, in K

! Local variables
!*************************************************************************************
    INTEGER                                 :: cpl_index, ig 
    INTEGER                                 :: error, sum_error
    CHARACTER(len = 25)                     :: modname = 'cpl_send_ocean_fields'
    CHARACTER(len = 80)                     :: abort_message

!*************************************************************************************
! Start calculation
! The ocean points are saved with second array index=1
!
!*************************************************************************************
    cpl_index = 1

!*************************************************************************************
! Reset fields to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
       cpl_sols(1:knon,cpl_index) = 0.0
       cpl_nsol(1:knon,cpl_index) = 0.0
       cpl_rain(1:knon,cpl_index) = 0.0
       cpl_snow(1:knon,cpl_index) = 0.0
       cpl_evap(1:knon,cpl_index) = 0.0
       cpl_tsol(1:knon,cpl_index) = 0.0
       cpl_fder(1:knon,cpl_index) = 0.0
       cpl_albe(1:knon,cpl_index) = 0.0
       cpl_taux(1:knon,cpl_index) = 0.0
       cpl_tauy(1:knon,cpl_index) = 0.0
       cpl_windsp(1:knon,cpl_index) = 0.0
       cpl_sens_rain(1:knon,cpl_index) = 0.0
       cpl_sens_snow(1:knon,cpl_index) = 0.0
       cpl_taumod(1:knon,cpl_index) = 0.0
       IF (carbon_cycle_cpl) cpl_atm_co2(1:knon,cpl_index) = 0.0

       if (activate_ocean_skin == 2) then
          cpl_delta_sst = 0.
          cpl_delta_sal = 0.
          cpl_dter = 0.
          cpl_dser = 0.
          cpl_dt_ds = 0.
       end if
    ENDIF
       
!*************************************************************************************
! Cumulate at each time-step
!
!*************************************************************************************    
    DO ig = 1, knon
       cpl_sols(ig,cpl_index) = cpl_sols(ig,cpl_index) + &
            swdown(ig)      / REAL(nexca)
       cpl_nsol(ig,cpl_index) = cpl_nsol(ig,cpl_index) + &
            (lwdown(ig) + fluxlat(ig) +fluxsens(ig)) / REAL(nexca)
       cpl_rain(ig,cpl_index) = cpl_rain(ig,cpl_index) + &
            precip_rain(ig) / REAL(nexca)
       cpl_snow(ig,cpl_index) = cpl_snow(ig,cpl_index) + &
            precip_snow(ig) / REAL(nexca)
       cpl_evap(ig,cpl_index) = cpl_evap(ig,cpl_index) + &
            evap(ig)        / REAL(nexca)
       cpl_tsol(ig,cpl_index) = cpl_tsol(ig,cpl_index) + &
            tsurf(ig)       / REAL(nexca)
       cpl_fder(ig,cpl_index) = cpl_fder(ig,cpl_index) + &
            fder(ig)        / REAL(nexca)
       cpl_albe(ig,cpl_index) = cpl_albe(ig,cpl_index) + &
            albsol(ig)      / REAL(nexca)
       cpl_taux(ig,cpl_index) = cpl_taux(ig,cpl_index) + &
            taux(ig)        / REAL(nexca)
       cpl_tauy(ig,cpl_index) = cpl_tauy(ig,cpl_index) + &
            tauy(ig)        / REAL(nexca)      
       cpl_windsp(ig,cpl_index) = cpl_windsp(ig,cpl_index) + &
            windsp(ig)      / REAL(nexca)
       cpl_sens_rain(ig,cpl_index) = cpl_sens_rain(ig,cpl_index) + &
            sens_prec_liq(ig)      / REAL(nexca)
       cpl_sens_snow(ig,cpl_index) = cpl_sens_snow(ig,cpl_index) + &
            sens_prec_sol(ig)      / REAL(nexca)
       cpl_taumod(ig,cpl_index) =   cpl_taumod(ig,cpl_index) + &
          SQRT ( taux(ig)*taux(ig)+tauy(ig)*tauy(ig) ) / REAL (nexca)

       IF (carbon_cycle_cpl) THEN
          cpl_atm_co2(ig,cpl_index) = cpl_atm_co2(ig,cpl_index) + &
               co2_send(knindex(ig))/ REAL(nexca) 
!!---OB: this is correct but why knindex ??
       ENDIF

       if (activate_ocean_skin == 2) then
          cpl_delta_sst(ig) = cpl_delta_sst(ig) + delta_sst(ig) / REAL(nexca)
          cpl_delta_sal(ig) = cpl_delta_sal(ig) + delta_sal(ig) / REAL(nexca)
          cpl_dter(ig) = cpl_dter(ig) + dter(ig) / REAL(nexca)
          cpl_dser(ig) = cpl_dser(ig) + dser(ig) / REAL(nexca)
          cpl_dt_ds(ig) = cpl_dt_ds(ig) + dt_ds(ig) / REAL(nexca)
       end if
     ENDDO

!*************************************************************************************
! If the time-step corresponds to the end of coupling period the 
! fields are transformed to the 2D grid. 
! No sending to the coupler (it is done from cpl_send_seaice_fields).
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 0) THEN

       IF (.NOT. ALLOCATED(cpl_sols2D)) THEN
          sum_error = 0
          ALLOCATE(cpl_sols2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_nsol2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_rain2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_snow2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_evap2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tsol2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_fder2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_albe2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taux2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tauy2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_windsp2D(nbp_lon,jj_nb), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_sens_rain2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_sens_snow2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taumod2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          
          IF (carbon_cycle_cpl) THEN
             ALLOCATE(cpl_atm_co22D(nbp_lon,jj_nb), stat=error)
             sum_error = sum_error + error
          ENDIF

          if (activate_ocean_skin == 2) then
             ALLOCATE(cpl_delta_sst_2D(nbp_lon, jj_nb), &
                  cpl_delta_sal_2D(nbp_lon, jj_nb), &
                  cpl_dter_2D(nbp_lon, jj_nb), cpl_dser_2D(nbp_lon, jj_nb), &
                  cpl_dt_ds_2D(nbp_lon, jj_nb), stat = error)
             sum_error = sum_error + error
          end if

          IF (sum_error /= 0) THEN
             abort_message='Pb allocation variables couplees pour l''ecriture'
             CALL abort_physic(modname,abort_message,1)
          ENDIF
       ENDIF
       

       CALL gath2cpl(cpl_sols(:,cpl_index), cpl_sols2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_nsol(:,cpl_index), cpl_nsol2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_rain(:,cpl_index), cpl_rain2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_snow(:,cpl_index), cpl_snow2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_evap(:,cpl_index), cpl_evap2D(:,:,cpl_index), &
            knon, knindex)

! cpl_tsol2D(:,:,:) not used!
       CALL gath2cpl(cpl_tsol(:,cpl_index), cpl_tsol2D(:,:, cpl_index), &
            knon, knindex)

! cpl_fder2D(:,:,1) not used, only cpl_fder(:,:,2)!
       CALL gath2cpl(cpl_fder(:,cpl_index), cpl_fder2D(:,:,cpl_index), &
            knon, knindex)

! cpl_albe2D(:,:,:) not used!
       CALL gath2cpl(cpl_albe(:,cpl_index), cpl_albe2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taux(:,cpl_index), cpl_taux2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_tauy(:,cpl_index), cpl_tauy2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_windsp(:,cpl_index), cpl_windsp2D(:,:), &
            knon, knindex)

       CALL gath2cpl(cpl_sens_rain(:,cpl_index), cpl_sens_rain2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_sens_snow(:,cpl_index), cpl_sens_snow2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taumod(:,cpl_index), cpl_taumod2D(:,:,cpl_index), &
            knon, knindex)

       IF (carbon_cycle_cpl) &
            CALL gath2cpl(cpl_atm_co2(:,cpl_index), cpl_atm_co22D(:,:), knon, knindex)
       if (activate_ocean_skin == 2) then
          CALL gath2cpl(cpl_delta_sst, cpl_delta_sst_2D, knon, knindex)
          CALL gath2cpl(cpl_delta_sal, cpl_delta_sal_2D, knon, knindex)
          CALL gath2cpl(cpl_dter, cpl_dter_2D, knon, knindex)
          CALL gath2cpl(cpl_dser, cpl_dser_2D, knon, knindex)
          CALL gath2cpl(cpl_dt_ds, cpl_dt_ds_2D, knon, knindex)
       end if
    ENDIF

  END SUBROUTINE cpl_send_ocean_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_seaice_fields(itime, dtime, knon, knindex, &
       pctsrf, lafin, rlon, rlat, &
       swdown, lwdown, fluxlat, fluxsens, &
       precip_rain, precip_snow, evap, tsurf, fder, albsol, taux, tauy,&
       sens_prec_liq, sens_prec_sol, lat_prec_liq, lat_prec_sol)
!
! This subroutine cumulates some fields for each time-step during a coupling 
! period. At last time-step in a coupling period the fields are transformed to the 
! grid accepted by the coupler. All fields for all types of surfaces are sent to
! the coupler.
!
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl
    USE indice_sol_mod
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                     :: itime
    INTEGER, INTENT(IN)                     :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)    :: knindex
    REAL, INTENT(IN)                        :: dtime
    REAL, DIMENSION(klon), INTENT(IN)       :: rlon, rlat
    REAL, DIMENSION(klon), INTENT(IN)       :: swdown, lwdown 
    REAL, DIMENSION(klon), INTENT(IN)       :: fluxlat, fluxsens
    REAL, DIMENSION(klon), INTENT(IN)       :: precip_rain, precip_snow
    REAL, DIMENSION(klon), INTENT(IN)       :: evap, tsurf, fder
    REAL, DIMENSION(klon), INTENT(IN)       :: albsol, taux, tauy
    REAL, DIMENSION(klon,nbsrf), INTENT(IN) :: pctsrf
    REAL, INTENT(IN):: sens_prec_liq(:), sens_prec_sol(:) ! (knon)
    REAL, DIMENSION(klon), INTENT(IN)       :: lat_prec_liq, lat_prec_sol
    LOGICAL, INTENT(IN)                     :: lafin

! Local variables
!*************************************************************************************
    INTEGER                                 :: cpl_index, ig 
    INTEGER                                 :: error, sum_error
    CHARACTER(len = 25)                     :: modname = 'cpl_send_seaice_fields'
    CHARACTER(len = 80)                     :: abort_message
    REAL, DIMENSION(klon)                   :: cpl_fder_tmp

!*************************************************************************************
! Start calulation
! The sea-ice points are saved with second array index=2
!
!*************************************************************************************
    cpl_index = 2

!*************************************************************************************
! Reset fields to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
       cpl_sols(1:knon,cpl_index) = 0.0
       cpl_nsol(1:knon,cpl_index) = 0.0
       cpl_rain(1:knon,cpl_index) = 0.0
       cpl_snow(1:knon,cpl_index) = 0.0
       cpl_evap(1:knon,cpl_index) = 0.0
       cpl_tsol(1:knon,cpl_index) = 0.0
       cpl_fder(1:knon,cpl_index) = 0.0
       cpl_albe(1:knon,cpl_index) = 0.0
       cpl_taux(1:knon,cpl_index) = 0.0
       cpl_tauy(1:knon,cpl_index) = 0.0
       cpl_sens_rain(1:knon,cpl_index) = 0.0
       cpl_sens_snow(1:knon,cpl_index) = 0.0
       cpl_taumod(1:knon,cpl_index) = 0.0
    ENDIF
       
!*************************************************************************************
! Cumulate at each time-step
!
!*************************************************************************************    
    DO ig = 1, knon
       cpl_sols(ig,cpl_index) = cpl_sols(ig,cpl_index) + &
            swdown(ig)      / REAL(nexca)
       cpl_nsol(ig,cpl_index) = cpl_nsol(ig,cpl_index) + &
            (lwdown(ig) + fluxlat(ig) +fluxsens(ig)) / REAL(nexca)
       cpl_rain(ig,cpl_index) = cpl_rain(ig,cpl_index) + &
            precip_rain(ig) / REAL(nexca)
       cpl_snow(ig,cpl_index) = cpl_snow(ig,cpl_index) + &
            precip_snow(ig) / REAL(nexca)
       cpl_evap(ig,cpl_index) = cpl_evap(ig,cpl_index) + &
            evap(ig)        / REAL(nexca)
       cpl_tsol(ig,cpl_index) = cpl_tsol(ig,cpl_index) + &
            tsurf(ig)       / REAL(nexca)
       cpl_fder(ig,cpl_index) = cpl_fder(ig,cpl_index) + &
            fder(ig)        / REAL(nexca)
       cpl_albe(ig,cpl_index) = cpl_albe(ig,cpl_index) + &
            albsol(ig)      / REAL(nexca)
       cpl_taux(ig,cpl_index) = cpl_taux(ig,cpl_index) + &
            taux(ig)        / REAL(nexca)
       cpl_tauy(ig,cpl_index) = cpl_tauy(ig,cpl_index) + &
            tauy(ig)        / REAL(nexca)     
       cpl_sens_rain(ig,cpl_index) = cpl_sens_rain(ig,cpl_index) + &
            sens_prec_liq(ig)      / REAL(nexca)
       cpl_sens_snow(ig,cpl_index) = cpl_sens_snow(ig,cpl_index) + &
            sens_prec_sol(ig)      / REAL(nexca)
       cpl_taumod(ig,cpl_index) = cpl_taumod(ig,cpl_index) + &
            SQRT ( taux(ig)*taux(ig)+tauy(ig)*tauy(ig) ) / REAL(nexca) 
    ENDDO

!*************************************************************************************
! If the time-step corresponds to the end of coupling period the 
! fields are transformed to the 2D grid and all fields are sent to coupler.
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 0) THEN
       IF (.NOT. ALLOCATED(cpl_sols2D)) THEN
          sum_error = 0
          ALLOCATE(cpl_sols2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_nsol2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_rain2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_snow2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_evap2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tsol2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_fder2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_albe2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taux2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_tauy2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_windsp2D(nbp_lon,jj_nb), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_sens_rain2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_sens_snow2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error
          ALLOCATE(cpl_taumod2D(nbp_lon,jj_nb,2), stat=error)
          sum_error = sum_error + error

          IF (carbon_cycle_cpl) THEN
             ALLOCATE(cpl_atm_co22D(nbp_lon,jj_nb), stat=error)
             sum_error = sum_error + error
          ENDIF

          IF (sum_error /= 0) THEN
             abort_message='Pb allocation variables couplees pour l''ecriture'
             CALL abort_physic(modname,abort_message,1)
          ENDIF
       ENDIF

       CALL gath2cpl(cpl_sols(:,cpl_index), cpl_sols2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_nsol(:,cpl_index), cpl_nsol2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_rain(:,cpl_index), cpl_rain2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_snow(:,cpl_index), cpl_snow2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_evap(:,cpl_index), cpl_evap2D(:,:,cpl_index), &
            knon, knindex)

! cpl_tsol2D(:,:,:) not used!
       CALL gath2cpl(cpl_tsol(:,cpl_index), cpl_tsol2D(:,:, cpl_index), &
            knon, knindex)

       ! Set default value and decompress before gath2cpl
       cpl_fder_tmp(:) = -20.
       DO ig = 1, knon
          cpl_fder_tmp(knindex(ig))=cpl_fder(ig,cpl_index)
       ENDDO
       CALL gath2cpl(cpl_fder_tmp(:), cpl_fder2D(:,:,cpl_index), &
            klon, unity)

! cpl_albe2D(:,:,:) not used!
       CALL gath2cpl(cpl_albe(:,cpl_index), cpl_albe2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taux(:,cpl_index), cpl_taux2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_tauy(:,cpl_index), cpl_tauy2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_sens_rain(:,cpl_index), cpl_sens_rain2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_sens_snow(:,cpl_index), cpl_sens_snow2D(:,:,cpl_index), &
            knon, knindex)

       CALL gath2cpl(cpl_taumod(:,cpl_index), cpl_taumod2D(:,:,cpl_index), &
            knon, knindex)

       ! Send all fields
       CALL cpl_send_all(itime, dtime, pctsrf, lafin, rlon, rlat)
    ENDIF

  END SUBROUTINE cpl_send_seaice_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_land_fields(itime, knon, knindex, rriv_in, rcoa_in)
!
! This subroutine cumulates some fields for each time-step during a coupling 
! period. At last time-step in a coupling period the fields are transformed to the 
! grid accepted by the coupler. No sending to the coupler will be done from here 
! (it is done in cpl_send_seaice_fields).
!
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat

! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                       :: itime
    INTEGER, INTENT(IN)                       :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex
    REAL, DIMENSION(klon), INTENT(IN)         :: rriv_in
    REAL, DIMENSION(klon), INTENT(IN)         :: rcoa_in

! Local variables
!*************************************************************************************
    REAL, DIMENSION(nbp_lon,jj_nb)             :: rriv2D
    REAL, DIMENSION(nbp_lon,jj_nb)             :: rcoa2D

!*************************************************************************************
! Rearrange fields in 2D variables 
! First initialize to zero to avoid unvalid points causing problems
!
!*************************************************************************************
!$OMP MASTER
    rriv2D(:,:) = 0.0
    rcoa2D(:,:) = 0.0
!$OMP END MASTER
    CALL gath2cpl(rriv_in, rriv2D, knon, knindex)
    CALL gath2cpl(rcoa_in, rcoa2D, knon, knindex)

!*************************************************************************************
! Reset cumulated fields to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
!$OMP MASTER
       cpl_rriv2D(:,:) = 0.0
       cpl_rcoa2D(:,:) = 0.0
!$OMP END MASTER
    ENDIF

!*************************************************************************************
! Cumulate : Following fields should be cumulated at each time-step
!
!*************************************************************************************    
!$OMP MASTER
    cpl_rriv2D(:,:) = cpl_rriv2D(:,:) + rriv2D(:,:) / REAL(nexca)
    cpl_rcoa2D(:,:) = cpl_rcoa2D(:,:) + rcoa2D(:,:) / REAL(nexca)
!$OMP END MASTER

  END SUBROUTINE cpl_send_land_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_landice_fields(itime, knon, knindex, rlic_in, rlic_in_frac)
! This subroutine cumulates the field for melting ice for each time-step 
! during a coupling period. This routine will not send to coupler. Sending 
! will be done in cpl_send_seaice_fields.
!

    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat

! Input varibales
!*************************************************************************************
    INTEGER, INTENT(IN)                       :: itime
    INTEGER, INTENT(IN)                       :: knon
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex
    REAL, DIMENSION(klon), INTENT(IN)         :: rlic_in
    REAL, DIMENSION(klon), INTENT(IN)         :: rlic_in_frac  ! Fraction for continental ice, can be equal to 
                                                               ! pctsrf(:,is_lic) or not, depending on landice_opt
    

! Local varibales
!*************************************************************************************
    REAL, DIMENSION(nbp_lon,jj_nb)             :: rlic2D

!*************************************************************************************
! Rearrange field in a 2D variable 
! First initialize to zero to avoid unvalid points causing problems
!
!*************************************************************************************
!$OMP MASTER
    rlic2D(:,:) = 0.0
!$OMP END MASTER
    CALL gath2cpl(rlic_in, rlic2D, knon, knindex)
    CALL gath2cpl(rlic_in_frac(:), rlic_in_frac2D(:,:), knon, knindex) 
!*************************************************************************************
! Reset field to zero in the beginning of a new coupling period 
!
!*************************************************************************************
    IF (MOD(itime, nexca) == 1) THEN
!$OMP MASTER
       cpl_rlic2D(:,:) = 0.0
!$OMP END MASTER
    ENDIF

!*************************************************************************************
! Cumulate : Melting ice should be cumulated at each time-step
!
!*************************************************************************************    
!$OMP MASTER
    cpl_rlic2D(:,:) = cpl_rlic2D(:,:) + rlic2D(:,:) / REAL(nexca)
!$OMP END MASTER

  END SUBROUTINE cpl_send_landice_fields

!
!*************************************************************************************
!

  SUBROUTINE cpl_send_all(itime, dtime, pctsrf, lafin, rlon, rlat)
! This routine will send fields for all different surfaces to the coupler.
! This subroutine should be executed after calculations by the last surface(sea-ice),
! all calculations at the different surfaces have to be done before. 
!    
    USE surface_data
    USE carbon_cycle_mod, ONLY : carbon_cycle_cpl
    USE indice_sol_mod
    USE print_control_mod, ONLY: lunout
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
    USE time_phylmdz_mod, ONLY: start_time, itau_phy
    USE config_ocean_skin_m, only: activate_ocean_skin
    USE lmdz_mpi

! Some includes
!    
! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                                  :: itime
    REAL, INTENT(IN)                                     :: dtime
    REAL, DIMENSION(klon), INTENT(IN)                    :: rlon, rlat
    REAL, DIMENSION(klon,nbsrf), INTENT(IN)              :: pctsrf
    LOGICAL, INTENT(IN)                                  :: lafin
    
! Local variables
!*************************************************************************************
    INTEGER                                              :: error, sum_error, i,j,k
    INTEGER                                              :: istart, iend
    INTEGER                                              :: wstart, wend
    INTEGER                                              :: itau_w
    INTEGER                                              :: time_sec
    INTEGER, DIMENSION(nbp_lon*(nbp_lat))                      :: ndexct
    REAL                                                 :: Up, Down
    REAL, DIMENSION(nbp_lon, jj_nb)                          :: tmp_lon, tmp_lat
    REAL, DIMENSION(nbp_lon, jj_nb, 4)                       :: pctsrf2D
    REAL, DIMENSION(nbp_lon, jj_nb)                          :: deno
    CHARACTER(len = 20)                                  :: modname = 'cpl_send_all'
    CHARACTER(len = 80)                                  :: abort_message
    LOGICAL                                              :: checkout

! Variables with fields to coupler
    REAL, DIMENSION(nbp_lon, jj_nb)                          :: tmp_taux
    REAL, DIMENSION(nbp_lon, jj_nb)                          :: tmp_tauy
    REAL, DIMENSION(nbp_lon, jj_nb)                          :: tmp_calv
    REAL, DIMENSION(nbp_lon*jj_nb,1)                         :: field
! Table with all fields to send to coupler
    REAL, DIMENSION(nbp_lon, jj_nb, maxsend_phys)            :: tab_flds
    REAL, DIMENSION(klon_mpi)                                :: rlon_mpi, rlat_mpi
    REAL  :: calving(nb_zone_calving)
    REAL  :: calving_glo(nb_zone_calving)
    
    INTEGER, DIMENSION(MPI_STATUS_SIZE)                  :: status

! End definitions
!*************************************************************************************
    


!*************************************************************************************
! All fields are stored in a table tab_flds(:,:,:)
! First store the fields which are already on the right format
!
!*************************************************************************************
!$OMP MASTER
    tab_flds(:,:,ids_windsp) = cpl_windsp2D(:,:)
    tab_flds(:,:,ids_shfice) = cpl_sols2D(:,:,2)
    tab_flds(:,:,ids_nsfice) = cpl_nsol2D(:,:,2)
    tab_flds(:,:,ids_dflxdt) = cpl_fder2D(:,:,2)
    tab_flds(:,:,ids_qraioc) = cpl_sens_rain2D(:,:,1)
    tab_flds(:,:,ids_qsnooc) = cpl_sens_snow2D(:,:,1)
    tab_flds(:,:,ids_qraiic) = cpl_sens_rain2D(:,:,2)
    tab_flds(:,:,ids_qsnoic) = cpl_sens_snow2D(:,:,2)

    if (activate_ocean_skin == 2) then
       tab_flds(:, :, ids_delta_sst) = cpl_delta_sst_2D
       tab_flds(:, :, ids_delta_sal) = cpl_delta_sal_2D
       tab_flds(:, :, ids_dter) = cpl_dter_2D
       tab_flds(:, :, ids_dser) = cpl_dser_2D
       tab_flds(:, :, ids_dt_ds) = cpl_dt_ds_2D
    end if
    
    IF (version_ocean=='nemo') THEN
       tab_flds(:,:,ids_liqrun) = (cpl_rriv2D(:,:) + cpl_rcoa2D(:,:))
       IF (carbon_cycle_cpl) tab_flds(:,:,ids_atmco2)=cpl_atm_co22D(:,:)
    ELSE IF (version_ocean=='opa8') THEN
       tab_flds(:,:,ids_shfoce) = cpl_sols2D(:,:,1)
       tab_flds(:,:,ids_nsfoce) = cpl_nsol2D(:,:,1)
       tab_flds(:,:,ids_icevap) = cpl_evap2D(:,:,2)
       tab_flds(:,:,ids_ocevap) = cpl_evap2D(:,:,1)
       tab_flds(:,:,ids_runcoa) = cpl_rcoa2D(:,:)
       tab_flds(:,:,ids_rivflu) = cpl_rriv2D(:,:)
    ENDIF

!*************************************************************************************
! Transform the fraction of sub-surfaces from 1D to 2D array
!
!*************************************************************************************
    pctsrf2D(:,:,:) = 0.
!$OMP END MASTER
    CALL gath2cpl(pctsrf(:,is_oce), pctsrf2D(:,:,is_oce), klon, unity)
    CALL gath2cpl(pctsrf(:,is_sic), pctsrf2D(:,:,is_sic), klon, unity)




!*************************************************************************************
! Calculate the average calving per latitude
! Store calving in tab_flds(:,:,19)
! 
!*************************************************************************************      
    IF (is_omp_root) THEN

      IF (cpl_old_calving) THEN   ! use old calving

        DO j = 1, jj_nb
           tmp_calv(:,j) = DOT_PRODUCT (cpl_rlic2D(1:nbp_lon,j), &
                rlic_in_frac2D(1:nbp_lon,j)) / REAL(nbp_lon)
        ENDDO
    
    
        IF (is_parallel) THEN
           IF (.NOT. is_north_pole_dyn) THEN
              CALL MPI_RECV(Up,1,MPI_REAL_LMDZ,mpi_rank-1,1234,COMM_LMDZ_PHY,status,error)
              CALL MPI_SEND(tmp_calv(1,1),1,MPI_REAL_LMDZ,mpi_rank-1,1234,COMM_LMDZ_PHY,error)
           ENDIF
       
           IF (.NOT. is_south_pole_dyn) THEN
              CALL MPI_SEND(tmp_calv(1,jj_nb),1,MPI_REAL_LMDZ,mpi_rank+1,1234,COMM_LMDZ_PHY,error)
              CALL MPI_RECV(down,1,MPI_REAL_LMDZ,mpi_rank+1,1234,COMM_LMDZ_PHY,status,error)
           ENDIF
         
           IF (.NOT. is_north_pole_dyn .AND. ii_begin /=1) THEN
              Up=Up+tmp_calv(nbp_lon,1)
              tmp_calv(:,1)=Up
           ENDIF
           
           IF (.NOT. is_south_pole_dyn .AND. ii_end /= nbp_lon) THEN
              Down=Down+tmp_calv(1,jj_nb)
              tmp_calv(:,jj_nb)=Down
           ENDIF
        ENDIF
        tab_flds(:,:,ids_calvin) = tmp_calv(:,:)

      ELSE
         ! cpl_old_calving=FALSE
         ! To be used with new method for calculation of coupling weights
         DO k=1,nb_zone_calving
            calving(k)=0
            DO j = 1, jj_nb
               calving(k)= calving(k)+DOT_PRODUCT(cpl_rlic2D(:,j)*area_calving(:,j,k),rlic_in_frac2D(:,j))
            ENDDO
         ENDDO
         
         CALL MPI_ALLREDUCE(calving, calving_glo, nb_zone_calving, MPI_REAL_LMDZ, MPI_SUM, COMM_LMDZ_PHY, error)
         
         tab_flds(:,:,ids_calvin) = 0
         DO k=1,nb_zone_calving
            IF (ind_calving(k)>0 ) THEN
               j=(ind_calving(k)-1)/nbp_lon + 1
               i=MOD(ind_calving(k)-1,nbp_lon)+1
               tab_flds(i,j,ids_calvin) = calving_glo(k)
            ENDIF
         ENDDO
         
      ENDIF
      
!*************************************************************************************
! Calculate total flux for snow, rain and wind with weighted addition using the 
! fractions of ocean and seaice.
!
!*************************************************************************************    
       ! fraction oce+seaice
       deno =  pctsrf2D(:,:,is_oce) + pctsrf2D(:,:,is_sic) 

       IF (version_ocean=='nemo') THEN
          tab_flds(:,:,ids_shftot)  = 0.0
          tab_flds(:,:,ids_nsftot) = 0.0
          tab_flds(:,:,ids_totrai) = 0.0
          tab_flds(:,:,ids_totsno) = 0.0
          tab_flds(:,:,ids_toteva) = 0.0
          tab_flds(:,:,ids_taumod) = 0.0
  
          tmp_taux(:,:)    = 0.0
          tmp_tauy(:,:)    = 0.0
          ! For all valid grid cells containing some fraction of ocean or sea-ice
          WHERE ( deno(:,:) /= 0 )
             tmp_taux = cpl_taux2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_taux2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tmp_tauy = cpl_tauy2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_tauy2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)

             tab_flds(:,:,ids_shftot) = cpl_sols2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_sols2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_nsftot) = cpl_nsol2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_nsol2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_totrai) = cpl_rain2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_rain2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_totsno) = cpl_snow2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_snow2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_toteva) = cpl_evap2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_evap2D(:,:,2)  * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_taumod) = cpl_taumod2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_taumod2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             
         ENDWHERE

          tab_flds(:,:,ids_icevap) = cpl_evap2D(:,:,2) 
          
       ELSE IF (version_ocean=='opa8') THEN
          ! Store fields for rain and snow in tab_flds(:,:,15) and tab_flds(:,:,16)
          tab_flds(:,:,ids_totrai) = 0.0
          tab_flds(:,:,ids_totsno) = 0.0
          tmp_taux(:,:)    = 0.0
          tmp_tauy(:,:)    = 0.0
          ! For all valid grid cells containing some fraction of ocean or sea-ice
          WHERE ( deno(:,:) /= 0 )
             tab_flds(:,:,ids_totrai) = cpl_rain2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_rain2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tab_flds(:,:,ids_totsno) = cpl_snow2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_snow2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             
             tmp_taux = cpl_taux2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_taux2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
             tmp_tauy = cpl_tauy2D(:,:,1) * pctsrf2D(:,:,is_oce) / deno(:,:) +    &
                  cpl_tauy2D(:,:,2) * pctsrf2D(:,:,is_sic) / deno(:,:)
          ENDWHERE
       ENDIF

    ENDIF ! is_omp_root
  
!*************************************************************************************
! Transform the wind components from local atmospheric 2D coordinates to geocentric 
! 3D coordinates. 
! Store the resulting wind components in tab_flds(:,:,1:6)
!*************************************************************************************

! Transform the longitudes and latitudes on 2D arrays
    
    CALL gather_omp(rlon,rlon_mpi)
    CALL gather_omp(rlat,rlat_mpi)
!$OMP MASTER
    CALL Grid1DTo2D_mpi(rlon_mpi,tmp_lon)
    CALL Grid1DTo2D_mpi(rlat_mpi,tmp_lat)
!$OMP END MASTER    

    IF (is_sequential) THEN
       IF (is_north_pole_dyn) tmp_lon(:,1)     = tmp_lon(:,2)
       IF (is_south_pole_dyn) tmp_lon(:,nbp_lat) = tmp_lon(:,nbp_lat-1)
    ENDIF
      
! NetCDF output of the wind before transformation of coordinate system
    IF (is_sequential) THEN
       ndexct(:) = 0
       itau_w = itau_phy + itime + start_time * day_step_phy
       CALL histwrite(nidct,'tauxe',itau_w,tmp_taux,nbp_lon*(nbp_lat),ndexct)
       CALL histwrite(nidct,'tauyn',itau_w,tmp_tauy,nbp_lon*(nbp_lat),ndexct)
       CALL histwrite(nidct,'tmp_lon',itau_w,tmp_lon,nbp_lon*(nbp_lat),ndexct)
       CALL histwrite(nidct,'tmp_lat',itau_w,tmp_lat,nbp_lon*(nbp_lat),ndexct)
    ENDIF

! Transform the wind from spherical atmospheric 2D coordinates to geocentric
! cartesian 3D coordinates 
!$OMP MASTER
    CALL atm2geo (nbp_lon, jj_nb, tmp_taux, tmp_tauy, tmp_lon, tmp_lat, &
         tab_flds(:,:,ids_tauxxu), tab_flds(:,:,ids_tauyyu), tab_flds(:,:,ids_tauzzu) )
    
    tab_flds(:,:,ids_tauxxv)  = tab_flds(:,:,ids_tauxxu)
    tab_flds(:,:,ids_tauyyv)  = tab_flds(:,:,ids_tauyyu)
    tab_flds(:,:,ids_tauzzv)  = tab_flds(:,:,ids_tauzzu)
!$OMP END MASTER

!*************************************************************************************
! NetCDF output of all fields just before sending to coupler.
!
!*************************************************************************************
    IF (is_sequential) THEN
        DO j=1,maxsend_phys
          IF (infosend(midcpl)%fld(j)%action) CALL histwrite(nidct,infosend(midcpl)%fld(j)%name, itau_w, &
             tab_flds(:,:,j),nbp_lon*(nbp_lat),ndexct)
        ENDDO
    ENDIF
!*************************************************************************************
! Send the table of all fields
!
!*************************************************************************************
    time_sec=(itime-1)*dtime
#ifdef CPP_COUPLE
!$OMP MASTER
! ======================================================================
! L. Fairhead (09/2003) adapted From L.Z.X Li: this subroutine provides the
! atmospheric coupling fields to the coupler with the psmile library.
! IF last time step, writes output fields to binary files.
! ======================================================================
    checkout=.FALSE.

    WRITE(lunout,*) ' '
    WRITE(lunout,*) 'Intocpl: sending fields to CPL, time_sec= ', time_sec
    WRITE(lunout,*) 'last = ', lafin
    WRITE(lunout,*)

    istart=ii_begin
    IF (is_south_pole_dyn) THEN
       iend=(jj_end-jj_begin)*nbp_lon+nbp_lon
    ELSE
       iend=(jj_end-jj_begin)*nbp_lon+ii_end
    ENDIF

    IF (checkout) THEN
       wstart=istart
       wend=iend
       IF (is_north_pole_dyn) wstart=istart+nbp_lon-1
       IF (is_south_pole_dyn) wend=iend-nbp_lon+1

       DO i = 1, maxsend_phys
          IF (infosend(midcpl)%fld(i)%action) THEN
             field(:,1) = RESHAPE(tab_flds(:,:,i),(/nbp_lon*jj_nb/))
             CALL writefield_phy(infosend(midcpl)%fld(i)%name,field(wstart:wend,1),1)
          END IF
       END DO
    END IF

    DO i = 1, maxsend_phys
      IF (infosend(midcpl)%fld(i)%action .AND. infosend(midcpl)%fld(i)%nid .NE. -1 ) THEN
          field(:,1) = RESHAPE(tab_flds(:,:,i),(/nbp_lon*jj_nb/))
          CALL cpl_snd(midcpl, i , time_sec, field(istart:iend,:))
      ENDIF
    END DO

!************************************************************************************
! Finalize PSMILE for the case is_sequential, if parallel finalization is done
! from Finalize_parallel in dyn3dpar/parallel.F90
!************************************************************************************

    IF (lafin) THEN
       IF (is_sequential) THEN
          CALL prism_terminate_proto(error)
          IF (error .NE. PRISM_Ok) THEN
             abort_message=' Problem in prism_terminate_proto '
             CALL abort_physic(modname,abort_message,1)
          ENDIF
       ENDIF
    ENDIF
!$OMP END MASTER
#endif

!*************************************************************************************
! Finish with some dellocate
!
!*************************************************************************************  
    sum_error=0
    DEALLOCATE(cpl_sols2D, cpl_nsol2D, cpl_rain2D, cpl_snow2D, stat=error )
    sum_error = sum_error + error
    DEALLOCATE(cpl_evap2D, cpl_tsol2D, cpl_fder2D, cpl_albe2D, stat=error )
    sum_error = sum_error + error
    DEALLOCATE(cpl_taux2D, cpl_tauy2D, cpl_windsp2D, cpl_taumod2D, stat=error )
    sum_error = sum_error + error
    DEALLOCATE(cpl_sens_rain2D, cpl_sens_snow2D, stat=error)
    sum_error = sum_error + error

    
    IF (carbon_cycle_cpl) THEN
       DEALLOCATE(cpl_atm_co22D, stat=error )
       sum_error = sum_error + error
    ENDIF

    if (activate_ocean_skin == 2) deallocate(cpl_delta_sst_2d, &
         cpl_delta_sal_2d, cpl_dter_2d, cpl_dser_2d, cpl_dt_ds_2d)

    IF (sum_error /= 0) THEN
       abort_message='Pb in deallocation of cpl_xxxx2D coupling variables'
       CALL abort_physic(modname,abort_message,1)
    ENDIF
    
  END SUBROUTINE cpl_send_all
!
!*************************************************************************************
!
  SUBROUTINE cpl2gath(champ_in, champ_out, knon, knindex)
  USE mod_phys_lmdz_para
! Cette routine transforme un champs de la grille 2D recu du coupleur sur la grille 
! 'gathered' (la grille physiq comprime).
!
! 
! input:         
!   champ_in     champ sur la grille 2D
!   knon         nombre de points dans le domaine a traiter
!   knindex      index des points de la surface a traiter
!
! output:
!   champ_out    champ sur la grille 'gatherd'
!
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat

! Input
    INTEGER, INTENT(IN)                       :: knon
    REAL, DIMENSION(nbp_lon,jj_nb), INTENT(IN)    :: champ_in
    INTEGER, DIMENSION(klon), INTENT(IN)      :: knindex

! Output
    REAL, DIMENSION(klon_mpi), INTENT(OUT)        :: champ_out

! Local
    INTEGER                                   :: i, ig
    REAL, DIMENSION(klon_mpi)                 :: temp_mpi
    REAL, DIMENSION(klon)                     :: temp_omp

!*************************************************************************************
!
    

! Transform from 2 dimensions (nbp_lon,jj_nb) to 1 dimension (klon)
!$OMP MASTER 
    CALL Grid2Dto1D_mpi(champ_in,temp_mpi)
!$OMP END MASTER

    CALL scatter_omp(temp_mpi,temp_omp)
    
! Compress from klon to knon
    DO i = 1, knon
       ig = knindex(i)
       champ_out(i) = temp_omp(ig)
    ENDDO

  END SUBROUTINE cpl2gath
!
!*************************************************************************************
!
  SUBROUTINE gath2cpl(champ_in, champ_out, knon, knindex)
  USE mod_phys_lmdz_para
! Cette routine ecrit un champ 'gathered' sur la grille 2D pour le passer
! au coupleur.
!
! input:         
!   champ_in     champ sur la grille gathere        
!   knon         nombre de points dans le domaine a traiter
!   knindex      index des points de la surface a traiter
!
! output:
!   champ_out    champ sur la grille 2D
!
    USE mod_grid_phy_lmdz, ONLY : nbp_lon, nbp_lat
    
! Input arguments
!*************************************************************************************
    INTEGER, INTENT(IN)                    :: knon
    REAL, DIMENSION(klon), INTENT(IN)      :: champ_in
    INTEGER, DIMENSION(klon), INTENT(IN)   :: knindex

! Output arguments
!*************************************************************************************
    REAL, DIMENSION(nbp_lon,jj_nb), INTENT(OUT) :: champ_out

! Local variables
!*************************************************************************************
    INTEGER                                :: i, ig
    REAL, DIMENSION(klon)                  :: temp_omp
    REAL, DIMENSION(klon_mpi)              :: temp_mpi
!*************************************************************************************

! Decompress from knon to klon
    temp_omp = 0.
    DO i = 1, knon
       ig = knindex(i)
       temp_omp(ig) = champ_in(i)
    ENDDO

! Transform from 1 dimension (klon) to 2 dimensions (nbp_lon,jj_nb)
    CALL gather_omp(temp_omp,temp_mpi)

!$OMP MASTER    
    CALL Grid1Dto2D_mpi(temp_mpi,champ_out)
    
    IF (is_north_pole_dyn) champ_out(:,1)=temp_mpi(1)
    IF (is_south_pole_dyn) champ_out(:,jj_nb)=temp_mpi(klon)
!$OMP END MASTER
    
  END SUBROUTINE gath2cpl
!
!*************************************************************************************
!
END MODULE cpl_mod

