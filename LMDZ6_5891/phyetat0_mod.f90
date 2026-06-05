! $Id: phyetat0_mod.f90 5776 2025-07-15 12:13:19Z evignon $

MODULE phyetat0_mod

  PRIVATE
  PUBLIC :: phyetat0

CONTAINS

SUBROUTINE phyetat0 (fichnom, clesphy0, tabcntr0)

  USE clesphys_mod_h
  USE dimphy, only: klon, zmasq, klev, nbtersrf, nbtsoildepths
  USE iophy, ONLY : init_iophy_new
  USE ocean_cpl_mod,    ONLY : ocean_cpl_init
  USE fonte_neige_mod,  ONLY : fonte_neige_init
  USE pbl_surface_mod,  ONLY : pbl_surface_init
!GG  USE surface_data,     ONLY : type_ocean, version_ocean
  USE surface_data,     ONLY : type_ocean, version_ocean, iflag_seaice, &
                                   iflag_seaice_alb, iflag_leads
!GG
  USE phyetat0_get_mod, ONLY : phyetat0_get, phyetat0_srf
  USE phys_state_var_mod, ONLY : ancien_ok, clwcon, detr_therm, phys_tstep, &
       qsol, fevap, z0m, z0h, agesno, &
       du_gwd_rando, du_gwd_front, entr_therm, f0, fm_therm, &
       falb_dir, falb_dif, prw_ancien, prlw_ancien, prsw_ancien, prbsw_ancien, &
       ftsol, pbl_tke, pctsrf, q_ancien, ql_ancien, qs_ancien, qbs_ancien, &
       cf_ancien, rvc_ancien, tke_ancien, radpas, radsol, rain_fall, ratqs, &
       rnebcon, rugoro, sig1, snow_fall, bs_fall, solaire_etat0, sollw, sollwdown, &
       solsw, solswfdiff, t_ancien, u_ancien, v_ancien, w01, wake_cstar, wake_deltaq, &
       wake_deltat, wake_delta_pbl_TKE, delta_tsurf, beta_aridity, wake_fip, wake_pe, &
       wake_s, awake_s, wake_dens, awake_dens, cv_gen, zgam, zmax0, zmea, zpic, zsig, &
       zstd, zthe, zval, ale_bl, ale_bl_trig, alp_bl, u10m, v10m, treedrg, &
       ale_wake, ale_bl_stat, ds_ns, dt_ns, delta_sst, delta_sal, dter, dser, &
!GG       dt_ds, ratqs_inter_
       dt_ds, ratqs_inter_, &
       hice, tice, bilg_cumul, &
!GG
       frac_tersrf, z0m_tersrf, ratio_z0m_z0h_tersrf, &
       albedo_tersrf, beta_tersrf, inertie_tersrf, alpha_soil_tersrf, &
       period_tersrf, hcond_tersrf, tsurfi_tersrf, tsoili_tersrf, tsoil_depth, &
       qsurf_tersrf, tsurf_tersrf, tsoil_tersrf, tsurf_new_tersrf, cdragm_tersrf, &
       cdragh_tersrf, swnet_tersrf, lwnet_tersrf, fluxsens_tersrf, fluxlat_tersrf
!FC
  USE geometry_mod,     ONLY: longitude_deg, latitude_deg
  USE iostart,          ONLY: close_startphy, get_field, get_var, open_startphy
  USE infotrac_phy,     ONLY: nqtot, nbtr, type_trac, tracers, new2oldH2O
  USE strings_mod,      ONLY: maxlen
  USE traclmdz_mod,     ONLY: traclmdz_from_restart
  USE carbon_cycle_mod, ONLY: carbon_cycle_init, carbon_cycle_cpl, carbon_cycle_tr, carbon_cycle_rad, co2_send, RCO2_glo
  USE indice_sol_mod,   ONLY: nbsrf, is_ter, epsfra, is_lic, is_oce, is_sic
  !GG USE ocean_slab_mod,   ONLY: nslay, tslab, seaice, tice, ocean_slab_init
  USE ocean_slab_mod,   ONLY: nslay, tslab, seaice, tice_slab, ocean_slab_init
  !GG
  USE time_phylmdz_mod, ONLY: init_iteration, pdtphys, itau_phy
  use wxios_mod, ONLY: missing_val_xios => missing_val, using_xios
  use netcdf, only: missing_val_netcdf => nf90_fill_real
  use config_ocean_skin_m, only: activate_ocean_skin
  USE surf_param_mod, ONLY: average_surf_var, interpol_tsoil !AM
  USE dimsoil_mod_h, ONLY: nsoilmx
  USE yomcst_mod_h
  USE alpale_mod
  USE compbl_mod_h
  USE oasis
  USE cpl_mod, ONLY: cpl_inca
  USE pycpl
  USE pyfld
IMPLICIT none
  !======================================================================
  ! Auteur(s) Z.X. Li (LMD/CNRS) date: 19930818
  ! Objet: Lecture de l'etat initial pour la physique
  !======================================================================
  CHARACTER*(*) fichnom

  ! les variables globales lues dans le fichier restart

  REAL tsoil(klon, nsoilmx, nbsrf)
  REAL qsurf(klon, nbsrf)
  REAL snow(klon, nbsrf)
  real fder(klon)
  REAL run_off_lic_0(klon)
  REAL fractint(klon)
  REAL trs(klon, nbtr)
  REAL zts(klon)
  ! pour drag arbres FC
  REAL drg_ter(klon,klev)

  CHARACTER*6 ocean_in
  LOGICAL ok_veget_in

  INTEGER        longcles
  PARAMETER    ( longcles = 20 )
  REAL clesphy0( longcles )

  REAL xmin, xmax

  INTEGER nid, nvarid
  INTEGER ierr, i, nsrf, isoil , k
  INTEGER length
  PARAMETER (length=100)
  INTEGER it, iq, isw
  REAL tab_cntrl(length), tabcntr0(length)
  CHARACTER*7 str7
  CHARACTER*2 str2
  LOGICAL :: found
  REAL :: lon_startphy(klon), lat_startphy(klon)
  CHARACTER(LEN=maxlen) :: tname, t(2)
  REAL :: missing_val

  IF (using_xios) THEN
    missing_val=missing_val_xios
  ELSE
    missing_val=missing_val_netcdf
  ENDIF
  
  ! FH1D
  !     real iolat(jjm+1)
  !real iolat(jjm+1-1/(iim*jjm))

  ! Ouvrir le fichier contenant l'etat initial:

  CALL open_startphy(fichnom)

  ! Lecture des parametres de controle:

  CALL get_var("controle", tab_cntrl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! FH 2008/05/09 On elimine toutes les clefs physiques dans la dynamique
  ! Les constantes de la physiques sont lues dans la physique seulement.
  ! Les egalites du type
  !             tab_cntrl( 5 )=clesphy0(1)
  ! sont remplacees par
  !             clesphy0(1)=tab_cntrl( 5 )
  ! On inverse aussi la logique.
  ! On remplit les tab_cntrl avec les parametres lus dans les .def
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i = 1, length
     tabcntr0( i ) = tab_cntrl( i )
  ENDDO

  tab_cntrl(1)=pdtphys
  tab_cntrl(2)=radpas

  ! co2_ppm : value from the previous time step

  ! co2_ppm0 : initial value of atmospheric CO2 (from create_etat0_limit.e .def)
  co2_ppm0 = 284.32
  ! when no initial value is available e.g., from a restart 
  ! this variable must be set  in a .def file which will then be 
  ! used by the conf_phys_m.F90 routine.
  ! co2_ppm0 = 284.32 (illustrative example on how to set the variable in .def
  ! file, for a pre-industrial CO2 concentration value)

  IF (carbon_cycle_tr .OR. carbon_cycle_cpl) THEN
     co2_ppm = tab_cntrl(3)
     RCO2    = co2_ppm * 1.0e-06 * RMCO2 / RMD
     IF (tab_cntrl(17) > 0. .AND. carbon_cycle_rad) THEN
           RCO2_glo = tab_cntrl(17)
       ELSE
           RCO2_glo    = co2_ppm0 * 1.0e-06 * RMCO2 / RMD
     ENDIF
     ! ELSE : keep value from .def
  ENDIF

  solaire_etat0      = tab_cntrl(4)
  tab_cntrl(5)=iflag_con
  tab_cntrl(6)=nbapp_rad

  IF (iflag_cycle_diurne.GE.1) tab_cntrl( 7) = iflag_cycle_diurne
  IF (soil_model) tab_cntrl( 8) =1.
  IF (liqice_in_radocond) tab_cntrl( 9) =1.
  IF (ok_orodr) tab_cntrl(10) =1.
  IF (ok_orolf) tab_cntrl(11) =1.
  IF (ok_limitvrai) tab_cntrl(12) =1.
  !GG
  tab_cntrl(18) =iflag_seaice
  tab_cntrl(19) =iflag_seaice_alb
  tab_cntrl(20) =iflag_leads
  !GG

  itau_phy = tab_cntrl(15)

  clesphy0(1)=tab_cntrl( 5 )
  clesphy0(2)=tab_cntrl( 6 )
  clesphy0(3)=tab_cntrl( 7 )
  clesphy0(4)=tab_cntrl( 8 )
  clesphy0(5)=tab_cntrl( 9 )
  clesphy0(6)=tab_cntrl( 10 )
  clesphy0(7)=tab_cntrl( 11 )
  clesphy0(8)=tab_cntrl( 12 )
  clesphy0(9)=tab_cntrl( 17 )

  ! set time iteration
   CALL init_iteration(itau_phy)

  ! read latitudes and make a sanity check (because already known from dyn)
  CALL get_field("latitude",lat_startphy)
  DO i=1,klon
    IF (ABS(lat_startphy(i)-latitude_deg(i))>=1) THEN
      WRITE(*,*) "phyetat0: Error! Latitude discrepancy wrt startphy file:",&
                 " i=",i," lat_startphy(i)=",lat_startphy(i),&
                 " latitude_deg(i)=",latitude_deg(i)
      ! This is presumably serious enough to abort run
      CALL abort_physic("phyetat0","discrepancy in latitudes!",1)
    ENDIF
    IF (ABS(lat_startphy(i)-latitude_deg(i))>=0.0001) THEN
      WRITE(*,*) "phyetat0: Warning! Latitude discrepancy wrt startphy file:",&
                 " i=",i," lat_startphy(i)=",lat_startphy(i),&
                 " latitude_deg(i)=",latitude_deg(i)
    ENDIF
  ENDDO

  ! read longitudes and make a sanity check (because already known from dyn)
  CALL get_field("longitude",lon_startphy)
  DO i=1,klon
    IF (ABS(lon_startphy(i)-longitude_deg(i))>=1) THEN
      IF (ABS(360-ABS(lon_startphy(i)-longitude_deg(i)))>=1) THEN
        WRITE(*,*) "phyetat0: Error! Longitude discrepancy wrt startphy file:",&
                   " i=",i," lon_startphy(i)=",lon_startphy(i),&
                   " longitude_deg(i)=",longitude_deg(i)
        ! This is presumably serious enough to abort run
        CALL abort_physic("phyetat0","discrepancy in longitudes!",1)
      ENDIF
    ENDIF
    IF (ABS(lon_startphy(i)-longitude_deg(i))>=1) THEN
      IF (ABS(360-ABS(lon_startphy(i)-longitude_deg(i))) > 0.0001) THEN
        WRITE(*,*) "phyetat0: Warning! Longitude discrepancy wrt startphy file:",&
                   " i=",i," lon_startphy(i)=",lon_startphy(i),&
                   " longitude_deg(i)=",longitude_deg(i)
      ENDIF 
    ENDIF
  ENDDO

  ! Lecture du masque terre mer

  CALL get_field("masque", zmasq, found)
  IF (.NOT. found) THEN
     PRINT*, 'phyetat0: Le champ <masque> est absent'
     PRINT *, 'fichier startphy non compatible avec phyetat0'
  ENDIF

  ! Lecture des fractions pour chaque sous-surface

  ! initialisation des sous-surfaces

  pctsrf = 0.

  ! fraction de terre

  CALL get_field("FTER", pctsrf(:, is_ter), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FTER> est absent'

  ! fraction de glace de terre

  CALL get_field("FLIC", pctsrf(:, is_lic), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FLIC> est absent'

  ! fraction d'ocean

  CALL get_field("FOCE", pctsrf(:, is_oce), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FOCE> est absent'

  ! fraction glace de mer

  CALL get_field("FSIC", pctsrf(:, is_sic), found)
  IF (.NOT. found) PRINT*, 'phyetat0: Le champ <FSIC> est absent'

  !  Verification de l'adequation entre le masque et les sous-surfaces

  fractint( 1 : klon) = pctsrf(1 : klon, is_ter)  &
       + pctsrf(1 : klon, is_lic)
  DO i = 1 , klon
     IF ( abs(fractint(i) - zmasq(i) ) .GT. EPSFRA ) THEN
        WRITE(*, *) 'phyetat0: attention fraction terre pas ',  &
             'coherente ', i, zmasq(i), pctsrf(i, is_ter) &
             , pctsrf(i, is_lic)
        WRITE(*, *) 'Je force la coherence zmasq=fractint'
        zmasq(i) = fractint(i)
     ENDIF
  ENDDO
  fractint (1 : klon) =  pctsrf(1 : klon, is_oce)  &
       + pctsrf(1 : klon, is_sic)
  DO i = 1 , klon
     IF ( abs( fractint(i) - (1. - zmasq(i))) .GT. EPSFRA ) THEN
        WRITE(*, *) 'phyetat0 attention fraction ocean pas ',  &
             'coherente ', i, zmasq(i) , pctsrf(i, is_oce) &
             , pctsrf(i, is_sic)
        WRITE(*, *) 'Je force la coherence zmasq=1.-fractint'
        zmasq(i) = 1. - fractint(i)
     ENDIF
  ENDDO

!===================================================================
! Lecture des temperatures du sol:
!===================================================================

  found=phyetat0_get(ftsol(:,1),"TS","Surface temperature",283.)
  IF (found) THEN
     DO nsrf=2,nbsrf
        ftsol(:,nsrf)=ftsol(:,1)
     ENDDO
  ELSE
     found=phyetat0_srf(ftsol,"TS","Surface temperature",283.)
  ENDIF

!===================================================================
  ! Lecture des albedo difus et direct
!===================================================================

  DO nsrf = 1, nbsrf
     DO isw=1, nsw
        IF (isw.GT.99) THEN
           PRINT*, "Trop de bandes SW"
           call abort_physic("phyetat0", "", 1)
        ENDIF
        WRITE(str2, '(i2.2)') isw
        found=phyetat0_srf(falb_dir(:, isw,:),"A_dir_SW"//str2//"srf","Direct Albedo",0.2)
        found=phyetat0_srf(falb_dif(:, isw,:),"A_dif_SW"//str2//"srf","Direct Albedo",0.2)
     ENDDO
  ENDDO

  found=phyetat0_srf(u10m,"U10M","u a 10m",0.)
  found=phyetat0_srf(v10m,"V10M","v a 10m",0.)

!===================================================================
! Lecture dans le cas iflag_pbl_surface =1
!===================================================================

   if ( iflag_physiq <= 1 ) then
!===================================================================
  ! Lecture des temperatures du sol profond:
!===================================================================

   DO isoil=1, nsoilmx
        IF (isoil.GT.99) THEN
           PRINT*, "Trop de couches "
           call abort_physic("phyetat0", "", 1)
        ENDIF
        WRITE(str2,'(i2.2)') isoil
        found=phyetat0_srf(tsoil(:, isoil,:),"Tsoil"//str2//"srf","Temp soil",0.)
        IF (.NOT. found) THEN
           PRINT*, "phyetat0: Le champ <Tsoil"//str7//"> est absent"
           PRINT*, "          Il prend donc la valeur de surface"
           tsoil(:, isoil, :)=ftsol(:, :)
        ENDIF
   ENDDO

!=======================================================================
! Lecture precipitation/evaporation
!=======================================================================

  found=phyetat0_srf(qsurf,"QS","Near surface hmidity",0.)
  found=phyetat0_get(qsol,"QSOL","Surface hmidity / bucket",0.)
  found=phyetat0_srf(snow,"SNOW","Surface snow",0.)
  found=phyetat0_srf(fevap,"EVAP","evaporation",0.)
  found=phyetat0_get(snow_fall,"snow_f","snow fall",0.)
  found=phyetat0_get(rain_fall,"rain_f","rain fall",0.)
  IF (ok_bs) THEN
     found=phyetat0_get(bs_fall,"bs_f","blowing snow fall",0.)
  ELSE
     bs_fall(:)=0.
  ENDIF
!=======================================================================
! Radiation
!=======================================================================

  found=phyetat0_get(solsw,"solsw","net SW radiation surf",0.)
  found=phyetat0_get(solswfdiff,"solswfdiff","fraction of SW radiation surf that is diffuse",1.)
  found=phyetat0_get(sollw,"sollw","net LW radiation surf",0.)
  found=phyetat0_get(sollwdown,"sollwdown","down LW radiation surf",0.)
  IF (.NOT. found) THEN
     sollwdown(:) = 0. ;  zts(:)=0.
     DO nsrf=1,nbsrf
        zts(:)=zts(:)+ftsol(:,nsrf)*pctsrf(:,nsrf)
     ENDDO
     sollwdown(:)=sollw(:)+RSIGMA*zts(:)**4
  ENDIF

  found=phyetat0_get(radsol,"RADS","Solar radiation",0.)
  found=phyetat0_get(fder,"fder","Flux derivative",0.) 


  ! Lecture de la longueur de rugosite 
  found=phyetat0_srf(z0m,"RUG","Z0m ancien",0.001)
  IF (found) THEN
     z0h(:,1:nbsrf)=z0m(:,1:nbsrf)
  ELSE
     found=phyetat0_srf(z0m,"Z0m","Roughness length, momentum ",0.001)
     found=phyetat0_srf(z0h,"Z0h","Roughness length, enthalpy ",0.001)
  ENDIF
!FC
  IF (ifl_pbltree>0) THEN
!CALL get_field("FTER", pctsrf(:, is_ter), found)
    treedrg(:,1:klev,1:nbsrf)= 0.0
    CALL get_field("treedrg_ter", drg_ter(:,:), found)
!  found=phyetat0_srf(treedrg,"treedrg","drag from vegetation" , 0.)
    !lecture du profile de freinage des arbres
    IF (.not. found ) THEN
      treedrg(:,1:klev,1:nbsrf)= 0.0
    ELSE
      treedrg(:,1:klev,is_ter)= drg_ter(:,:)
!     found=phyetat0_get(treedrg,"treedrg","freinage arbres",0.)
    ENDIF
  ELSE
    ! initialize treedrg(), because it will be written in restartphy.nc
    treedrg(:,:,:) = 0.0 
  ENDIF

  IF (iflag_hetero_surf .GT. 0) THEN
    found=phyetat0_srf(frac_tersrf,"frac_tersrf","fraction of continental sub-surfaces",0.)
    found=phyetat0_srf(z0m_tersrf,"z0m_tersrf","roughness length for momentum of continental sub-surfaces",0.)
    found=phyetat0_srf(ratio_z0m_z0h_tersrf,"ratio_z0m_z0h_tersrf","ratio of heat to momentum roughness length of continental sub-surfaces",0.)
    found=phyetat0_srf(albedo_tersrf,"albedo_tersrf","albedo of continental sub-surfaces",0.)
    found=phyetat0_srf(beta_tersrf,"beta_tersrf","evapotranspiration coef of continental sub-surfaces",0.)
    found=phyetat0_srf(inertie_tersrf,"inertie_tersrf","soil thermal inertia of continental sub-surfaces",0.)
    found=phyetat0_srf(hcond_tersrf,"hcond_tersrf","heat conductivity of continental sub-surfaces",0.)
    found=phyetat0_srf(tsurfi_tersrf,"tsurfi_tersrf","initial surface temperature of continental sub-surfaces",0.)
    !
    ! Check if the sum of the sub-surface fractions is equal to 1
    DO it=1,klon
      IF (SUM(frac_tersrf(it,:)) .NE. 1.) THEN
        PRINT*, 'SUM(frac_tersrf) = ', SUM(frac_tersrf(it,:))
        CALL abort_physic('conf_phys', 'the sum of fractions of heterogeneous land subsurfaces must be equal &
                          & to 1 for iflag_hetero_surf = 1 and 2',1)
      ENDIF
    ENDDO
    !
    ! Initialisation of surface and soil temperatures (potentially different initial temperatures between sub-surfaces)
    DO iq=1,nbtersrf
      DO it=1,klon
        tsurf_tersrf(it,iq) = tsurfi_tersrf(it,iq)
      ENDDO
    ENDDO
    !
    DO isoil=1, nbtsoildepths 
      IF (isoil.GT.99) THEN
        PRINT*, "Trop de couches "
        CALL abort_physic("phyetat0", "", 1)
      ENDIF
      WRITE(str2,'(i2.2)') isoil
      found=phyetat0_srf(tsoil_depth(:,isoil,:),"tsoil_depth"//str2//"srf","soil depth of continental sub-surfaces",0.)
      found=phyetat0_srf(tsoili_tersrf(:,isoil,:),"Tsoili"//str2//"srf","initial soil temperature of continental sub-surfaces",0.)
      IF (.NOT. found) THEN
        PRINT*, "phyetat0: Le champ <Tsoili"//str2//"> est absent"
        PRINT*, "          Il prend donc la valeur de surface"
        tsoili_tersrf(:, isoil, :) = tsurfi_tersrf(:, :)
      ENDIF
    ENDDO
    !
    tsoil_tersrf = interpol_tsoil(klon, nbtersrf, nsoilmx, nbtsoildepths, alpha_soil_tersrf, period_tersrf, &
                   inertie_tersrf, hcond_tersrf, tsoil_depth, tsurf_tersrf, tsoili_tersrf)
    !
    ! initialise also average surface and soil temperatures
    ftsol(:,is_ter) = average_surf_var(klon, nbtersrf, tsurf_tersrf, frac_tersrf, 'ARI')
    DO k=1, nsoilmx
      tsoil(:,k,is_ter) = average_surf_var(klon, nbtersrf, tsoil_tersrf(:,k,:), frac_tersrf, 'ARI')
    ENDDO
    !
  ENDIF ! iflag_hetero_surf > 0

  endif ! iflag_physiq <= 1

  ! Lecture de l'age de la neige:
  found=phyetat0_srf(agesno,"AGESNO","SNOW AGE",0.001)

  ancien_ok=.true.
  ancien_ok=ancien_ok.AND.phyetat0_get(t_ancien,"TANCIEN","TANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(q_ancien,"QANCIEN","QANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(ql_ancien,"QLANCIEN","QLANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(qs_ancien,"QSANCIEN","QSANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(u_ancien,"UANCIEN","UANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(v_ancien,"VANCIEN","VANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(prw_ancien,"PRWANCIEN","PRWANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(prlw_ancien,"PRLWANCIEN","PRLWANCIEN",0.)
  ancien_ok=ancien_ok.AND.phyetat0_get(prsw_ancien,"PRSWANCIEN","PRSWANCIEN",0.)

  ! cas specifique des variables de la neige soufflee
  IF (ok_bs) THEN
     ancien_ok=ancien_ok.AND.phyetat0_get(qbs_ancien,"QBSANCIEN","QBSANCIEN",0.)
     ancien_ok=ancien_ok.AND.phyetat0_get(prbsw_ancien,"PRBSWANCIEN","PRBSWANCIEN",0.)
  ELSE
     qbs_ancien(:,:)=0.
     prbsw_ancien(:)=0.
  ENDIF
  
  ! cas specifique des variables de la sursaturation par rapport a la glace
  IF ( ok_ice_supersat ) THEN
    ancien_ok=ancien_ok.AND.phyetat0_get(cf_ancien,"CFANCIEN","CFANCIEN",0.)
    ancien_ok=ancien_ok.AND.phyetat0_get(rvc_ancien,"RVCANCIEN","RVCANCIEN",0.)
  ELSE
    cf_ancien(:,:)=0.
    rvc_ancien(:,:)=0.
  ENDIF

  ! Ehouarn: addtional tests to check if t_ancien, q_ancien contain
  !          dummy values (as is the case when generated by ce0l,
  !          or by iniaqua)
  IF ( (maxval(q_ancien).EQ.minval(q_ancien))       .OR. &
       (maxval(ql_ancien).EQ.minval(ql_ancien))     .OR. &
       (maxval(qs_ancien).EQ.minval(qs_ancien))     .OR. &
       (maxval(prw_ancien).EQ.minval(prw_ancien))   .OR. &
       (maxval(prlw_ancien).EQ.minval(prlw_ancien)) .OR. &
       (maxval(prsw_ancien).EQ.minval(prsw_ancien)) .OR. &
       (maxval(t_ancien).EQ.minval(t_ancien)) ) THEN
    ancien_ok=.false.
  ENDIF

  IF (ok_bs) THEN
    IF ( (maxval(qbs_ancien).EQ.minval(qbs_ancien))       .OR. &
         (maxval(prbsw_ancien).EQ.minval(prbsw_ancien)) ) THEN
       ancien_ok=.false.
    ENDIF
  ENDIF

  IF ( ok_ice_supersat ) THEN
    IF ( (maxval(cf_ancien).EQ.minval(cf_ancien))     .OR. &
         (maxval(rvc_ancien).EQ.minval(rvc_ancien)) ) THEN
       ancien_ok=.false.
     ENDIF
  ENDIF


  found=phyetat0_get(clwcon,"CLWCON","CLWCON",0.)
  found=phyetat0_get(rnebcon,"RNEBCON","RNEBCON",0.)
  found=phyetat0_get(ratqs,"RATQS","RATQS",0.)

  found=phyetat0_get(run_off_lic_0,"RUNOFFLIC0","RUNOFFLIC0",0.)

!==================================
!  TKE
!==================================
!
  ! cas specifique de l'advection de TKE
  IF (ok_advtke) THEN
       ancien_ok=ancien_ok.AND.phyetat0_get(tke_ancien,"TKEANCIEN","TKEANCIEN",0.)
  ELSE
    tke_ancien(:,:)=0.
  ENDIF

  IF (ok_advtke) THEN
    IF ( (maxval(tke_ancien).EQ.minval(tke_ancien))) THEN
       ancien_ok=.false.
    ENDIF
  ENDIF

  IF ((iflag_pbl>1)) then
     found=phyetat0_srf(pbl_tke,"TKE","Turb. Kinetic. Energ. ",1.e-8)
  ENDIF

  IF (iflag_pbl>1 .AND. iflag_wake>=1  .AND. iflag_pbl_split >=1 ) then
    found=phyetat0_srf(wake_delta_pbl_tke,"DELTATKE","Del TKE wk/env",0.)
!!    found=phyetat0_srf(delta_tsurf,"DELTA_TSURF","Delta Ts wk/env ",0.)
    found=phyetat0_srf(delta_tsurf,"DELTATS","Delta Ts wk/env ",0.)
!!    found=phyetat0_srf(beta_aridity,"BETA_S","Aridity factor ",1.)
    found=phyetat0_srf(beta_aridity,"BETAS","Aridity factor ",1.)
  ENDIF   !(iflag_pbl>1 .AND. iflag_wake>=1 .AND. iflag_pbl_split >=1 )

!==================================
!  thermiques, poches, convection
!==================================

! Emanuel
  found=phyetat0_get(sig1,"sig1","sig1",0.)
  found=phyetat0_get(w01,"w01","w01",0.)

! Wake
  found=phyetat0_get(wake_deltat,"WAKE_DELTAT","Delta T wake/env",0.)
  found=phyetat0_get(wake_deltaq,"WAKE_DELTAQ","Delta hum. wake/env",0.)
  found=phyetat0_get(wake_s,"WAKE_S","Wake frac. area",0.)
  found=phyetat0_get(awake_s,"AWAKE_S","Active Wake frac. area",0.)
!jyg<
!  Set wake_dens to -1000. when there is no restart so that the actual
!  initialization is made in calwake.
!!  found=phyetat0_get(1,wake_dens,"WAKE_DENS","Wake num. /unit area",0.)
  found=phyetat0_get(wake_dens,"WAKE_DENS","Wake num. /unit area",-1000.)
  found=phyetat0_get(awake_dens,"AWAKE_DENS","Active Wake num. /unit area",0.)
  found=phyetat0_get(cv_gen,"CV_GEN","CB birth rate",0.)
!>jyg
  found=phyetat0_get(wake_cstar,"WAKE_CSTAR","WAKE_CSTAR",0.)
  found=phyetat0_get(wake_pe,"WAKE_PE","WAKE_PE",0.)
  found=phyetat0_get(wake_fip,"WAKE_FIP","WAKE_FIP",0.)

! Thermiques
  found=phyetat0_get(zmax0,"ZMAX0","ZMAX0",40.)
  found=phyetat0_get(f0,"F0","F0",1.e-5)
  found=phyetat0_get(fm_therm,"FM_THERM","Thermals mass flux",0.)
  found=phyetat0_get(entr_therm,"ENTR_THERM","Thermals Entrain.",0.)
  found=phyetat0_get(detr_therm,"DETR_THERM","Thermals Detrain.",0.)

! ALE/ALP
  found=phyetat0_get(ale_bl,"ALE_BL","ALE BL",0.)
  found=phyetat0_get(ale_bl_trig,"ALE_BL_TRIG","ALE BL_TRIG",0.)
  found=phyetat0_get(alp_bl,"ALP_BL","ALP BL",0.)
  found=phyetat0_get(ale_wake,"ALE_WAKE","ALE_WAKE",0.)
  found=phyetat0_get(ale_bl_stat,"ALE_BL_STAT","ALE_BL_STAT",0.)

! fisrtilp/Clouds 0.002 could be ratqsbas. But can stay like this as well
  found=phyetat0_get(ratqs_inter_,"RATQS_INTER","Relative width of the lsc sugrid scale water",0.002)

!===========================================
  ! Read and send field trs to traclmdz
!===========================================

!--OB now this is for co2i - ThL: and therefore also for inco
  IF (ANY(type_trac == ['co2i','inco'])) THEN
     IF (carbon_cycle_cpl) THEN
        ALLOCATE(co2_send(klon), stat=ierr)
        IF (ierr /= 0) CALL abort_physic('phyetat0', 'pb allocation co2_send', 1)
        found=phyetat0_get(co2_send,"co2_send","co2 send",co2_ppm0)
     ENDIF
  ELSE IF (type_trac == 'lmdz') THEN
     it = 0
     DO iq = 1, nqtot
        IF(.NOT.tracers(iq)%isInPhysics) CYCLE
        it = it+1
        tname = tracers(iq)%name
        t(1) = 'trs_'//TRIM(tname); t(2) = 'trs_'//TRIM(new2oldH2O(tname))
        found = phyetat0_get(trs(:,it), t(:), "Surf trac"//TRIM(tname), 0.)
     END DO
     CALL traclmdz_from_restart(trs)
  ENDIF


!===========================================
!  ondes de gravite / relief
!===========================================

!  ondes de gravite non orographiques
  IF (ok_gwd_rando) found = &
       phyetat0_get(du_gwd_rando,"du_gwd_rando","du_gwd_rando",0.)
  IF (.NOT. ok_hines .AND. ok_gwd_rando) found &
       = phyetat0_get(du_gwd_front,"du_gwd_front","du_gwd_front",0.)

!  prise en compte du relief sous-maille
  found=phyetat0_get(zmea,"ZMEA","sub grid orography",0.)
  found=phyetat0_get(zstd,"ZSTD","sub grid orography",0.)
  found=phyetat0_get(zsig,"ZSIG","sub grid orography",0.)
  found=phyetat0_get(zgam,"ZGAM","sub grid orography",0.)
  found=phyetat0_get(zthe,"ZTHE","sub grid orography",0.)
  found=phyetat0_get(zpic,"ZPIC","sub grid orography",0.)
  found=phyetat0_get(zval,"ZVAL","sub grid orography",0.)
  found=phyetat0_get(zmea,"ZMEA","sub grid orography",0.)
  found=phyetat0_get(rugoro,"RUGSREL","sub grid orography",0.)

!===========================================
! Initialize ocean
!===========================================

  IF ( type_ocean == 'slab' ) THEN
      CALL ocean_slab_init(phys_tstep, pctsrf)
      IF (nslay.EQ.1) THEN
        found=phyetat0_get(tslab,["tslab01","tslab  "],"tslab",0.)
      ELSE
          DO i=1,nslay
            WRITE(str2,'(i2.2)') i
            found=phyetat0_get(tslab(:,i),"tslab"//str2,"tslab",0.)  
          ENDDO
      ENDIF
      IF (.NOT. found) THEN 
          PRINT*, "phyetat0: Le champ <tslab> est absent"
          PRINT*, "Initialisation a tsol_oce"
          DO i=1,nslay
              tslab(:,i)=MAX(ftsol(:,is_oce),271.35)
          ENDDO
      ENDIF 

      ! Sea ice variables
      IF (version_ocean == 'sicINT') THEN
          found=phyetat0_get(tice_slab,"slab_tice","slab_tice",0.)
  !GG        found=phyetat0_get(tice,"slab_tice","slab_tice",0.)
          IF (.NOT. found) THEN 
  !GG            PRINT*, "phyetat0: Le champ <tice> est absent"
              PRINT*, "phyetat0: Le champ <tice_slab> est absent"
              PRINT*, "Initialisation a tsol_sic"
  !GG                tice(:)=ftsol(:,is_sic)
                  tice_slab(:)=ftsol(:,is_sic)
          ENDIF 
          found=phyetat0_get(seaice,"seaice","seaice",0.)
          IF (.NOT. found) THEN
              PRINT*, "phyetat0: Le champ <seaice> est absent"
              PRINT*, "Initialisation a 0/1m suivant fraction glace"
              seaice(:)=0.
              WHERE (pctsrf(:,is_sic).GT.EPSFRA)
                  seaice=917.
              ENDWHERE
          ENDIF
      ENDIF !sea ice INT
  ENDIF ! Slab        

  if (activate_ocean_skin >= 1) then
     if (activate_ocean_skin == 2 .and. type_ocean == 'couple') then
        found = phyetat0_get(delta_sal, "delta_sal", &
             "ocean-air interface salinity minus bulk salinity", 0.)
        found = phyetat0_get(delta_sst, "delta_SST", &
             "ocean-air interface temperature minus bulk SST", 0.)
        found = phyetat0_get(dter, "dter", &
             "ocean-air interface temperature minus subskin temperature", 0.)
        found = phyetat0_get(dser, "dser", &
             "ocean-air interface salinity minus subskin salinity", 0.)
        found = phyetat0_get(dt_ds, "dt_ds", "(tks / tkt) * dTer", 0.)

        where (pctsrf(:, is_oce) == 0.)
           delta_sst = missing_val
           delta_sal = missing_val
           dter = missing_val
           dser = missing_val
           dt_ds = missing_val
        end where
     end if
     
     found = phyetat0_get(ds_ns, "dS_ns", "delta salinity near surface", 0.)
     found = phyetat0_get(dt_ns, "dT_ns", "delta temperature near surface", &
          0.)

     where (pctsrf(:, is_oce) == 0.)
        ds_ns = missing_val
        dt_ns = missing_val
        delta_sst = missing_val
        delta_sal = missing_val
     end where
  end if

  !GG
  ! Sea ice
  !IF (iflag_seaice == 2) THEN

  found=phyetat0_get(hice,"hice","Ice thickness",0.)
  IF (.NOT. found) THEN
       PRINT*, "phyetat0: Le champ <hice> est absent"
       PRINT*, "Initialisation a hice=1m "
       hice(:)=1.0
  END IF
  found=phyetat0_get(tice,"tice","Sea Ice temperature",0.)
  IF (.NOT. found) THEN
       PRINT*, "phyetat0: Le champ <tice> est absent"
       PRINT*, "Initialisation a tsol_sic"
       tice(:)=ftsol(:,is_sic)
  END IF
  found=phyetat0_get(bilg_cumul,"bilg_cumul","Flux conductivite + transmit sea-ice",0.)
  IF (.NOT. found) THEN
       PRINT*, "phyetat0: Le champ <bilg_cumul> est absent"
       PRINT*, "Initialisation a zero"
       bilg_cumul(:)=0.0
  END IF

  !END IF
  !GG
  ! on ferme le fichier
  CALL close_startphy

  ! Initialize module pbl_surface_mod 

  if ( iflag_physiq <= 1 ) then
  !GG CALL pbl_surface_init(fder, snow, qsurf, tsoil)
  CALL pbl_surface_init(fder, snow, qsurf, tsoil, hice, tice, bilg_cumul)
  !GG
  endif

  ! --------- COUPLING SECTION --------
  ! Initialize coupling module
  IF ( type_ocean == 'couple' .OR. lk_pycpl ) CALL cpl_def_domain()

  ! Initialize module ocean_cpl_mod for the case of coupled ocean
  IF ( type_ocean == 'couple' ) THEN
     CALL ocean_cpl_init(phys_tstep, longitude_deg, latitude_deg)
  ENDIF

  ! Initialize python communication module
  IF ( lk_pycpl ) THEN
     CALL pyfld_alloc() 
     CALL init_python_coupling() 
  ENDIF

  ! Finalize coupling definition
  IF ( type_ocean == 'couple' .OR. lk_pycpl ) CALL cpl_enddef()
  IF ( type_ocean == 'couple' ) CALL cpl_inca()
  ! -----------------------------------
  
  !  CALL init_iophy_new(latitude_deg, longitude_deg)

  ! Initilialize module fonte_neige_mod      
  CALL fonte_neige_init(run_off_lic_0)

END SUBROUTINE phyetat0

END MODULE phyetat0_mod

