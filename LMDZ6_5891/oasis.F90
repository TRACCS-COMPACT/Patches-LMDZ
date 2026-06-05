MODULE oasis
!
! This module contains subroutines for initialization, sending and receiving 
! towards the coupler OASIS3. It also contains some parameters for the coupling. 
!
! This module should always be compiled. With the coupler OASIS3 available the cpp key
! CPP_COUPLE should be set and the entier of this file will then be compiled. 
! In a forced mode CPP_COUPLE should not be defined and the compilation ends before 
! the CONTAINS, without compiling the subroutines.
!
  USE dimphy 
  USE mod_phys_lmdz_para
  USE write_field_phy

#ifdef CPP_COUPLE
! Use of Oasis-MCT coupler 
#if defined CPP_OMCT
  USE mod_prism
! Use of Oasis3 coupler 
#else 
  USE mod_prism_proto
  USE mod_prism_def_partition_proto
  USE mod_prism_get_proto
  USE mod_prism_put_proto
#endif
#endif
  
  IMPLICIT NONE
 
  ! Ids for oasis configuration
  INTEGER   :: il_part_id

  INTEGER, PUBLIC, PARAMETER ::   midmax=2    ! Maximum number of identified modules   
  INTEGER, PUBLIC, PARAMETER ::   midcpl=1    ! module ID #1 : ocean coupling
  INTEGER, PUBLIC, PARAMETER ::   midpycpl=2  ! module ID #2 : Eophis-equipped external python scripts

  TYPE, PUBLIC ::   FLD_CPL            ! Type for coupling field information
     CHARACTER(len = 8) ::   name      ! Name of the coupling field   
     LOGICAL            ::   action    ! To be exchanged or not
     INTEGER            ::   nid       ! Id of the field
     INTEGER            ::   nlvl      ! Number of level to exchange
  END TYPE FLD_CPL

  TYPE, PUBLIC :: MOD_FLD_CPL             !: Type to sort coupling field information between calling modules
     TYPE(FLD_CPL), DIMENSION(:), ALLOCATABLE,  PUBLIC :: fld
  END TYPE MOD_FLD_CPL

  TYPE(MOD_FLD_CPL), DIMENSION(midmax), PUBLIC ::  inforecv, infosend   !: Coupling field informations
!$OMP THREADPRIVATE(infosend,inforecv)

#ifdef CPP_COUPLE

CONTAINS

  SUBROUTINE cpl_def_domain

     USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, grid_type, unstructured, regular_lonlat
     USE print_control_mod, ONLY: lunout
     USE geometry_mod, ONLY: ind_cell_glo

     ! -----------------------------------------------------------
     ! I/O
     ! Local variables
     INTEGER                            :: ierror, il_commlocal, comp_id
     INTEGER, ALLOCATABLE               :: ig_paral(:)
     INTEGER                            :: jf
     CHARACTER (len = 6)                :: clmodnam
     CHARACTER (len = 20)               :: modname = 'cpl_env_init'
     CHARACTER (len = 80)               :: abort_message 
     INTEGER, DIMENSION(klon_mpi)       :: ind_cell_glo_mpi
     ! -----------------------------------------------------------
     
     ! Define the model name
     IF (grid_type==unstructured) THEN
        clmodnam = 'icosa'                 ! as in $NBMODEL in Cpl/Nam/namcouple.tmp
     ELSE IF (grid_type==regular_lonlat) THEN
        clmodnam = 'LMDZ'                  ! as in $NBMODEL in Cpl/Nam/namcouple.tmp
     ELSE
        abort_message='Pb : type of grid unknown'
        CALL abort_physic(modname,abort_message,1)
     ENDIF
     
     !************************************************************************************
     !  psmile initialisation if not done by mpi module
     !************************************************************************************
     IF (is_sequential) THEN
        CALL prism_init_comp_proto (comp_id, clmodnam, ierror)
       
        IF (ierror .NE. PRISM_Ok) THEN
           abort_message=' Probleme init dans prism_init_comp '
           CALL abort_physic(modname,abort_message,1)
        ELSE
           WRITE(lunout,*) 'inicma : init psmile ok '
        ENDIF
     ENDIF

     CALL prism_get_localcomm_proto (il_commlocal, ierror)
  
     !************************************************************************************
     ! Gather global index to be used for oasis decomposition
     !************************************************************************************
     CALL gather_omp(ind_cell_glo,ind_cell_glo_mpi)

     !************************************************************************************
     ! Domain decomposition
     !************************************************************************************
     IF (grid_type==unstructured) THEN

        ALLOCATE( ig_paral(klon_mpi_para_nb(mpi_rank) + 2) ) 

        ig_paral(1) = 4                                      ! points partition for //
        ig_paral(2) = klon_mpi_para_nb(mpi_rank)             ! nb of local cells

        DO jf=1, klon_mpi_para_nb(mpi_rank)
           ig_paral(2+jf) = ind_cell_glo_mpi(jf)
        ENDDO

     ELSE IF (grid_type==regular_lonlat) THEN

        ALLOCATE( ig_paral(3) )

        ig_paral(1) = 1                            ! apple partition for //
        ig_paral(2) = (jj_begin-1)*nbp_lon+ii_begin-1  ! offset
        ig_paral(3) = (jj_end*nbp_lon+ii_end) - (jj_begin*nbp_lon+ii_begin) + 1

        IF (mpi_rank==mpi_size-1) ig_paral(3)=ig_paral(3)+nbp_lon-1
     ELSE
        abort_message='Pb : type of grid unknown'
        CALL abort_physic(modname,abort_message,1)
     ENDIF

     WRITE(lunout,*) mpi_rank,'ig_paral--->',ig_paral(2),ig_paral(3)
    
     ierror=PRISM_Ok
     CALL prism_def_partition_proto (il_part_id, ig_paral, ierror)

     IF (ierror .NE. PRISM_Ok) THEN
        abort_message=' Probleme dans prism_def_partition '
        CALL abort_physic(modname,abort_message,1)
     ELSE
        WRITE(lunout,*) 'inicma : decomposition domaine psmile ok '
     ENDIF

  END SUBROUTINE cpl_def_domain


  SUBROUTINE cpl_vardef( kmod )

     USE print_control_mod, ONLY: lunout
     USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat

     ! -----------------------------------------------------------
     ! I/O
     INTEGER, INTENT(IN)                :: kmod ! calling module ID
     ! Local variables
     INTEGER                            :: ierror, jf
     INTEGER                            :: il_var_type
     INTEGER, DIMENSION(4)              :: il_var_actual_shape
     CHARACTER (len = 20)               :: modname = 'cpl_vardef'
     CHARACTER (len = 80)               :: abort_message 
     ! -----------------------------------------------------------
     
!$OMP MASTER
     il_var_actual_shape(1) = 1                  ! min of 1st dimension (always 1)
     il_var_actual_shape(2) = nbp_lon            ! max of 1st dimension
     il_var_actual_shape(3) = 1                  ! min of 2nd dimension (always 1)
     il_var_actual_shape(4) = nbp_lat            ! max of 2nd dimension
   
     il_var_type = PRISM_Real
     
     !************************************************************************************
     ! Fields to receive
     ! Loop over all possible variables
     !************************************************************************************

     DO jf=1, SIZE(inforecv(kmod)%fld)
        IF (inforecv(kmod)%fld(jf)%action) THEN
           CALL prism_def_var_proto(inforecv(kmod)%fld(jf)%nid, inforecv(kmod)%fld(jf)%name, il_part_id, &
                (/2,inforecv(kmod)%fld(jf)%nlvl/), PRISM_In, il_var_actual_shape, il_var_type, &
                ierror)
           IF (ierror .NE. PRISM_Ok) THEN
              WRITE(lunout,*) 'cpl_vardef : Problem with prism_def_var_proto for field : ',&
                   inforecv(kmod)%fld(jf)%name
              abort_message=' Problem in call to prism_def_var_proto for fields to receive'
              CALL abort_physic(modname,abort_message,1)
           ENDIF
        ENDIF
     END DO

     !************************************************************************************
     ! Fields to send
     ! Loop over all possible variables
     !************************************************************************************
     DO jf=1,SIZE(infosend(kmod)%fld)
        IF (infosend(kmod)%fld(jf)%action) THEN
           CALL prism_def_var_proto(infosend(kmod)%fld(jf)%nid, infosend(kmod)%fld(jf)%name, il_part_id, &
                (/2,infosend(kmod)%fld(jf)%nlvl/), PRISM_Out, il_var_actual_shape, il_var_type, &
                ierror)
           IF (ierror .NE. PRISM_Ok) THEN
              WRITE(lunout,*) 'clp_vardef : Problem with prism_def_var_proto for field : ',&
                   infosend(kmod)%fld(jf)%name
              abort_message=' Problem in call to prism_def_var_proto for fields to send'
              CALL abort_physic(modname,abort_message,1)
           ENDIF
        ENDIF
     END DO
!$OMP END MASTER
 
  END SUBROUTINE cpl_vardef


  SUBROUTINE cpl_enddef

     USE lmdz_xios
     USE print_control_mod, ONLY: lunout

     ! -----------------------------------------------------------
     ! I/O
     ! Local variables
     INTEGER                       :: ierror
     CHARACTER (len = 20)          :: modname = 'cpl_enddef'
     CHARACTER (len = 80)          :: abort_message 
     ! -----------------------------------------------------------

     !************************************************************************************
     ! End definition
     !************************************************************************************
!$OMP MASTER
     IF (using_xios) CALL xios_oasis_enddef()

     CALL prism_enddef_proto(ierror)
     IF (ierror .NE. PRISM_Ok) THEN
        abort_message=' Problem in call to prism_endef_proto'
        CALL abort_physic(modname,abort_message,1)
     ELSE
        WRITE(lunout,*) 'inicma : endef psmile ok '
     ENDIF
!$OMP END MASTER
     !
  END SUBROUTINE cpl_enddef


  SUBROUTINE cpl_snd( kmod, kid, kstep, pdata )
     
     USE print_control_mod, ONLY: lunout
     
     !!----------------------------------------------------------------------
     ! I/O 
     INTEGER                   , INTENT(in   ) ::   kmod      ! calling module ID
     INTEGER                   , INTENT(in   ) ::   kid       ! variable index in the array
     INTEGER                   , INTENT(in   ) ::   kstep     ! ocean time-step in seconds
     REAL, DIMENSION(:,:)      , INTENT(in   ) ::   pdata
     ! local variables
     INTEGER                                   ::   ierror     ! OASIS3 info argument
     INTEGER                                   ::   jc,jm     ! local loop index
     LOGICAL                                   ::   ll3D      ! flag for 3D coupling
     CHARACTER (len = 80)                      ::   abort_message 
     CHARACTER (len = 20)                      ::   modname = 'cpl_snd'
     !!--------------------------------------------------------------------
     !
     ! Default values
     ll3D = .FALSE.
     IF( infosend(kmod)%fld(kid)%nlvl > 1 ) ll3D = .TRUE.
     !
     ! snd data to OASIS3
     IF( .NOT. ll3D ) THEN   ! send 2D or 3D fields
        CALL prism_put_proto( infosend(kmod)%fld(kid)%nid, kstep, pdata(:,1), ierror )
     ELSE
        CALL prism_put_proto( infosend(kmod)%fld(kid)%nid, kstep, pdata(:,1:infosend(kmod)%fld(kid)%nlvl), ierror )
     ENDIF
     !
     ! Check status
     IF (ierror .NE. PRISM_Ok .AND. ierror.NE.PRISM_Sent .AND. ierror.NE.PRISM_ToRest &
         .AND. ierror.NE.PRISM_LocTrans .AND. ierror.NE.PRISM_Output .AND. &
         ierror.NE.PRISM_SentOut .AND. ierror.NE.PRISM_ToRestOut) THEN
         WRITE (lunout,*) 'Error with sending field :', infosend(kmod)%fld(kid)%name, kstep
         abort_message=' Problem in prism_put_proto '
         CALL abort_physic(modname,abort_message,1)
     ENDIF
     !
  END SUBROUTINE cpl_snd


  SUBROUTINE cpl_rcv( kmod, kid, kstep, pdata )
     
     USE print_control_mod, ONLY: lunout
     
     !!--------------------------------------------------------------------
     ! I/O
     INTEGER                     , INTENT(in   )           ::   kmod      ! calling module ID
     INTEGER                     , INTENT(in   )           ::   kid       ! variable index in the array
     INTEGER                     , INTENT(in   )           ::   kstep     ! ocean time-step in seconds
     REAL      , DIMENSION(:,:)  , INTENT(inout)           ::   pdata     ! IN to keep the value if nothing is done
     ! local variables
     INTEGER                                             ::   ierror     ! OASIS3 info argument
     LOGICAL                                             ::   llaction, ll3D
     CHARACTER (len = 80)                                ::   abort_message 
     CHARACTER (len = 20)                                ::   modname = 'cpl_rcv'
     !!--------------------------------------------------------------------
     !
     ! Default values
     ll3D = .FALSE.
     IF( inforecv(kmod)%fld(kid)%nlvl > 1 ) ll3D = .TRUE.
    
     ! receive data from OASIS3
     IF( .NOT. ll3D ) THEN   ! send 2D or 3D fields
        CALL prism_get_proto( inforecv(kmod)%fld(kid)%nid, kstep, pdata(:,1), ierror )
     ELSE
        CALL prism_get_proto( inforecv(kmod)%fld(kid)%nid, kstep, pdata(:,1:inforecv(kmod)%fld(kid)%nlvl), ierror )
     ENDIF
     !
     ! Check status
     IF (ierror .NE. PRISM_Ok .AND. ierror.NE.PRISM_Recvd .AND. &
        ierror.NE.PRISM_FromRest &
        .AND. ierror.NE.PRISM_Input .AND. ierror.NE.PRISM_RecvOut &
        .AND. ierror.NE.PRISM_FromRestOut) THEN
        WRITE (lunout,*)  'Error with receiving field: ', inforecv(kmod)%fld(kid)%name, kstep
        abort_message=' Problem in prism_get_proto '
        CALL abort_physic(modname,abort_message,1)
     ENDIF
     !
  END SUBROUTINE cpl_rcv

#endif

END MODULE oasis
