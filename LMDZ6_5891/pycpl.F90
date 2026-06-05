MODULE pycpl
   !!====================================================================
   !!                       ***  MODULE  pycpl  ***
   !! Python coupling module : interface for communicating with coupled Python scripts
   !! The module makes no assumptions about configuration of the coupling libraries
   !!====================================================================
   !! History :  LMDZ6  ! 2026-05  (A. Barge)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    'key_eophis'    coupled NEMO/Python-scripts via OASIS3-MCT and Eophis
   !!    'key_...'       coupled NEMO/Python-scripts via another method
   !!----------------------------------------------------------------------

   !!------------------------------ MODULE API ----------------------------
   !!   init_python_coupling     : Initialize coupling with Python
   !!   send_to_python           : send fields to external Python model
   !!   receive_from_python      : receive fields from external Python model
   !!   finalize_python_coupling : Free memory
   !!----------------------------------------------------------------------
   USE eophis_def
   USE oasis 
   USE mod_phys_lmdz_mpi_data
   USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lev
   USE phys_state_var_mod, ONLY: phys_tstep
   USE print_control_mod, ONLY: lunout

   IMPLICIT NONE
   PUBLIC

#if defined key_eophis
   LOGICAL, PUBLIC :: lk_pycpl = .TRUE.
#else
   LOGICAL, PUBLIC :: lk_pycpl = .FALSE.
#endif
   INTEGER, PRIVATE :: kstart, kend

   INTERFACE send_to_python
      MODULE PROCEDURE send_to_python_3d, send_to_python_2d
   END INTERFACE send_to_python

   INTERFACE receive_from_python
      MODULE PROCEDURE receive_from_python_3d, receive_from_python_2d
   END INTERFACE receive_from_python

CONTAINS

   SUBROUTINE init_python_coupling()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE init_python_coupling  ***
      !!
      !! ** Purpose :   Initialisation of the python coupling
      !!
      !! ** Method  :   * Read eophis namelist if used
      !!                * Define exchanges
      !!                * Configure coupling layer
      !!----------------------------------------------------------------------
      ! I/O
      ! local variables
      INTEGER :: ios, jpexch
      INTEGER :: jsnd = 1, jrcv = 1
      TYPE(eophis_var), POINTER :: curr_var
      !!----------------------------------------------------------------------
      !
      ! ===============
      !    Initialize
      ! ===============
      !
      IF (is_mpi_root) THEN    ! control print 
         WRITE(lunout,*) 
         WRITE(lunout,*) 'init_python_coupling: Setting Python models'
         WRITE(lunout,*) '~~~~~~~~~~~~~~~~~~~~'
      END IF
      !
#if defined key_eophis
!$OMP MASTER
      !
      IF (is_mpi_root) WRITE(lunout,*) '      Reading Eophis namelist'
      !
      CALL build_eophis_list(COMM_LMDZ_PHY)
      jpexch = count_eophis_var()
      !
      ALLOCATE( infosend(midpycpl)%fld(jpexch), inforecv(midpycpl)%fld(jpexch) )
      !
      ! ========================================= !
      !     Configure meta-array for coupling     !
      ! ========================================= !
      !
      IF( is_mpi_root ) WRITE(lunout,*) '      Configure coupling layer for pycpl module'
      ! default definitions of infosend(midpycpl)%fld and inforecv(midpycpl)%fld
      infosend(midpycpl)%fld(:)%action = .FALSE.  ;  infosend(midpycpl)%fld(:)%name = ''  ;  infosend(midpycpl)%fld(:)%nlvl = 1
      inforecv(midpycpl)%fld(:)%action = .FALSE.  ;  inforecv(midpycpl)%fld(:)%name = ''  ;  inforecv(midpycpl)%fld(:)%nlvl = 1
      !
      CALL first_eophis_var(curr_var)
      DO WHILE (associated(curr_var))
         IF(.NOT.curr_var%in) THEN
            infosend(midpycpl)%fld(jsnd)%name = curr_var%alias
            infosend(midpycpl)%fld(jsnd)%action = .TRUE.
            infosend(midpycpl)%fld(jsnd)%nlvl = curr_var%nlvl
            curr_var%idx = jsnd
            jsnd = jsnd + 1
         ELSE
            inforecv(midpycpl)%fld(jrcv)%name = curr_var%alias
            inforecv(midpycpl)%fld(jrcv)%action = .TRUE.
            inforecv(midpycpl)%fld(jrcv)%nlvl = curr_var%nlvl
            curr_var%idx = jrcv
            jrcv = jrcv + 1
         ENDIF
         CALL eophis_next_var(curr_var)
      END DO
      !
      ! Array bounds
      kstart = ii_begin
      IF (is_south_pole_dyn) THEN
          kend = (jj_end-jj_begin)*nbp_lon + nbp_lon
      ELSE
          kend = (jj_end - jj_begin)*nbp_lon + ii_end
      ENDIF 
      !
      ! ============================== !
      !    Configure coupling layer    !
      ! ============================== !
      CALL cpl_vardef(midpycpl)
#endif
!$OMP END MASTER
      !
   END SUBROUTINE init_python_coupling


   SUBROUTINE send_to_python_3d(varname,to_send,kt)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE send_to_python ***
      !!
      !! ** Purpose :   Proceed coupler sending from coupling definition
      !!
      !! ** Arguments : CHAR varname : name of the field to send
      !!                REAL(:,:,:) to_send  : Array to send
      !!                INT kt : ocean time step
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt             ! ocean time step
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:,:,:), INTENT(in) ::  to_send
      ! local variables
      INTEGER :: isec, ilvl
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,nbp_lev) :: zbuf
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
!$OMP MASTER
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      ! Get Eophis variable
      CALL find_eophis_var(varname,curr_var)
      IF (.NOT.associated(curr_var)) THEN
         CALL abort_physic( 'send_to_python', ' unrecognized variable name '//TRIM(varname) )
      END IF
      !
      ! Coupling layer
      IF (curr_var%in) THEN
         CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
      ELSE
         DO ilvl = 1, infosend(midpycpl)%fld(curr_var%idx)%nlvl
            zbuf(:,ilvl) = RESHAPE(to_send(:,:,ilvl),(/nbp_lon*jj_nb/))
         END DO
         CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
      END IF
!$OMP END MASTER
#endif
      !
   END SUBROUTINE send_to_python_3d


   SUBROUTINE send_to_python_2d(varname,to_send,kt)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE send_to_python ***
      !!
      !! ** Purpose :   Proceed coupler sending from coupling definition
      !!
      !! ** Arguments : CHAR varname : name of the field to send
      !!                REAL(:,:) to_send  : Array to send
      !!                INT kt : ocean time step
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt             ! ocean time step
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:,:), INTENT(in) ::  to_send
      ! local variables
      INTEGER :: isec
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,1) :: zbuf
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
!$OMP MASTER
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      ! Get Eophis variable
      CALL find_eophis_var(varname,curr_var)
      IF (.NOT.associated(curr_var)) THEN
         CALL abort_physic( 'send_to_python' , ' unrecognized variable name '//TRIM(varname) )
      END IF
      !
      ! Coupling layer
      IF (curr_var%in) THEN
         CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
      ELSE
         zbuf(:,1) = RESHAPE(to_send(:,:),(/nbp_lon*jj_nb/))
         CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
      END IF
!$OMP END MASTER
#endif
      !
   END SUBROUTINE send_to_python_2d


   SUBROUTINE receive_from_python_3d(varname,to_rcv,kt)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE receive_from_python  ***
      !!
      !! ** Purpose :   Proceed coupler receiving from coupling definition
      !!
      !! ** Arguments : CHAR varname : name of the field to receive
      !!                REAL(:,:,:) to_rcv : Array in which store received field
      !!                INT kt : ocean time step
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:,:,:), INTENT(inout) ::  to_rcv
      ! local variables
      INTEGER :: isec, ilvl
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,nbp_lev) :: zbuf
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
!$OMP MASTER
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      ! Get Eophis variable
      CALL find_eophis_var(varname,curr_var)
      IF (.NOT.associated(curr_var)) THEN
         CALL abort_physic( 'receive_from_python' , ' unrecognized variable name '//TRIM(varname) )
      END IF
      !
      ! Coupling layer
      IF (.NOT. curr_var%in) THEN
         CALL abort_physic( 'receive_from_python' , ' function called for outcoming variable '//TRIM(varname) )
      ELSE
         ! save value if nothing is done
         DO ilvl = 1, inforecv(midpycpl)%fld(curr_var%idx)%nlvl
            zbuf(:,ilvl) = RESHAPE(to_rcv(:,:,ilvl),(/nbp_lon*jj_nb/))
         END DO
         CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
         DO ilvl = 1, inforecv(midpycpl)%fld(curr_var%idx)%nlvl
            to_rcv(:,:,ilvl) = RESHAPE(zbuf(:,ilvl),(/nbp_lon,jj_nb/))
         END DO
      END IF
!$OMP END MASTER
#endif
      !
   END SUBROUTINE receive_from_python_3d


   SUBROUTINE receive_from_python_2d(varname,to_rcv,kt)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE receive_from_python  ***
      !!
      !! ** Purpose :   Proceed coupler receiving from coupling definition
      !!
      !! ** Arguments : CHAR varname : name of the field to receive
      !!                REAL(:,:) to_rcv : Array in which store received field
      !!                INT kt : ocean time step
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:,:), INTENT(inout) ::  to_rcv
      ! local variables
      INTEGER :: isec
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,1) :: zbuf
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
!$OMP MASTER
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      ! Get Eophis variable
      CALL find_eophis_var(varname,curr_var)
      IF (.NOT.associated(curr_var)) THEN
         CALL abort_physic( 'receive_from_python' , ' unrecognized variable name '//TRIM(varname) )
      END IF
      !
      ! Coupling layer
      IF (.NOT. curr_var%in) THEN
         CALL abort_physic( 'receive_from_python' , ' function called for outcoming variable '//TRIM(varname) )
      ELSE
         ! save value if nothing is done
         zbuf(:,1) = RESHAPE(to_rcv(:,:),(/nbp_lon*jj_nb/))
         CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
         to_rcv(:,:) = RESHAPE(zbuf(:,1),(/nbp_lon,jj_nb/))
      END IF
!$OMP END MASTER
#endif
      !
   END SUBROUTINE receive_from_python_2d


   SUBROUTINE finalize_python_coupling
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE finalize_python_coupling  ***
      !!
      !! ** Purpose :   Free memory used for Python coupling
      !!
      !! ** Method  :   * Deallocate arrays
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
!$OMP MASTER
      DEALLOCATE(infosend(midpycpl)%fld,inforecv(midpycpl)%fld)
      CALL purge_eophis()
!$OMP END MASTER
#endif
      !
   END SUBROUTINE finalize_python_coupling

END MODULE pycpl

