MODULE pycpl
   !!====================================================================
   !!                       ***  MODULE  pycpl  ***
   !! Python coupling module : interface for communicating with coupled Python scripts
   !! The module makes no assumptions about configuration of the coupling libraries
   !!====================================================================
   !! History :  LMDZ6  ! 2026-05  (A. Barge)  Original code
   !!            LMDZ6  ! 2026-09  (A. Barge)  Generic grids: physics and dynamics
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
   !!
   !!   Exchanged fields live on the coupling grid (nbp_lon, jj_nb, nlvl).
   !!   Three families of layouts:
   !!     - coupling grid (nbp_lon, jj_nb, nbp_lev) : no flag
   !!     - physics grid  (klon, klev) / (klon)     : flag pycpl_phys
   !!     - dynamics grids, flattened (ij, llm)     : flags pycpl_phys_dyn_u / pycpl_phys_dyn_v
   !!
   !!   kt argument is the physics time step number, whatever the calling context.
   !!
   !!   This module is a pure transition layer: it only performs grid
   !!   remapping (and the grid completion that the caller cannot provide:
   !!   fake 90S row and duplicated last v row on pycpl_phys_dyn_v, wrap column
   !!   iip1 on the dynamics grids). Values are sent / received AS PROVIDED
   !!----------------------------------------------------------------------
   USE eophis_def
   USE oasis
   USE mod_phys_lmdz_mpi_data
   USE mod_grid_phy_lmdz, ONLY: nbp_lon, nbp_lat, nbp_lev
   USE phys_state_var_mod, ONLY: phys_tstep
   USE print_control_mod, ONLY: lunout
   USE cpl_mod, ONLY: gath2cpl, cpl2gath
   USE mod_phys_lmdz_para, ONLY: bcast_omp
   USE dimphy, ONLY: klon, klev

   IMPLICIT NONE
   PUBLIC

#if defined key_eophis
   LOGICAL, PUBLIC :: lk_pycpl = .TRUE.
#else
   LOGICAL, PUBLIC :: lk_pycpl = .FALSE.
#endif
   INTEGER, PRIVATE :: kstart, kend

   ! Identity index array for physics-to-coupling grid transformation.
   INTEGER, PRIVATE, ALLOCATABLE, SAVE, DIMENSION(:) :: pycpl_unity
!$OMP THREADPRIVATE(pycpl_unity)

   ! Grid selectors
   INTEGER, PUBLIC, PARAMETER :: pycpl_phys  = 1  ! physics grid (klon,klev) / (klon)
   INTEGER, PUBLIC, PARAMETER :: pycpl_phys_dyn_u = 2  ! dynamics u-grid flattened on the physics band
   INTEGER, PUBLIC, PARAMETER :: pycpl_phys_dyn_v = 3  ! dynamics v-grid flattened on the physics band

   INTERFACE send_to_python
      MODULE PROCEDURE send_to_python_3d, send_to_python_2d, &
                       send_to_python_gen_3d, send_to_python_gen_2d
   END INTERFACE send_to_python

   INTERFACE receive_from_python
      MODULE PROCEDURE receive_from_python_3d, receive_from_python_2d, &
                       receive_from_python_gen_3d, receive_from_python_gen_2d
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
      INTEGER :: ios, jpexch, ig
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
      !
      ! Identity index array for physics-to-coupling grid transformation.
      ! every OMP thread allocates and fills its own copy
      ALLOCATE(pycpl_unity(klon))
      DO ig = 1, klon
          pycpl_unity(ig) = ig
      ENDDO
      !
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
!$OMP END MASTER
#endif
      !
   END SUBROUTINE init_python_coupling


   SUBROUTINE send_to_python_3d(varname,to_send,kt)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE send_to_python ***
      !!
      !! ** Purpose :   Proceed coupler sending from coupling definition
      !!
      !! ** Arguments : CHAR varname : name of the field to send
      !!                REAL(:,:,:) to_send  : Array to send on the coupling grid
      !!                INT kt : time step
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
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
      ! Check
      IF (curr_var%in) THEN
         CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
      END IF
      !
      ! coupling layer
      DO ilvl = 1, infosend(midpycpl)%fld(curr_var%idx)%nlvl
          zbuf(:,ilvl) = RESHAPE(to_send(:,:,ilvl),(/nbp_lon*jj_nb/))
      END DO
      CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
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
      !!                REAL(:,:) to_send  : Array to send on the coupling grid
      !!                INT kt : time step
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
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
      ! Check
      IF (curr_var%in) THEN
         CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
      END IF
      !
      ! Coupling layer
      zbuf(:,1) = RESHAPE(to_send(:,:),(/nbp_lon*jj_nb/))
      CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
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
      !!                INT kt : time step
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
      ! Check
      IF (.NOT. curr_var%in) THEN
         CALL abort_physic( 'receive_from_python' , ' function called for outgoing variable '//TRIM(varname) )
      END IF
      !
      ! save value if nothing is done
      DO ilvl = 1, inforecv(midpycpl)%fld(curr_var%idx)%nlvl
         zbuf(:,ilvl) = RESHAPE(to_rcv(:,:,ilvl),(/nbp_lon*jj_nb/))
      END DO
      !
      ! Coupling layer
      CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
      DO ilvl = 1, inforecv(midpycpl)%fld(curr_var%idx)%nlvl
         to_rcv(:,:,ilvl) = RESHAPE(zbuf(:,ilvl),(/nbp_lon,jj_nb/))
      END DO
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
      !!                INT kt : time step
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
      ! Check
      IF (.NOT. curr_var%in) THEN
         CALL abort_physic( 'receive_from_python' , ' function called for outgoing variable '//TRIM(varname) )
      END IF
      !
      ! Save value if nothing is done
      zbuf(:,1) = RESHAPE(to_rcv(:,:),(/nbp_lon*jj_nb/))
      !
      ! Coupling layer
      CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
      to_rcv(:,:) = RESHAPE(zbuf(:,1),(/nbp_lon,jj_nb/))
!$OMP END MASTER
#endif
      !
   END SUBROUTINE receive_from_python_2d


   SUBROUTINE send_to_python_gen_3d(varname,to_send,kt,grid)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE send_to_python_gen_3d  ***
      !!
      !! ** Purpose :   Send a 3D field to the Python coupler from either the
      !!                physics grid or one of the dynamics grids.
      !!
      !! ** Arguments : CHAR varname   : name of the field to send
      !!                REAL(:,:) to_send : 3D field, layout depends on the grid:
      !!                   pycpl_phys     : (klon, klev), physics grid
      !!                   pycpl_phys_dyn_u : flattened (ij, llm) u-grid, indexed by the
      !!                                global flattened dynamics index; must cover
      !!                                the coupling band (distrib_physic buffers do)
      !!                   pycpl_phys_dyn_v : flattened (ij, llm) v-grid, indexed by the
      !!                                global flattened dynamics index; must cover
      !!                                the coupling band (distrib_physic buffers do)
      !!                INT kt        : physics time step (same numbering as physiq itap)
      !!                INT grid      : grid selector
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
      INTEGER, INTENT(in)           ::  grid
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:,:), INTENT(in) ::  to_send
      ! local variables
      INTEGER :: isec, ilvl, nlvl, i, j, l, nrow, ij, ij0, ij_lo, ij_hi, j_lo, j_hi
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,nbp_lev) :: zbuf
      REAL, DIMENSION(nbp_lon,jj_nb,nbp_lev) :: field_3d
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      SELECT CASE (grid)
      !
      CASE (pycpl_phys)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'send_to_python', ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (curr_var%in) THEN
            CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
         END IF
         !
         ! Number of levels to send
         nlvl = infosend(midpycpl)%fld(curr_var%idx)%nlvl
!$OMP END MASTER
         !
         ! Distribute the level count to all threads
         CALL bcast_omp(nlvl)
         !
         ! Gather from physics to coupling grid
         DO ilvl = 1, nlvl
            CALL gath2cpl(to_send(:,ilvl), field_3d(:,:,ilvl), klon, pycpl_unity)
         END DO
         !
!$OMP MASTER
         ! Coupling layer
         zbuf(:,1:nlvl) = RESHAPE(field_3d(:,:,1:nlvl),(/nbp_lon*jj_nb,nlvl/))
         CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,1:nlvl))
!$OMP END MASTER
         !
      CASE (pycpl_phys_dyn_u, pycpl_phys_dyn_v)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'send_to_python', ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (curr_var%in) THEN
            CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
         END IF
         !
         ! Number of levels to send
         nlvl = infosend(midpycpl)%fld(curr_var%idx)%nlvl
         !
         ! Band of the coupling grid flattened on the dynamics grid
         ij_lo = LBOUND(to_send,1)
         ij_hi = UBOUND(to_send,1)
         j_lo  = (ij_lo-1)/(nbp_lon+1) + 1
         j_hi  = ij_hi/(nbp_lon+1)
         !
         ! v-grid has no pole row. Fill last coupling row with a duplicate of the last v row
         IF (grid == pycpl_phys_dyn_v .AND. is_south_pole_dyn) THEN
            nrow = (nbp_lat-1) - j_lo + 1
         ELSE
            nrow = j_hi - j_lo + 1
         END IF
         !
         ! Grid remapping
         zbuf = 0.
         DO j = j_lo, j_hi
            DO i = 1, nbp_lon
               ij  = (j-1)*(nbp_lon+1) + i
               ij0 = ij - ij_lo + 1
               DO l = 1, nlvl
                  zbuf(i+(j-j_lo)*nbp_lon, l) = to_send(ij0,l)
               END DO
            END DO
         END DO
         !
         ! Duplicate last v row
         IF (grid == pycpl_phys_dyn_v .AND. is_south_pole_dyn .AND. nrow >= 1) THEN
            DO l = 1, nlvl
               zbuf((jj_nb-1)*nbp_lon+1:jj_nb*nbp_lon, l) = zbuf((jj_nb-2)*nbp_lon+1:(jj_nb-1)*nbp_lon, l)
            END DO
         END IF
         !
         ! Coupling layer
         CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,1:nlvl))
!$OMP END MASTER
         !
      CASE DEFAULT
         CALL abort_physic( 'send_to_python', ' unsupported grid selector for '//TRIM(varname) )
         !
      END SELECT
#endif
      !
   END SUBROUTINE send_to_python_gen_3d


   SUBROUTINE send_to_python_gen_2d(varname,to_send,kt,grid)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE send_to_python_gen_2d  ***
      !!
      !! ** Purpose :   Send a 2D field to the Python coupler from either the
      !!                physics grid or the dynamics u-grid.
      !!
      !! ** Arguments : CHAR varname : name of the field to send
      !!                REAL(:) to_send : 2D field, layout depends on the grid:
      !!                   pycpl_phys     : (klon), physics grid
      !!                   pycpl_phys_dyn_u : flattened (ij), u-grid, indexed by the
      !!                                global flattened dynamics index; must cover
      !!                                the coupling band (distrib_physic buffers do)
      !!                INT kt          : physics time step (same numbering as physiq itap)
      !!                INT grid        : grid selector
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
      INTEGER, INTENT(in)           ::  grid
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:), INTENT(in) ::  to_send
      ! local variables
      INTEGER :: isec, i, j, nrow, ij, ij0, ij_lo, ij_hi, j_lo, j_hi
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,1) :: zbuf
      REAL, DIMENSION(nbp_lon,jj_nb) :: field_2d
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      SELECT CASE (grid)
      !
      CASE (pycpl_phys)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'send_to_python', ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (curr_var%in) THEN
            CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
         END IF
!$OMP END MASTER
         !
         ! Gather from physics to coupling grid
         CALL gath2cpl(to_send, field_2d, klon, pycpl_unity)
         !
!$OMP MASTER
         ! Coupling layer
         zbuf(:,1) = RESHAPE(field_2d,(/nbp_lon*jj_nb/))
         CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
!$OMP END MASTER
         !
      CASE (pycpl_phys_dyn_u)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'send_to_python', ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (curr_var%in) THEN
            CALL abort_physic( 'send_to_python' , ' function called for incoming variable '//TRIM(varname) )
         END IF
         !
         ! Band of the coupling grid flattened on the dynamics grid
         ij_lo = LBOUND(to_send,1)
         ij_hi = UBOUND(to_send,1)
         j_lo  = (ij_lo-1)/(nbp_lon+1) + 1
         j_hi  = ij_hi/(nbp_lon+1)
         !
         ! Grid remapping
         DO j = j_lo, j_hi
            DO i = 1, nbp_lon
               ij  = (j-1)*(nbp_lon+1) + i
               ij0 = ij - ij_lo + 1
               zbuf(i+(j-j_lo)*nbp_lon, 1) = to_send(ij0)
            END DO
         END DO
         !
         ! Coupling layer
         CALL cpl_snd(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
!$OMP END MASTER
         !
      CASE DEFAULT
         CALL abort_physic( 'send_to_python', ' unsupported grid selector for '//TRIM(varname) )
         !
      END SELECT
#endif
      !
   END SUBROUTINE send_to_python_gen_2d


   SUBROUTINE receive_from_python_gen_3d(varname,to_rcv,kt,grid)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE receive_from_python_gen_3d  ***
      !!
      !! ** Purpose :   Receive a 3D field from the Python coupler on either the
      !!                physics grid or one of the dynamics grids.
      !!
      !! ** Arguments : CHAR varname   : name of the field to receive
      !!                REAL(:,:) to_rcv : 3D field, layout depends on the grid:
      !!                   pycpl_phys     : (klon, klev), physics grid
      !!                   pycpl_phys_dyn_u : flattened (ij, llm) u-grid, indexed by the
      !!                                global flattened dynamics index; must cover
      !!                                the coupling band (distrib_physic buffers do)
      !!                   pycpl_phys_dyn_v : flattened (ij, llm) v-grid, indexed by the
      !!                                global flattened dynamics index; must cover
      !!                                the coupling band (distrib_physic buffers do)
      !!                INT kt        : physics time step (same numbering as physiq itap)
      !!                INT grid      : grid selector
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
      INTEGER, INTENT(in)           ::  grid
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:,:), INTENT(inout) ::  to_rcv
      ! local variables
      INTEGER :: isec, ilvl, nlvl, i, j, l, nrow, ij, ij0, ij_lo, ij_hi, j_lo, j_hi
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,nbp_lev) :: zbuf
      REAL, DIMENSION(nbp_lon,jj_nb,nbp_lev) :: field_3d
      REAL, DIMENSION(klon_mpi) :: gath_buf
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      SELECT CASE (grid)
      !
      CASE (pycpl_phys)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'receive_from_python' , ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (.NOT. curr_var%in) THEN
            CALL abort_physic( 'receive_from_python' , ' function called for outgoing variable '//TRIM(varname) )
         END IF
         !
         ! Number of levels to receive
         nlvl = inforecv(midpycpl)%fld(curr_var%idx)%nlvl
!$OMP END MASTER
         !
         ! Distribute the level count to all threads
         CALL bcast_omp(nlvl)
         !
         ! save value if nothing is done
         DO ilvl = 1, nlvl
            CALL gath2cpl(to_rcv(:,ilvl), field_3d(:,:,ilvl), klon, pycpl_unity)
         END DO
         !
!$OMP MASTER
         ! Coupling layer
         zbuf(:,1:nlvl) = RESHAPE(field_3d(:,:,1:nlvl),(/nbp_lon*jj_nb,nlvl/))
         CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,1:nlvl))
         field_3d(:,:,1:nlvl) = RESHAPE(zbuf(:,1:nlvl),(/nbp_lon,jj_nb,nlvl/))
!$OMP END MASTER
         !
         ! Scatter from coupling to physics grid
         DO ilvl = 1, nlvl
            CALL cpl2gath(field_3d(:,:,ilvl), gath_buf, klon, pycpl_unity)
            to_rcv(:,ilvl) = gath_buf(1:klon)
         END DO
         !
      CASE (pycpl_phys_dyn_u, pycpl_phys_dyn_v)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'receive_from_python' , ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (.NOT. curr_var%in) THEN
            CALL abort_physic( 'receive_from_python' , ' function called for outgoing variable '//TRIM(varname) )
         END IF
         !
         ! Number of levels to receive
         nlvl = inforecv(midpycpl)%fld(curr_var%idx)%nlvl
         !
         ! Band of the coupling grid flattened on the dynamics grid
         ij_lo = LBOUND(to_rcv,1)
         ij_hi = UBOUND(to_rcv,1)
         j_lo  = (ij_lo-1)/(nbp_lon+1) + 1
         j_hi  = ij_hi/(nbp_lon+1)
         !
         ! Drop last fake S row: v-grid has no pole row
         IF (grid == pycpl_phys_dyn_v .AND. is_south_pole_dyn) THEN
            nrow = (nbp_lat-1) - j_lo + 1
         ELSE
            nrow = j_hi - j_lo + 1
         END IF
         !
         ! save value if nothing is done
         zbuf = 0.
         DO j = j_lo, j_hi
            DO i = 1, nbp_lon
               ij  = (j-1)*(nbp_lon+1) + i
               ij0 = ij - ij_lo + 1
               DO l = 1, nlvl
                  zbuf(i+(j-j_lo)*nbp_lon, l) = to_rcv(ij0,l)
               END DO
            END DO
         END DO
         !
         ! Coupling layer
         CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,1:nlvl))
         !
         ! Scatter from coupling to dynamics grid
         DO j = j_lo, j_hi
            DO i = 1, nbp_lon
               ij  = (j-1)*(nbp_lon+1) + i
               ij0 = ij - ij_lo + 1
               DO l = 1, nlvl
                  to_rcv(ij0,l) = zbuf(i+(j-j_lo)*nbp_lon, l)
               END DO
            END DO
            ! column iip1 duplicates column 1 (dynamics periodicity)
            ij = (j-1)*(nbp_lon+1)
            DO l = 1, nlvl
               to_rcv(ij+nbp_lon+1-ij_lo+1, l) = to_rcv(ij+1-ij_lo+1, l)
            END DO
         END DO
!$OMP END MASTER
         !
      CASE DEFAULT
         CALL abort_physic( 'receive_from_python', ' unsupported grid selector for '//TRIM(varname) )
         !
      END SELECT
#endif
      !
   END SUBROUTINE receive_from_python_gen_3d


   SUBROUTINE receive_from_python_gen_2d(varname,to_rcv,kt,grid)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE receive_from_python_gen_2d  ***
      !!
      !! ** Purpose :   Receive a 2D field from the Python coupler on the
      !!                physics grid (klon).
      !!
      !! ** Arguments : CHAR varname : name of the field to receive
      !!                REAL(:) to_rcv : 2D field on the physics grid (klon)
      !!                INT kt          : physics time step (same numbering as physiq itap)
      !!                INT grid        : grid selector (pycpl_phys only)
      !!----------------------------------------------------------------------
      ! I/O
      INTEGER, INTENT(in)           ::  kt
      INTEGER, INTENT(in)           ::  grid
      CHARACTER(len=*), INTENT(in)  :: varname
      REAL, DIMENSION(:), INTENT(inout) ::  to_rcv
      ! local variables
      INTEGER :: isec
      TYPE(eophis_var), POINTER :: curr_var
      REAL, DIMENSION(nbp_lon*jj_nb,1) :: zbuf
      REAL, DIMENSION(nbp_lon,jj_nb) :: field_2d
      REAL, DIMENSION(klon_mpi) :: gath_buf
      !!----------------------------------------------------------------------
      !
#if defined key_eophis
      ! Date of exchange
      isec = ( kt - 1 ) * phys_tstep
      !
      SELECT CASE (grid)
      !
      CASE (pycpl_phys)
         !
!$OMP MASTER
         ! Get Eophis variable
         CALL find_eophis_var(varname,curr_var)
         IF (.NOT.associated(curr_var)) THEN
            CALL abort_physic( 'receive_from_python' , ' unrecognized variable name '//TRIM(varname) )
         END IF
         !
         ! Check
         IF (.NOT. curr_var%in) THEN
            CALL abort_physic( 'receive_from_python' , ' function called for outgoing variable '//TRIM(varname) )
         END IF
!$OMP END MASTER
         !
         ! save value if nothing is done
         CALL gath2cpl(to_rcv, field_2d, klon, pycpl_unity)
         !
!$OMP MASTER
         ! Coupling layer
         zbuf(:,1) = RESHAPE(field_2d,(/nbp_lon*jj_nb/))
         CALL cpl_rcv(midpycpl, curr_var%idx, isec, zbuf(kstart:kend,:))
         field_2d = RESHAPE(zbuf(:,1),(/nbp_lon,jj_nb/))
!$OMP END MASTER
         !
         ! Scatter from coupling to physics grid
         CALL cpl2gath(field_2d, gath_buf, klon, pycpl_unity)
         to_rcv(:) = gath_buf(1:klon)
         !
      CASE DEFAULT
         CALL abort_physic( 'receive_from_python', ' 2D reception is only available on the physics grid' )
         !
      END SELECT
#endif
      !
   END SUBROUTINE receive_from_python_gen_2d


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
      !
      IF (ALLOCATED(pycpl_unity)) DEALLOCATE(pycpl_unity)
#endif
      !
   END SUBROUTINE finalize_python_coupling

END MODULE pycpl
