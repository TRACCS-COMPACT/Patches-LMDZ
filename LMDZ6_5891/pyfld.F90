MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE  pyfld  ***
   !! Python module : fields exchanged with Python scripts stored in core memory
   !!======================================================================
   !! History :  LMDZ6  ! 2026-06  (A. Barge)  Original code
   !!----------------------------------------------------------------------
   USE pycpl
   USE dimensions_mod, ONLY: llm, jjm
   USE mod_phys_lmdz_mpi_data, ONLY: jj_begin, jj_end
   USE mod_grid_phy_lmdz, ONLY: nbp_lon

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!          Seasonal cycle scalars on the coupling grid
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: nn_cosday, nn_sinday

   !!----------------------------------------------------------------------
   !!          NN inputs on the flattened dynamics grids
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: nn_u_dyn, nn_v_dyn
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    :: nn_topo_dyn

   !!----------------------------------------------------------------------
   !!          NN outputss on the flattened dynamics grids
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: py_du_dyn, py_dv_dyn

   ! Flattened dynamics bounds of the coupling band
   INTEGER, PUBLIC, SAVE :: ij_lo_u, ij_hi_u, ij_lo_v, ij_hi_v

CONTAINS

   SUBROUTINE pyfld_alloc()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pyfld_alloc  ***
      !!
      !! ** Purpose :   Initialisation of the Python coupling working arrays
      !!
      !! ** Method  :   * Allocate arrays for Python fields
      !!----------------------------------------------------------------------
      !
 !$OMP MASTER
      ! Allocate arrays
      IF ( lk_pycpl ) THEN
         ! Coupling grid
         ALLOCATE( nn_cosday(nbp_lon,jj_nb,nbp_lev), nn_sinday(nbp_lon,jj_nb,nbp_lev) )
         ! Dynamics grids without halos
         ij_lo_u = (jj_begin-1)*(nbp_lon+1) + 1
         ij_hi_u =  jj_end   *(nbp_lon+1)
         ij_lo_v = ij_lo_u
         ij_hi_v = MIN(jj_end, nbp_lat-1)*(nbp_lon+1)
         ALLOCATE( nn_u_dyn(ij_lo_u:ij_hi_u, llm) )
         ALLOCATE( nn_v_dyn(ij_lo_v:ij_hi_v, llm) )
         ALLOCATE( nn_topo_dyn(ij_lo_u:ij_hi_u) )
         ALLOCATE( py_du_dyn(ij_lo_u:ij_hi_u, llm) )
         ALLOCATE( py_dv_dyn(ij_lo_v:ij_hi_v, llm) )
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_alloc


   SUBROUTINE pyfld_dealloc()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE finalize_python_fields  ***
      !!
      !! ** Purpose :   Free memory used by Python fields
      !!
      !! ** Method  :   * deallocate arrays for Python fields
      !!----------------------------------------------------------------------
      !
 !$OMP MASTER
      ! Free memory
      IF ( lk_pycpl ) THEN
         DEALLOCATE( nn_cosday, nn_sinday )
         IF ( ALLOCATED(nn_u_dyn) ) DEALLOCATE( nn_u_dyn )
         IF ( ALLOCATED(nn_v_dyn) ) DEALLOCATE( nn_v_dyn )
         IF ( ALLOCATED(nn_topo_dyn) ) DEALLOCATE( nn_topo_dyn )
         IF ( ALLOCATED(py_du_dyn) ) DEALLOCATE( py_du_dyn )
         IF ( ALLOCATED(py_dv_dyn) ) DEALLOCATE( py_dv_dyn )
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_dealloc


   SUBROUTINE pycpl_dyn_fields(ucov, vcov, teta, ps, phis, itau, iphysiq)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pycpl_exchange  ***
      !!
      !! ** Purpose :   Perform a full Python coupling exchange of dynamic fields
      !!
      !! ** Arguments : REAL ucov(ijb_u:ije_u, llm)   : zonal covariant wind (caldyn, IN)
      !!                REAL vcov(ijb_v:ije_v, llm)   : meridional covariant wind (caldyn, IN)
      !!                REAL teta(ijb_u:ije_u, llm)   : potential temperature (caldyn, IN)
      !!                REAL ps(ijb_u:ije_u)         : surface pressure (caldyn, IN)
      !!                REAL phis(ijb_u:ije_u)        : surface geopotential (caldyn, IN)
      !!                INT itau    : current dynamics time step
      !!                INT iphysiq : physics period in dynamics steps
      !!----------------------------------------------------------------------
      USE paramet_mod_h
      USE parallel_lmdz
      USE comgeom_mod_h, ONLY: cu, cv
      USE phys_state_var_mod, ONLY: phys_tstep
      USE phys_cal_mod, ONLY: days_elapsed, jH_cur, year_len
      USE comconst_mod, ONLY: daysec, pi
      USE lmdz_xios, ONLY: using_xios, xios_set_current_context
      USE mod_xios_dyn3dmem, ONLY: dyn3d_ctx_handle, writefield_dyn2d_u, writefield_dyn2d_v
      USE control_mod, ONLY: ok_dyn_xios
      ! I/O
      REAL, INTENT(IN)    :: ucov(ijb_u:ije_u, llm)
      REAL, INTENT(IN)    :: vcov(ijb_v:ije_v, llm)
      REAL, INTENT(IN)    :: teta(ijb_u:ije_u, llm)
      REAL, INTENT(IN)    :: ps(ijb_u:ije_u)
      REAL, INTENT(IN)    :: phis(ijb_u:ije_u)
      INTEGER, INTENT(IN) :: itau
      INTEGER, INTENT(IN) :: iphysiq
      ! Local variables
      INTEGER :: ktphy
      INTEGER :: j, l, ij, ijb, ije
      INTEGER :: nstep, nyear, istep
      REAL :: xcos, xsin
      !!----------------------------------------------------------------------
      !
      ! ===============
      !  1. Time step conversion
      ! ===============
      ktphy = itau / iphysiq
      !
      ! ===============
      !  2. Fill derived fields
      ! ===============
      ! bounds
      ijb = (jj_begin-1)*iip1 + 1
      ije = jj_end*iip1

      ! u = ucov / cu, pole rows zeroed (cu vanishes at the poles)
      nn_u_dyn(ijb:ije,:) = 0.
      DO j = jj_begin, jj_end
         IF (j == 1 .OR. j == jjp1) CYCLE   ! pole rows stay zero
         ij = (j-1)*iip1
         DO l = 1,llm
            nn_u_dyn(ij+1:ij+iip1,l) = ucov(ij+1:ij+iip1,l) / cu(ij+1:ij+iip1)
         END DO
      END DO
      !
      ! v = vcov / cv (no pole row on the v-grid)
      DO j = jj_begin, MIN(jj_end, jjm)
         ij = (j-1)*iip1
         DO l = 1,llm
            nn_v_dyn(ij+1:ij+iip1,l) = vcov(ij+1:ij+iip1,l) / cv(ij+1:ij+iip1)
         END DO
      END DO
      !
      ! topo = phis / 9.81 (orography in meters, as in the training data)
      nn_topo_dyn(ijb:ije) = phis(ijb:ije) / 9.81
      !
      ! Seasonal cycle scalars (same recipe as the training data, saison.py)
      nstep = NINT( daysec / phys_tstep )
      nyear = year_len * nstep
      istep = MOD( days_elapsed*nstep + INT(jh_cur*nstep), nyear ) + 1
      xcos = COS( 2.*pi*istep/nyear ) * SQRT(2.)
      xsin = SIN( 2.*pi*istep/nyear ) * SQRT(2.)
      nn_cosday = xcos
      nn_sinday = xsin
      !
      ! ===============
      !  3. OASIS exchange
      ! ===============
      ! Send inputs
      CALL send_to_python('u',       nn_u_dyn,    ktphy, pycpl_dyn_u)
      CALL send_to_python('v',       nn_v_dyn,    ktphy, pycpl_dyn_v)
      CALL send_to_python('tpot',    teta(ijb:ije,:),   ktphy, pycpl_dyn_u)
      CALL send_to_python('psol',    ps(ijb:ije),       ktphy, pycpl_dyn_u)       
      CALL send_to_python('topo',    nn_topo_dyn, ktphy, pycpl_dyn_u)
      CALL send_to_python('cos_day', nn_cosday,   ktphy)
      CALL send_to_python('sin_day', nn_sinday,   ktphy)
      ! Receive results
      CALL receive_from_python('du', py_du_dyn, ktphy, pycpl_dyn_u)
      CALL receive_from_python('dv', py_dv_dyn, ktphy, pycpl_dyn_v)
      !
      ! ===============
      !  4. Output of the received corrections
      ! ===============
      !IF (using_xios .AND. ok_dyn_xios) THEN
      !   CALL xios_set_current_context(dyn3d_ctx_handle)
      !   CALL writefield_dyn2d_u('py_du', buf_du_phys(ij_begin:ij_end,:))
      !   CALL writefield_dyn2d_v('py_dv', buf_dv_phys(ij_begin:ij_end,:))
      !ENDIF
      !
   END SUBROUTINE pycpl_dyn_fields

END MODULE pyfld
