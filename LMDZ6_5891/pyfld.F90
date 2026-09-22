MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE  pyfld  ***
   !! Python module : fields exchanged with Python scripts stored in core memory
   !! Also contains the dynamics <-> python coupling exchange routine
   !!======================================================================
   !! History :  2026-06  (A. Barge)  Original code
   !!            2026-09  (A. Barge)  Swap-based exchange via distrib_physic
   !!----------------------------------------------------------------------
   USE pycpl
   USE Bands, ONLY: distrib_caldyn, distrib_physic
   USE dimensions_mod, ONLY: llm
   USE mod_phys_lmdz_mpi_data, ONLY: phy_jj_begin => jj_begin, phy_jj_end => jj_end
   USE comgeom_mod_h, ONLY: cu, cv
   USE phys_state_var_mod, ONLY: phys_tstep
   USE phys_cal_mod, ONLY: days_elapsed, jH_cur, year_len
   USE comconst_mod, ONLY: daysec, pi

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!          Seasonal cycle scalars on the coupling grid
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: nn_cosday, nn_sinday

   !!----------------------------------------------------------------------
   !!          Coupling exchange buffers on the physics distribution
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: buf_u_phys, buf_v_phys
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: buf_teta_phys
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    :: buf_ps_phys, buf_topo_phys
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: buf_du_phys, buf_dv_phys

   !!----------------------------------------------------------------------
   !!          NN wind corrections on the dynamics grid (distrib_caldyn)
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: py_du_dyn, py_dv_dyn

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
      IF ( lk_pycpl ) THEN
         ! Coupling grid
         ALLOCATE( nn_cosday(nbp_lon,jj_nb,nbp_lev), nn_sinday(nbp_lon,jj_nb,nbp_lev) )
         ! Dyanmics arrays on physics distribution buffers
         ALLOCATE( buf_u_phys(distrib_physic%ijb_u:distrib_physic%ije_u, llm) )
         ALLOCATE( buf_v_phys(distrib_physic%ijb_v:distrib_physic%ije_v, llm) )
         ALLOCATE( buf_teta_phys(distrib_physic%ijb_u:distrib_physic%ije_u, llm) )
         ALLOCATE( buf_ps_phys(distrib_physic%ijb_u:distrib_physic%ije_u) )
         ALLOCATE( buf_topo_phys(distrib_physic%ijb_u:distrib_physic%ije_u) )
         ALLOCATE( buf_du_phys(distrib_physic%ijb_u:distrib_physic%ije_u, llm) )
         ALLOCATE( buf_dv_phys(distrib_physic%ijb_v:distrib_physic%ije_v, llm) )
         ! Dynamics arrays on dyn distrib
         ALLOCATE( py_du_dyn(distrib_caldyn%ijb_u:distrib_caldyn%ije_u, llm) )
         ALLOCATE( py_dv_dyn(distrib_caldyn%ijb_v:distrib_caldyn%ije_v, llm) )
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
      IF ( lk_pycpl ) THEN
         DEALLOCATE( nn_cosday, nn_sinday )
         IF ( ALLOCATED(buf_u_phys) ) DEALLOCATE( buf_u_phys )
         IF ( ALLOCATED(buf_v_phys) ) DEALLOCATE( buf_v_phys )
         IF ( ALLOCATED(buf_teta_phys) ) DEALLOCATE( buf_teta_phys )
         IF ( ALLOCATED(buf_ps_phys) ) DEALLOCATE( buf_ps_phys )
         IF ( ALLOCATED(buf_topo_phys) ) DEALLOCATE( buf_topo_phys )
         IF ( ALLOCATED(buf_du_phys) ) DEALLOCATE( buf_du_phys )
         IF ( ALLOCATED(buf_dv_phys) ) DEALLOCATE( buf_dv_phys )
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
      !!                swap dynamics fields to the physics distribution,
      !!                send / receive via OASIS, swap results back.
      !!
      !! ** Arguments : REAL ucov(ijb_u:ije_u, llm)   : zonal covariant wind (caldyn, IN)
      !!                REAL vcov(ijb_v:ije_v, llm)   : meridional covariant wind (caldyn, IN)
      !!                REAL teta(ijb_u:ije_u, llm)   : potential temperature (caldyn, IN)
      !!                REAL ps(ijb_u:ije_u)         : surface pressure (caldyn, IN)
      !!                REAL phis(ijb_u:ije_u)        : surface geopotential (caldyn, IN)
      !!                INT itau    : current dynamics time step
      !!                INT iphysiq : physics period in dynamics steps
      !!
      !! ** OMP: must be called by the whole thread team (swaps are collective).
      !!----------------------------------------------------------------------
      USE parallel_lmdz
      USE mod_hallo
      USE dimensions_mod, ONLY: iim, jjm, ndm
      USE paramet_mod_h
      USE lmdz_xios, ONLY: using_xios, xios_set_current_context
      USE mod_xios_dyn3dmem, ONLY: dyn3d_ctx_handle, writefield_dyn2d_u, writefield_dyn2d_v
      USE control_mod, ONLY: ok_dyn_xios
      !
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
      INTEGER :: j, l, ij, lo_u, hi_u, lo_v, hi_v
      INTEGER :: nstep, nyear, istep
      REAL :: xcos, xsin
      TYPE(Request), SAVE :: Request_pycpl
      !$OMP THREADPRIVATE(Request_pycpl)
      !!----------------------------------------------------------------------
      !
      ! ===============
      !  1. Time step conversion and skip
      ! ===============
      ! At the first apphys (itau = iphysiq-1): ktphy = 0 -> skip
      ! (OASIS is being initialized inside physiq during this very step).
      ktphy = itau / iphysiq
      IF (ktphy < 1) RETURN
      !
      ! ===============
      !  2. Swap caldyn -> distrib_physic
      ! ===============
      CALL SetTag(Request_pycpl, 800)
      CALL Register_SwapField_u(ucov,  buf_u_phys,    distrib_physic, Request_pycpl, up=2, down=2)
      CALL Register_SwapField_v(vcov,  buf_v_phys,    distrib_physic, Request_pycpl, up=2, down=2)
      CALL Register_SwapField_u(teta,  buf_teta_phys, distrib_physic, Request_pycpl, up=2, down=2)
      CALL Register_SwapField_u(ps,    buf_ps_phys,   distrib_physic, Request_pycpl, up=2, down=2)
      CALL Register_SwapField_u(phis,  buf_topo_phys, distrib_physic, Request_pycpl, up=2, down=2)
      !
      CALL SendRequest(Request_pycpl)
      !$OMP BARRIER
      CALL WaitRequest(Request_pycpl)
      !$OMP BARRIER
      !
      ! ===============
      !  3. Fill derived fields on the physics band
      ! ===============
      ! u = ucov / cu, pole rows zeroed
      DO j = phy_jj_begin, phy_jj_end
         IF (j == 1 .OR. j == jjp1) CYCLE
         ij = (j-1)*iip1
         DO l = 1, llm
            buf_u_phys(ij+1:ij+iip1, l) = buf_u_phys(ij+1:ij+iip1, l) / cu(ij+1:ij+iip1)
         END DO
      END DO
      ! Zero the pole rows if they are in the physics band
      IF (phy_jj_begin == 1) THEN
         DO l = 1, llm
            buf_u_phys(1:iip1, l) = 0.
         END DO
      ENDIF
      IF (phy_jj_end == jjp1) THEN
         DO l = 1, llm
            buf_u_phys((jjp1-1)*iip1+1:jjp1*iip1, l) = 0.
         END DO
      ENDIF
      !
      ! v = vcov / cv, no pole row on the v-grid
      DO j = phy_jj_begin, MIN(phy_jj_end, jjm)
         ij = (j-1)*iip1
         DO l = 1, llm
            buf_v_phys(ij+1:ij+iip1, l) = buf_v_phys(ij+1:ij+iip1, l) / cv(ij+1:ij+iip1)
         END DO
      END DO
      !
      ! orography in meters
      buf_topo_phys = buf_topo_phys / 9.81
      !
      ! Seasonal cycle scalars
      nstep = NINT( daysec / phys_tstep )
      nyear = year_len * nstep
      istep = MOD( days_elapsed*nstep + INT(jH_cur*nstep), nyear ) + 1
      xcos = COS( 2.*pi*istep/nyear ) * SQRT(2.)
      xsin = SIN( 2.*pi*istep/nyear ) * SQRT(2.)
      !
      ! ===============
      !  4. OASIS exchange (MASTER only)
      ! ===============
      ! Slice the distrib_physic buffers on the local coupling band
      lo_u = (phy_jj_begin-1)*iip1 + 1
      hi_u =  phy_jj_end  *iip1
      lo_v = lo_u
      hi_v = MIN(phy_jj_end, jjm)*iip1
      !
      !$OMP MASTER
      nn_cosday = xcos
      nn_sinday = xsin
      !
      ! Send inputs (from distrib_physic buffers, flattened dynamics grids)
      CALL send_to_python('u',       buf_u_phys(lo_u:hi_u,:),     ktphy, pycpl_phys_dyn_u)
      CALL send_to_python('v',       buf_v_phys(lo_v:hi_v,:),     ktphy, pycpl_phys_dyn_v)
      CALL send_to_python('tpot',    buf_teta_phys(lo_u:hi_u,:),  ktphy, pycpl_phys_dyn_u)
      CALL send_to_python('psol',    buf_ps_phys(lo_u:hi_u),      ktphy, pycpl_phys_dyn_u)
      CALL send_to_python('topo',    buf_topo_phys(lo_u:hi_u),   ktphy, pycpl_phys_dyn_u)
      CALL send_to_python('cos_day', nn_cosday,     ktphy)
      CALL send_to_python('sin_day', nn_sinday,     ktphy)
      !
      ! Receive results
      CALL receive_from_python('du', buf_du_phys(lo_u:hi_u,:), ktphy, pycpl_phys_dyn_u)
      CALL receive_from_python('dv', buf_dv_phys(lo_v:hi_v,:), ktphy, pycpl_phys_dyn_v)
      !$OMP END MASTER
      !$OMP BARRIER
      !
      ! ===============
      !  5. Swap back
      ! ===============
      !$OMP MASTER
      CALL Set_Distrib(distrib_physic)
      !$OMP END MASTER
      !$OMP BARRIER
      !
      ! Output of the received corrections
      IF (using_xios .AND. ok_dyn_xios) THEN
         !$OMP MASTER
         CALL xios_set_current_context(dyn3d_ctx_handle)
         !$OMP END MASTER
         !$OMP BARRIER
         CALL writefield_dyn2d_u('py_du', buf_du_phys(ij_begin:ij_end,:))
         CALL writefield_dyn2d_v('py_dv', buf_dv_phys(ij_begin:ij_end,:))
      ENDIF
      !$OMP BARRIER
      !
      CALL SetTag(Request_pycpl, 800)
      CALL Register_SwapField_u(buf_du_phys, py_du_dyn, distrib_caldyn, Request_pycpl)
      CALL Register_SwapField_v(buf_dv_phys, py_dv_dyn, distrib_caldyn, Request_pycpl)
      CALL SendRequest(Request_pycpl)
      !$OMP BARRIER
      CALL WaitRequest(Request_pycpl)
      !$OMP BARRIER
      ! Restore the caldyn distribution
      !$OMP MASTER
      CALL Set_Distrib(distrib_caldyn)
      !$OMP END MASTER
      !$OMP BARRIER
      !
   END SUBROUTINE pycpl_dyn_fields
   !
END MODULE pyfld
