MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE  pyfld  ***
   !! Python module : fields exchanged with Python scripts stored in core memory
   !!======================================================================
   !! History :  LMDZ6  ! 2026-06  (A. Barge)  Original code
   !!----------------------------------------------------------------------
   USE pycpl
   USE dimensions_mod, ONLY: llm
   USE mod_phys_lmdz_mpi_data, ONLY: jj_begin, jj_end
   USE mod_grid_phy_lmdz, ONLY: nbp_lon

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!          Seasonal cycle scalars on the coupling grid
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: nn_cosday, nn_sinday

   !!----------------------------------------------------------------------
   !!          NN inputs on the dynamics grids
   !!          flattened (ij, llm) / (ij), same bounds as ucov / ps
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: nn_u_dyn, nn_v_dyn
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:)    :: nn_topo_dyn

   !!----------------------------------------------------------------------
   !!          NN wind corrections on the dynamics grids
   !!          flattened (ij, llm), same bounds as ucov / vcov
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

END MODULE pyfld
