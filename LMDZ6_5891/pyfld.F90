MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE  pyfld  ***
   !! Python module : fields exchanged with Python scripts stored in core memory
   !!======================================================================
   !! History :  LMDZ6  ! 2026-06  (A. Barge)  Original code
   !!----------------------------------------------------------------------
   USE pycpl
   USE Bands, ONLY: distrib_caldyn
   USE dimensions_mod, ONLY: llm

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
         ! Dynamics grids (same bounds as ucov / vcov / ps, halo included)
         ALLOCATE( nn_u_dyn(distrib_caldyn%ijb_u:distrib_caldyn%ije_u, llm) )
         ALLOCATE( nn_v_dyn(distrib_caldyn%ijb_v:distrib_caldyn%ije_v, llm) )
         ALLOCATE( nn_topo_dyn(distrib_caldyn%ijb_u:distrib_caldyn%ije_u) )
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
