MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE pyfld  ***
   !! Python module : fields returned by Python script stored in core memory
   !!======================================================================
   !! History :  LMDZ6  ! 2026-06  (A. Barge)  Original code
   !!----------------------------------------------------------------------
   !!
   !!----------------------------------------------------------------------
   USE pycpl
   USE dimphy, ONLY: klon, klev

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!                    3D Python coupling Module fields
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: nn_u, nn_v, nn_tpot, nn_cosday, nn_sinday
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: py_du, py_dv

   !!----------------------------------------------------------------------
   !!                    2D Python coupling Module fields
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: nn_psol, nn_topo

CONTAINS

   SUBROUTINE pyfld_alloc()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pyfld_alloc  ***
      !!
      !! ** Purpose :   Initialisation of the Python-computed fields
      !!
      !! ** Method  :   * Allocate arrays for Python fields
      !!----------------------------------------------------------------------
      !
      ! Allocate arrays
 !$OMP MASTER
      IF ( lk_pycpl ) THEN
         ALLOCATE( nn_u(nbp_lon,jj_nb,nbp_lev), nn_v(nbp_lon,jj_nb,nbp_lev), &
                 & nn_tpot(nbp_lon,jj_nb,nbp_lev), nn_cosday(nbp_lon,jj_nb,nbp_lev) )
         ALLOCATE( nn_sinday(nbp_lon,jj_nb,nbp_lev), nn_psol(nbp_lon,jj_nb), nn_topo(nbp_lon,jj_nb) )
         ALLOCATE( py_du(nbp_lon,jj_nb,nbp_lev), py_dv(nbp_lon,jj_nb,nbp_lev) )
         nn_cosday = 0.3 ! /
         nn_sinday = -0.3 ! /
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
      ! Free memory
 !$OMP MASTER
      IF ( lk_pycpl ) THEN
         DEALLOCATE( nn_u, nn_v, nn_tpot, nn_cosday, nn_sinday, nn_psol, nn_topo)
         DEALLOCATE( py_du, py_dv )
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_dealloc

END MODULE pyfld

