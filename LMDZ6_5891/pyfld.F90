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
   !!                    2D Python coupling Module fields
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: fld_a, fld_b, fld_res_a, fld_res_b  !: dummy field to store 2D fields

   !!----------------------------------------------------------------------
   !!                    3D Python coupling Module fields
   !!----------------------------------------------------------------------
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: fld_3D  !: dummy field to store 3D fields

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
         ALLOCATE( fld_a(nbp_lon,jj_nb), fld_b(nbp_lon,jj_nb), fld_res_a(nbp_lon,jj_nb), fld_res_b(nbp_lon,jj_nb) )
         fld_a = 5.0
         fld_b = -5.0
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
         DEALLOCATE( fld_a, fld_b, fld_res_a, fld_res_b )
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_dealloc

END MODULE pyfld

