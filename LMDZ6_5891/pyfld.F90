MODULE pyfld
   !!======================================================================
   !!                       ***  MODULE  pyfld  ***
   !! Python module : fields exchanged with Python scripts stored in core memory
   !!======================================================================
   !! History :  LMDZ6  ! 2026-06  (A. Barge)  Original code
   !!----------------------------------------------------------------------
   USE pycpl
   USE dimphy, ONLY: klon
   USE dimensions_mod, ONLY: iim, llm, jjm
   USE mod_phys_lmdz_mpi_data, ONLY: jj_begin, jj_end, klon_mpi
   USE mod_phys_lmdz_omp_data, ONLY: klon_omp
   USE mod_grid_phy_lmdz, ONLY: nbp_lon
   USE mod_interface_dyn_phys, ONLY: index_i, index_j
   USE Bands, ONLY: distrib_physic, distrib_caldyn

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!                   Fields on the coupling grid
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: nn_cosday, nn_sinday
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)  :: template_cpl

   !!----------------------------------------------------------------------
   !!                   Fields on the physics grid
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: py_du_phys, py_dv_phys
   !$OMP THREADPRIVATE(py_du_phys,py_dv_phys)
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: template_phys
   !!$OMP THREADPRIVATE(template_phys)

   !!----------------------------------------------------------------------
   !!                   Fields on the dynamic grid
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: py_du_dyn, py_dv_dyn
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   :: py_du_caldyn, py_dv_caldyn

   !!-----------------------------------------------------------------------
   !!        OMP working buffers - phys --> dynamic grid transformation
   !!-----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zbuf_u, zbuf_v   
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zbuf_u2, zbuf_v2

CONTAINS

   SUBROUTINE pyfld_alloc()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pyfld_alloc  ***
      !!
      !! ** Purpose :   Initialisation of the Python coupling working arrays
      !!
      !! ** Method  :   * Allocate arrays for Python fields
      !!----------------------------------------------------------------------
      USE paramet_mod_h, ONLY: iip1
      !
      ! Allocate arrays
      IF ( lk_pycpl ) THEN
 !$OMP MASTER
         ! Coupling grid
         ALLOCATE( nn_cosday(nbp_lon,jj_nb,nbp_lev), nn_sinday(nbp_lon,jj_nb,nbp_lev) )
         !ALLOCATE( template_cpl(nbp_lon,jj_nb,nbp_lev) )
 !$OMP END MASTER
         !
         ! Physics grid
         ALLOCATE( py_du_phys(klon, klev) )
         ALLOCATE( py_dv_phys(klon, klev) )
         !ALLOCATE( template_phys(klon, klev) )
         !
 !$OMP MASTER
         ! Dynamics : physics band
         ALLOCATE( py_du_dyn(iip1, distrib_physic%jjb_u:distrib_physic%jje_u, llm) )
         ALLOCATE( py_dv_dyn(iip1, distrib_physic%jjb_v:distrib_physic%jje_v, llm) )
         ! Dynamics : caldyn band
         ALLOCATE( py_du_caldyn(distrib_caldyn%ijb_u:distrib_caldyn%ije_u, llm) )
         ALLOCATE( py_dv_caldyn(distrib_caldyn%ijb_v:distrib_caldyn%ije_v, llm) )
         py_du_dyn = 0.0
         py_dv_dyn = 0.0
         py_du_caldyn = 0.0
         py_dv_caldyn = 0.0
         ! Shared OMP working buffers
         ALLOCATE( zbuf_u(klon_mpi, llm), zbuf_v(klon_mpi, llm) )
         ALLOCATE( zbuf_u2(klon_mpi+iim, llm), zbuf_v2(klon_mpi+iim, llm) )
 !$OMP END MASTER
      END IF
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
         DEALLOCATE( py_du_dyn, py_dv_dyn )
         DEALLOCATE( py_du_caldyn, py_dv_caldyn )
         DEALLOCATE( zbuf_u, zbuf_v, zbuf_u2, zbuf_v2 )
         !DEALLOCATE( template_phys, template_cpl)
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_dealloc


   SUBROUTINE pyfld_swap_ph2cd(fld_du_dyn, fld_dv_dyn, fld_du_caldyn, fld_dv_caldyn)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pyfld_swap_ph2cd  ***
      !!
      !! ** Purpose :   Swap U and V wind fields from the physics distribution
      !!                (distrib_physic) to the dynamics distribution (distrib_caldyn)
      !!
      !! ** Method  :   Mimics the exit swap block of call_calfis()
      !!----------------------------------------------------------------------
      USE parallel_lmdz
      USE mod_hallo
      ! I/O
      REAL, INTENT(IN)  :: fld_du_dyn(distrib_physic%ijb_u:distrib_physic%ije_u, llm)
      REAL, INTENT(IN)  :: fld_dv_dyn(distrib_physic%ijb_v:distrib_physic%ije_v, llm)
      REAL, INTENT(OUT) :: fld_du_caldyn(distrib_caldyn%ijb_u:distrib_caldyn%ije_u, llm)
      REAL, INTENT(OUT) :: fld_dv_caldyn(distrib_caldyn%ijb_v:distrib_caldyn%ije_v, llm)
      ! Local variables
      TYPE(Request), SAVE :: Req_swap
      !$OMP THREADPRIVATE(Req_swap)
      !!----------------------------------------------------------------------
      !
      CALL SetTag(Req_swap, 800)
      CALL Register_SwapField_u(fld_du_dyn, fld_du_caldyn, distrib_caldyn, Req_swap)
      CALL Register_SwapField_v(fld_dv_dyn, fld_dv_caldyn, distrib_caldyn, Req_swap)
      CALL SendRequest(Req_swap)
      !$OMP BARRIER
      CALL WaitRequest(Req_swap)
      !$OMP BARRIER
      !
   END SUBROUTINE pyfld_swap_ph2cd


   SUBROUTINE phys_to_dyn(fld_phys_u, fld_phys_v, fld_dyn_u, fld_dyn_v)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE phys_to_dyn  ***
      !!
      !! ** Purpose :   Re-index wind fields from the compact physics grid
      !!                onto the distrib_physic dynamics grid. Mimics calfis_loc()
      !!
      !! ** Arguments : REAL fld_phys_u(klon, llm)         : zonal wind (IN)
      !!                REAL fld_phys_v(klon, llm)         : meridional wind (IN)
      !!                REAL fld_dyn_u(iip1, jjb_u:jje_u, llm) : zonal wind (OUT)
      !!                REAL fld_dyn_v(iip1, jjb_v:jje_v, llm) : meridional wind (OUT)
      !!----------------------------------------------------------------------
      USE paramet_mod_h
      USE parallel_lmdz
      USE lmdz_mpi
      USE comgeom2_mod_h
      USE mod_hallo
      ! I/O
      REAL, INTENT(IN)              :: fld_phys_u(klon, llm)
      REAL, INTENT(IN)              :: fld_phys_v(klon, llm)
      REAL, INTENT(OUT)             :: fld_dyn_u(iip1,distrib_physic%jjb_u:distrib_physic%jje_u, llm)
      REAL, INTENT(OUT)             :: fld_dyn_v(iip1,distrib_physic%jjb_v:distrib_physic%jje_v, llm)
      ! Local variables
      INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: Status
      INTEGER, DIMENSION(4) :: Req
      INTEGER :: i, j ,l, istart, iend, ig0, ierr
      REAL,SAVE,DIMENSION(1:iim,1:llm) :: du_send, du_recv, dv_send, dv_recv
      !!----------------------------------------------------------------------
      !
!$OMP MASTER
      zbuf_u = 0.0
      zbuf_v = 0.0
      zbuf_u2 = 0.0
      zbuf_v2 = 0.0
!$OMP END MASTER
!$OMP BARRIER
      !
      DO l = 1, llm
         DO i = 1, klon_omp
            zbuf_u(klon_omp_begin-1+i,l) = fld_phys_u(i,l)
            zbuf_v(klon_omp_begin-1+i,l) = fld_phys_v(i,l)
         ENDDO
      ENDDO
!$OMP BARRIER
      IF (using_mpi) THEN
         IF (MPI_rank>0) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l = 1, llm
               du_send(1:iim,l) = zbuf_u(1:iim,l)
               dv_send(1:iim,l) = zbuf_v(1:iim,l)
            ENDDO
!$OMP END DO NOWAIT
!$OMP BARRIER
!$OMP MASTER
!$OMP CRITICAL (MPI)
            CALL MPI_ISSEND(du_send,iim*llm,MPI_REAL8,MPI_Rank-1,401,COMM_LMDZ,Req(1),ierr)
            CALL MPI_ISSEND(dv_send,iim*llm,MPI_REAL8,MPI_Rank-1,402,COMM_LMDZ,Req(2),ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
         ENDIF

         IF (MPI_Rank<MPI_Size-1) THEN
!$OMP BARRIER
!$OMP MASTER
!$OMP CRITICAL (MPI)
            CALL MPI_IRECV(du_recv,iim*llm,MPI_REAL8,MPI_Rank+1,401,COMM_LMDZ,Req(3),ierr)
            CALL MPI_IRECV(dv_recv,iim*llm,MPI_REAL8,MPI_Rank+1,402,COMM_LMDZ,Req(4),ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
         ENDIF
!$OMP BARRIER

!$OMP MASTER
!$OMP CRITICAL (MPI)
         IF (MPI_Rank>0 .AND. MPI_rank < MPI_Size-1) THEN
            CALL MPI_WAITALL(4,Req(1),Status,ierr)
         ELSE IF (MPI_rank>0) THEN
            CALL MPI_WAITALL(2,Req(1),Status,ierr)
         ELSE IF (MPI_rank < MPI_Size-1) THEN
            CALL MPI_WAITALL(2,Req(3),Status,ierr)
         ENDIF
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
!$OMP BARRIER
      ENDIF

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
         zbuf_u2(1:klon_mpi,l) = zbuf_u(1:klon_mpi,l)
         zbuf_u2(klon_mpi+1:klon_mpi+iim,l) = du_recv(1:iim,l)
         zbuf_v2(1:klon_mpi,l) = zbuf_v(1:klon_mpi,l)
         zbuf_v2(klon_mpi+1:klon_mpi+iim,l) = dv_recv(1:iim,l)
      ENDDO
!$OMP END DO NOWAIT
      !
      ! ==========
      !   U grid
      ! ==========
      istart = 1
      iend = klon_mpi
      !
      IF (is_north_pole_dyn) istart = 2
      IF (is_south_pole_dyn) iend = klon_mpi-1
      !
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
         DO ig0 = istart, iend
            i = index_i(ig0)
            j = index_j(ig0)

            IF (i.NE.iim) THEN
               fld_dyn_u(i,j,l) = 0.5*(zbuf_u2(ig0,l) + zbuf_u2(ig0+1,l))
            ENDIF

            IF (i.EQ.1) THEN
               fld_dyn_u(iim,j,l) = 0.5*(zbuf_u2(ig0,l) + zbuf_u2(ig0+iim-1,l))
               fld_dyn_u(iip1,j,l) = 0.5*(zbuf_u2(ig0,l) + zbuf_u2(ig0+1,l))
            ENDIF
         ENDDO
         !
         ! Pole rows
         IF (is_north_pole_dyn) THEN
            DO i=1,iip1
               fld_dyn_u(i,1,l) = 0.0
           ENDDO
         ENDIF
         !
         IF (is_south_pole_dyn) THEN
            DO i=1,iip1
               fld_dyn_u(i,jjp1,l) = 0.0
            ENDDO
         ENDIF
      ENDDO
!$OMP END DO NOWAIT
      !
      ! ==========
      !   V grid
      ! ==========
      IF (is_north_pole_dyn) istart = 2
      IF (is_south_pole_dyn) iend = klon_mpi-1-iim
      !
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1,llm
         DO ig0 = istart, iend
            i = index_i(ig0)
            j = index_j(ig0)

            fld_dyn_v(i,j,l) = 0.5*(zbuf_v2(ig0,l) + zbuf_v2(ig0+iim,l))

            IF (i.EQ.1) THEN
               fld_dyn_v(iip1,j,l) = 0.5*(zbuf_v2(ig0,l) + zbuf_v2(ig0+iim,l))
            ENDIF
         ENDDO
      ENDDO
!$OMP END DO NOWAIT
      !
      ! Poles
      IF (is_north_pole_dyn) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l = 1,llm
            DO i = 1, iim
               fld_dyn_v(i,1,l) = zbuf_u2(1,l)*COS(rlonv(i)) +  zbuf_v2(1,l)*SIN(rlonv(i))
               fld_dyn_v(i,1,l) = 0.5*(fld_dyn_v(i,1,l) + zbuf_v2(i+1,l))
            ENDDO
            fld_dyn_v(iip1,1,l) = fld_dyn_v(1,1,l)
         ENDDO
!$OMP END DO NOWAIT
      ENDIF
      !
      IF (is_south_pole_dyn) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l = 1,llm
            DO i =1,iim
               fld_dyn_v(i,jjm,l) = zbuf_u2(klon_mpi,l)*COS(rlonv(i)) + zbuf_v2(klon_mpi,l)*SIN(rlonv(i))
               fld_dyn_v(i,jjm,l) = 0.5*(fld_dyn_v(i,jjm,l) + zbuf_v2(klon_mpi-iip1,l))
            ENDDO
            fld_dyn_v(iip1,jjm,l) = fld_dyn_v(1,jjm,l)
         ENDDO
!$OMP END DO NOWAIT
      ENDIF
!$OMP BARRIER
      !
   END SUBROUTINE phys_to_dyn

END MODULE pyfld
