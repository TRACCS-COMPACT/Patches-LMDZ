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
   USE mod_phys_lmdz_mpi_data, ONLY: jj_begin, jj_end, klon_mpi, ii_begin
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
   !REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)  :: template_u_dyn, template_v_dyn

   !!----------------------------------------------------------------------
   !!        Shared OMP working buffers
   !!----------------------------------------------------------------------
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zbuf_u, zbuf_v   
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zbuf_u2, zbuf_v2
   REAL, PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) :: zbuf 

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
         ALLOCATE( zbuf_u2(1-iim:klon_mpi+iim, llm), zbuf_v2(1-iim:klon_mpi+iim, llm) )
         ALLOCATE( zbuf(klon_mpi, llm) )
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
         DEALLOCATE( zbuf_u, zbuf_v, zbuf_u2, zbuf_v2, zbuf )
         !DEALLOCATE( template_phys, template_dyn)
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_dealloc


   SUBROUTINE pyfld_swap_ph2cd(fld_du_dyn, fld_dv_dyn, fld_du_caldyn, fld_dv_caldyn)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pyfld_swap_ph2cd  ***
      !!
      !! ** Purpose :   Swap U and V wind fields from distrib_physic to distrib_caldyn
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
      !!----------------------------------------------------------------------
      CALL SetTag(Req_swap, 800)
      CALL Register_SwapField_u(fld_du_dyn, fld_du_caldyn, distrib_caldyn, Req_swap)
      CALL Register_SwapField_v(fld_dv_dyn, fld_dv_caldyn, distrib_caldyn, Req_swap)
      CALL SendRequest(Req_swap)
      !$OMP BARRIER
      CALL WaitRequest(Req_swap)
      !$OMP BARRIER
      !
   END SUBROUTINE pyfld_swap_ph2cd


   SUBROUTINE pyfld_ghosts()
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE pyfld_ghosts  ***
      !!
      !! ** Purpose :   Fill the ghost columns of the shared buffers
      !!
      !! ** Method  :   One posting sequence with a request counter and a
      !!                single MPI_WAITALL : no per-rank case logic.
      !!
      !! ** OMP : collective, must be called by all the threads
      !!----------------------------------------------------------------------
      USE parallel_lmdz
      USE lmdz_mpi
      ! Local variables
      INTEGER, DIMENSION(8) :: Req_all
      INTEGER, DIMENSION(MPI_STATUS_SIZE,8) :: Status
      INTEGER :: l, nb, ierr
      REAL,SAVE,DIMENSION(1:iim,1:llm) :: du_send_south, du_recv_south, dv_send_south, dv_recv_south
      REAL,SAVE,DIMENSION(1:iim,1:llm) :: du_send_north, du_recv_north, dv_send_north, dv_recv_north
      !!----------------------------------------------------------------------
      !
      ! Extract columns
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
         du_send_south(1:iim,l) = zbuf_u(1:iim,l)
         dv_send_south(1:iim,l) = zbuf_v(1:iim,l)
         du_send_north(1:iim,l) = zbuf_u(klon_mpi-iim+1:klon_mpi,l)
         dv_send_north(1:iim,l) = zbuf_v(klon_mpi-iim+1:klon_mpi,l)
      ENDDO
!$OMP END DO NOWAIT
!$OMP BARRIER
      !
      ! Communicate columns between processes 
      nb = 0
!$OMP MASTER
!$OMP CRITICAL (MPI)
      IF (MPI_Rank>0) THEN
         nb = nb+2
         CALL MPI_ISSEND(du_send_south,iim*llm,MPI_REAL8,MPI_Rank-1,401,COMM_LMDZ,Req_all(nb-1),ierr)
         CALL MPI_ISSEND(dv_send_south,iim*llm,MPI_REAL8,MPI_Rank-1,402,COMM_LMDZ,Req_all(nb),ierr)
      ENDIF
      IF (MPI_Rank<MPI_Size-1) THEN
         nb = nb+2
         CALL MPI_IRECV(du_recv_south,iim*llm,MPI_REAL8,MPI_Rank+1,401,COMM_LMDZ,Req_all(nb-1),ierr)
         CALL MPI_IRECV(dv_recv_south,iim*llm,MPI_REAL8,MPI_Rank+1,402,COMM_LMDZ,Req_all(nb),ierr)
      ENDIF
      IF (MPI_Rank<MPI_Size-1) THEN
         nb = nb+2
         CALL MPI_ISSEND(du_send_north,iim*llm,MPI_REAL8,MPI_Rank+1,403,COMM_LMDZ,Req_all(nb-1),ierr)
         CALL MPI_ISSEND(dv_send_north,iim*llm,MPI_REAL8,MPI_Rank+1,404,COMM_LMDZ,Req_all(nb),ierr)
      ENDIF
      IF (MPI_Rank>0) THEN
         nb = nb+2
         CALL MPI_IRECV(du_recv_north,iim*llm,MPI_REAL8,MPI_Rank-1,403,COMM_LMDZ,Req_all(nb-1),ierr)
         CALL MPI_IRECV(dv_recv_north,iim*llm,MPI_REAL8,MPI_Rank-1,404,COMM_LMDZ,Req_all(nb),ierr)
      ENDIF
      CALL MPI_WAITALL(nb,Req_all,Status,ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
!$OMP BARRIER
      !
      ! Add received columns
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
         zbuf_u2(1:klon_mpi,l) = zbuf_u(1:klon_mpi,l)
         zbuf_u2(klon_mpi+1:klon_mpi+iim,l) = du_recv_south(1:iim,l)
         zbuf_u2(1-iim:0,l) = du_recv_north(1:iim,l)
         zbuf_v2(1:klon_mpi,l) = zbuf_v(1:klon_mpi,l)
         zbuf_v2(klon_mpi+1:klon_mpi+iim,l) = dv_recv_south(1:iim,l)
         zbuf_v2(1-iim:0,l) = dv_recv_north(1:iim,l)
      ENDDO
!$OMP END DO NOWAIT
      !
   END SUBROUTINE pyfld_ghosts


   SUBROUTINE phys_to_dyn(fld_phys_u, fld_phys_v, fld_dyn_u, fld_dyn_v)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE phys_to_dyn_u  ***
      !!
      !! ** Purpose :   Re-index U and V fields from the physics to dynamic grid
      !!                on distrib_physic
      !!
      !! ** Arguments : REAL fld_phys_u(klon, llm)         : U field on physics grid (IN)
      !!                REAL fld_phys_v(klon, llm)         : V field on physics grid (IN)
      !!                REAL fld_dyn_u(iip1, jjb_u:jje_u, llm) : U field on dyn grid (OUT)
      !!                REAL fld_dyn_v(iip1, jjb_v:jje_v, llm) : V field on dyn grid (OUT)
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
      INTEGER :: i, j ,l, istart, iend, ig0
      !!----------------------------------------------------------------------
      !
      ! Fill OMP shared buffers
      DO l = 1, llm
         DO i = 1, klon_omp
            zbuf_u(klon_omp_begin-1+i,l) = fld_phys_u(i,l)
            zbuf_v(klon_omp_begin-1+i,l) = fld_phys_v(i,l)
         ENDDO
      ENDDO
!$OMP BARRIER
      ! Complete columns
      CALL pyfld_ghosts()
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
         IF (ii_begin > 1) THEN
            DO i = 1, ii_begin-2
               fld_dyn_u(i,jj_begin,l) = 0.5*(zbuf_u2(i-ii_begin+1,l) + zbuf_u2(i-ii_begin+2,l))
            ENDDO
            fld_dyn_u(ii_begin-1,jj_begin,l) = 0.5*(zbuf_u2(0,l) + zbuf_u2(1,l))
            fld_dyn_u(iim ,jj_begin,l) = 0.5*(zbuf_u2(2-ii_begin,l) + zbuf_u2(iim-ii_begin+1,l))
            fld_dyn_u(iip1,jj_begin,l) = 0.5*(zbuf_u2(2-ii_begin,l) + zbuf_u2(3-ii_begin,l))
         ENDIF
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
         !
         IF (ii_begin > 1) THEN
            DO i = 1, ii_begin-1
               fld_dyn_v(i,jj_begin,l) = 0.5*(zbuf_v2(i-ii_begin+1,l) + zbuf_v2(i-ii_begin+1+iim,l))
            ENDDO
            fld_dyn_v(iip1,jj_begin,l) = fld_dyn_v(1,jj_begin,l)
         ENDIF
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
