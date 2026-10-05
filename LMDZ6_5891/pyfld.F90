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
   USE Bands, ONLY: distrib_physic

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
   !!$OMP THREADPRIVATE(ptemplate_phys)

   !!----------------------------------------------------------------------
   !!                   Fields on the dynamic grid
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
      ! Allocate arrays
      IF ( lk_pycpl ) THEN
 !$OMP MASTER
      ! Allocate arrays
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
         ! Dynamics
         ALLOCATE( py_du_dyn(distrib_physic%ijb_u:distrib_physic%ije_u, llm) )
         ALLOCATE( py_dv_dyn(distrib_physic%ijb_v:distrib_physic%ije_v, llm) )
         !CALL allocate_u(py_du_dyn,llm,d)
         !CALL allocate_v(py_dv_dyn,llm,d)
         !CALL allocate_u(template_u_dyn,llm)
         !CALL allocate_v(template_v_dyn,llm)
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
         !DEALLOCATE( template_phys, template_dyn)
      END IF
 !$OMP END MASTER
      !
   END SUBROUTINE pyfld_dealloc


   SUBROUTINE phys_to_dyn_u(fld_phys, fld_dyn)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE phys_to_dyn_u  ***
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
      USE lmdz_mpi
      ! I/O
      REAL, INTENT(OUT)             :: fld_dyn(iip1,distrib_physic%jjb_u:distrib_physic%jje_u, llm)
      REAL, INTENT(IN)              :: fld_phys(klon, llm)
      ! Local variables
      INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: Status
      INTEGER, DIMENSION(2) :: Req
      REAL, DIMENSION(klon_mpi,llm) :: zbuf
      REAL, DIMENSION(klon_mpi+iim,llm) :: zbuf2
      INTEGER :: i, j ,l, istart, iend, ig0, ierr
      REAL,SAVE,DIMENSION(1:iim,1:llm) :: du_send, du_recv
      !!----------------------------------------------------------------------
      !
      zbuf = 0.0
      zbuf2 = 0.0
      DO l = 1, llm
         DO i = 1, klon_omp
            zbuf(klon_omp_begin-1+i,l) = fld_phys(i,l)
         ENDDO
      ENDDO
!$OMP BARRIER
      IF (using_mpi) THEN
         IF (MPI_rank>0) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l = 1, llm
               du_send(1:iim,l) = zbuf(1:iim,l)
            ENDDO
!$OMP END DO NOWAIT
!$OMP BARRIER
!$OMP MASTER
!$OMP CRITICAL (MPI)
            CALL MPI_ISSEND(du_send,iim*llm,MPI_REAL8,MPI_Rank-1,401,COMM_LMDZ,Req(1),ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
         ENDIF

         IF (MPI_Rank<MPI_Size-1) THEN
!$OMP BARRIER
!$OMP MASTER
!$OMP CRITICAL (MPI)
            CALL MPI_IRECV(du_recv,iim*llm,MPI_REAL8,MPI_Rank+1,401,COMM_LMDZ,Req(2),ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
         ENDIF
!$OMP BARRIER

!$OMP MASTER
!$OMP CRITICAL (MPI)
         IF (MPI_Rank>0 .AND. MPI_rank < MPI_Size-1) THEN
            CALL MPI_WAITALL(2,Req(1),Status,ierr)
         ELSE IF (MPI_rank>0) THEN
            CALL MPI_WAITALL(1,Req(1),Status,ierr)
         ELSE IF (MPI_rank < MPI_Size-1) THEN
            CALL MPI_WAITALL(1,Req(2),Status,ierr)
         ENDIF
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
!$OMP BARRIER
      ENDIF

!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
         zbuf2(1:klon_mpi,l) = zbuf(1:klon_mpi,l)
         zbuf2(klon+1:klon+iim,l) = du_recv(1:iim,l)
      ENDDO
!$OMP END DO NOWAIT
      !
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
               fld_dyn(i,j,l) = 0.5*(zbuf2(ig0,l) + zbuf2(ig0+1,l))
            ENDIF

            IF (i.EQ.1) THEN
               fld_dyn(iim,j,l) = 0.5*(zbuf2(ig0,l) + zbuf2(ig0+iim-1,l))
               fld_dyn(iip1,j,l) = 0.5*(zbuf2(ig0,l) + zbuf2(ig0+1,l))
            ENDIF
         ENDDO
         !
         ! Pole rows
         IF (is_north_pole_dyn) THEN
            DO i=1,iip1
               fld_dyn(i,1,l) = 0.0
           ENDDO
         ENDIF
         !
         IF (is_south_pole_dyn) THEN
            DO i=1,iip1
               fld_dyn(i,jjp1,l) = 0.0
            ENDDO
         ENDIF
      ENDDO
!$OMP END DO NOWAIT
      !
      ! Dyn distrib
!      CALL SetTag(Request_physic,800)
!      CALL Register_SwapField_u(fld_dyn,fld_caldyn,distrib_caldyn,Request_physic)
!      CALL SendRequest(Requet_Physic)
!!$OMP BARRIER
!      CALL WaitRequest(Requet_Physic)
!!$OMP BARRIER
      !
   END SUBROUTINE phys_to_dyn_u


   SUBROUTINE phys_to_dyn_v(fld_phys, fld_dyn)
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE phys_to_dyn_v  ***
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
      USE lmdz_mpi
      ! I/O
      REAL, INTENT(OUT)             :: fld_dyn(iip1,distrib_physic%jjb_v:distrib_physic%jje_v, llm)
      REAL, INTENT(IN)              :: fld_phys(klon, llm)
      ! Local variables
      INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: Status
      INTEGER, DIMENSION(2) :: Req
      REAL, DIMENSION(klon_mpi,llm) :: zbuf
      REAL, DIMENSION(klon_mpi+iim,llm) :: zbuf2
      INTEGER :: i, j ,l, istart, iend, ig0, ierr
      REAL,SAVE,DIMENSION(1:iim,1:llm) :: dv_send, dv_recv
      !!----------------------------------------------------------------------
      !
      zbuf = 0.0
      zbuf2 = 0.0
      DO l = 1, llm
         DO i = 1, klon_omp
            zbuf(klon_omp_begin-1+i,l) = fld_phys(i,l)
         ENDDO
      ENDDO
      !
      IF (using_mpi) THEN
         IF (MPI_rank>0) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
            DO l = 1, llm
               dv_send(1:iim,l) = zbuf(1:iim,l)
            ENDDO
!$OMP END DO NOWAIT
!$OMP BARRIER
!$OMP MASTER
!$OMP CRITICAL (MPI)
            CALL MPI_ISSEND(dv_send,iim*llm,MPI_REAL8,MPI_Rank-1,402,COMM_LMDZ,Req(1),ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
         ENDIF

         IF (MPI_Rank<MPI_Size-1) THEN
!$OMP BARRIER
!$OMP MASTER
!$OMP CRITICAL (MPI)
            CALL MPI_IRECV(dv_recv,iim*llm,MPI_REAL8,MPI_Rank+1,402,COMM_LMDZ,Req(2),ierr)
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
         ENDIF
!$OMP BARRIER

!$OMP MASTER
!$OMP CRITICAL (MPI)
         IF (MPI_Rank>0 .AND. MPI_rank < MPI_Size-1) THEN
            CALL MPI_WAITALL(2,Req(1),Status,ierr)
         ELSE IF (MPI_rank>0) THEN
            CALL MPI_WAITALL(1,Req(1),Status,ierr)
         ELSE IF (MPI_rank < MPI_Size-1) THEN
            CALL MPI_WAITALL(1,Req(2),Status,ierr)
         ENDIF
!$OMP END CRITICAL (MPI)
!$OMP END MASTER
!$OMP BARRIER
      ENDIF
      !
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1, llm
         zbuf2(1:klon_mpi,l) = zbuf(1:klon_mpi,l)
         zbuf2(klon+1:klon+iim,l) = dv_recv(1:iim,l)
      ENDDO
!$OMP END DO NOWAIT
      !
      istart = 1
      iend = klon_mpi
      !
      IF (is_north_pole_dyn) istart = 2
      IF (is_south_pole_dyn) iend = klon_mpi-1-iim
      !
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
      DO l = 1,llm
         DO ig0 = istart, iend
            i = index_i(ig0)
            j = index_j(ig0)

            fld_dyn(i,j,l) = 0.5*(zbuf2(ig0,l) + zbuf2(ig0+iim,l))

            IF (i.EQ.1) THEN
               fld_dyn(iip1,j,l) = 0.5*(zbuf2(ig0,l) + zbuf2(ig0+iim,l))
            ENDIF
         ENDDO
      ENDDO
!$OMP END DO NOWAIT
         !
         ! Poles
      IF (is_north_pole_dyn) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l = 1,llm
            ! DUPLIQUE POUR LE MOMENT CAR POUR LE VENT IL FAUT LES DEUX COMPOSANTES
            DO i = 1, iim
               fld_dyn(i,1,l) = zbuf2(i+1,l)
               fld_dyn(i,1,l) = 0.5*(fld_dyn(i,1,l) + zbuf2(i+1,l))
            ENDDO
            fld_dyn(iip1,1,l) = fld_dyn(1,1,l)
         ENDDO
!$OMP END DO NOWAIT
      ENDIF
      !
      IF (is_south_pole_dyn) THEN
!$OMP DO SCHEDULE(STATIC,OMP_CHUNK)
         DO l = 1,llm
            ! DUPLIQUE POUR LE MOMENT CAR POUR LE VENT IL FAUT LES DEUX COMPOSANTES
            DO i =1,iim
               fld_dyn(i,jjm,l) = zbuf2(klon_mpi-iip1,l)
               fld_dyn(i,jjm,l) = 0.5*(fld_dyn(i,jjm,l) + zbuf2(klon_mpi-iip1,l))
            ENDDO
            fld_dyn(iip1,jjm,l) = fld_dyn(1,jjm,l)
         ENDDO
!$OMP END DO NOWAIT
      ENDIF
      ! Dyn distrib
!      CALL SetTag(Request_physic,800)
!      CALL Register_SwapField_v(fld_dyn,fld_caldyn,distrib_caldyn,Request_physic)
!      CALL SendRequest(Requet_Physic)
!!$OMP BARRIER
!      CALL WaitRequest(Requet_Physic)
!!$OMP BARRIER
      !
    END SUBROUTINE phys_to_dyn_v

END MODULE pyfld
