! 
! $Id: mod_const_para.F90 1279 2009-12-10 09:02:56Z fairhead $
!
MODULE mod_const_mpi
  IMPLICIT NONE
  INTEGER,SAVE :: COMM_LMDZ
  INTEGER,SAVE :: MPI_REAL_LMDZ
 

CONTAINS 

  SUBROUTINE Init_const_mpi
    USE lmdz_mpi
    USE IOIPSL, ONLY: getin

! Use of Oasis-MCT coupler 
#ifdef CPP_OMCT
    USE mod_prism
#endif
    use wxios_mod, only: wxios_init, using_xios
    USE pycpl, ONLY: lk_pycpl
    IMPLICIT NONE

    INTEGER             :: ierr
    INTEGER             :: comp_id
    INTEGER             :: thread_required
    INTEGER             :: thread_provided
    CHARACTER(len = 6)  :: type_ocean

!$OMP MASTER
    type_ocean = 'force '
    CALL getin('type_ocean', type_ocean)
!$OMP END MASTER
!$OMP BARRIER

    IF (using_mpi) THEN
      IF (type_ocean=='couple' .OR. lk_pycpl) THEN
#ifdef CPP_COUPLE
!$OMP MASTER
        IF (using_xios) THEN
          CALL prism_init_comp_proto (comp_id, 'LMDZ', ierr)
          CALL prism_get_localcomm_proto(COMM_LMDZ,ierr)
          CALL wxios_init("LMDZ", locom=COMM_LMDZ, outcom=COMM_LMDZ, type_ocean=type_ocean)
        ELSE
          CALL prism_init_comp_proto (comp_id, 'LMDZ', ierr)
          CALL prism_get_localcomm_proto(COMM_LMDZ,ierr)
        ENDIF
!$OMP END MASTER
#endif
        MPI_REAL_LMDZ=MPI_REAL8
      ELSE
        CALL init_mpi
      ENDIF
    ENDIF
  END SUBROUTINE Init_const_mpi
  
  SUBROUTINE Init_mpi
    USE lmdz_mpi
    use wxios_mod, only: wxios_init, using_xios

  IMPLICIT NONE
    INTEGER             :: ierr
    INTEGER             :: thread_required
    INTEGER             :: thread_provided

!$OMP MASTER
      thread_required=MPI_THREAD_SERIALIZED

      CALL MPI_INIT_THREAD(thread_required,thread_provided,ierr)
      IF (thread_provided < thread_required) THEN
        PRINT *,'Warning : The multithreaded level of MPI librairy do not provide the requiered level',  &
                ' in mod_const_mpi::Init_const_mpi'
      ENDIF
      COMM_LMDZ=MPI_COMM_WORLD
      MPI_REAL_LMDZ=MPI_REAL8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialisation de XIOS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (using_xios) THEN
        WRITE(*,*)'IN Init_mpi call wxios_init'
        CALL wxios_init("LMDZ", outcom=COMM_LMDZ)
      ENDIF
!$OMP END MASTER

END SUBROUTINE Init_mpi
    
END MODULE mod_const_mpi
