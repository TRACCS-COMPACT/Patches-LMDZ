! 
! $Id$
!
  MODULE parallel_lmdz
  USE mod_const_mpi
  USE lmdz_mpi, ONLY : using_mpi
      use IOIPSL

    INTEGER,PARAMETER :: halo_max=3
    
    LOGICAL,SAVE :: using_omp ! .true. if using OpenMP
    LOGICAL,SAVE :: is_master ! .true. if the core is both MPI & OpenMP master
!$OMP THREADPRIVATE(is_master)
    
    integer, save :: mpi_size
    integer, save :: mpi_rank
    integer, save :: jj_begin
    integer, save :: jj_end
    integer, save :: jj_nb
    integer, save :: ij_begin
    integer, save :: ij_end
    logical, save :: pole_nord
    logical, save :: pole_sud

    integer,save  :: jjb_u
    integer,save  :: jje_u
    integer,save  :: jjnb_u
    integer,save  :: jjb_v
    integer,save  :: jje_v
    integer,save  :: jjnb_v    

    integer,save  :: ijb_u
    integer,save  :: ije_u
    integer,save  :: ijnb_u    
    
    integer,save  :: ijb_v
    integer,save  :: ije_v
    integer,save  :: ijnb_v    
     
    
    integer, allocatable, save, dimension(:) :: jj_begin_para
    integer, allocatable, save, dimension(:) :: jj_end_para
    integer, allocatable, save, dimension(:) :: jj_nb_para
    integer, save :: OMP_CHUNK
    integer, save :: omp_rank
    integer, save :: omp_size  
!$OMP THREADPRIVATE(omp_rank)

    TYPE distrib
      integer :: jj_begin
      integer :: jj_end
      integer :: jj_nb
      integer :: ij_begin
      integer :: ij_end

      integer  :: jjb_u
      integer  :: jje_u
      integer  :: jjnb_u
      integer  :: jjb_v
      integer  :: jje_v
      integer  :: jjnb_v    
 
      integer  :: ijb_u
      integer  :: ije_u
      integer  :: ijnb_u    
    
      integer  :: ijb_v
      integer  :: ije_v
      integer  :: ijnb_v    
     
    
      integer, pointer :: jj_begin_para(:) => NULL()
      integer, pointer :: jj_end_para(:) => NULL()
      integer, pointer :: jj_nb_para(:) => NULL()
    END TYPE distrib  
    
    INTERFACE ASSIGNMENT (=)
      MODULE PROCEDURE copy_distrib
    END INTERFACE
    TYPE(distrib),SAVE :: current_dist
    
 contains
 
    subroutine init_parallel
    USE iniprint_mod_h
      USE vampir
    USE lmdz_mpi
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
implicit none



      integer :: ierr
      integer :: i,j
      integer :: type_size
      integer, dimension(3) :: blocklen,type
      integer :: comp_id
      character(len=4)  :: num
      character(len=20) :: filename
 
#ifdef CPP_OMP    
      INTEGER :: OMP_GET_NUM_THREADS
      EXTERNAL OMP_GET_NUM_THREADS
      INTEGER :: OMP_GET_THREAD_NUM
      EXTERNAL OMP_GET_THREAD_NUM
#endif  

#ifdef CPP_OMP
       using_OMP=.TRUE.
#else
       using_OMP=.FALSE.
#endif
      
      call InitVampir
      
      IF (using_mpi) THEN
        call MPI_COMM_SIZE(COMM_LMDZ,mpi_size,ierr)
        call MPI_COMM_RANK(COMM_LMDZ,mpi_rank,ierr)
      ELSE
        mpi_size=1
        mpi_rank=0
      ENDIF


! Open text output file with mpi_rank in suffix of file name 
      IF (lunout /= 5 .and. lunout /= 6) THEN
         WRITE(num,'(I4.4)') mpi_rank
         filename='lmdz.out_'//num
         IF (mpi_rank .NE. 0) THEN
            OPEN(UNIT=lunout,FILE=TRIM(filename),ACTION='write', &
               STATUS='unknown',FORM='formatted',IOSTAT=ierr) 
         ENDIF
      ENDIF

      
      allocate(jj_begin_para(0:mpi_size-1))
      allocate(jj_end_para(0:mpi_size-1))
      allocate(jj_nb_para(0:mpi_size-1))
      
      do i=0,mpi_size-1
        jj_nb_para(i)=(jjm+1)/mpi_size
        if ( i < MOD((jjm+1),mpi_size) ) jj_nb_para(i)=jj_nb_para(i)+1
        
        if (jj_nb_para(i) <= 1 ) then
          
         write(lunout,*)"Arret : le nombre de bande de lattitude par process est trop faible (<2)."
         write(lunout,*)" ---> diminuez le nombre de CPU ou augmentez la taille en lattitude"
          
          IF (using_mpi) call MPI_ABORT(COMM_LMDZ,-1, ierr)

        endif
        
      enddo
      
!      jj_nb_para(0)=11
!      jj_nb_para(1)=25
!      jj_nb_para(2)=25
!      jj_nb_para(3)=12      

      j=1
      
      do i=0,mpi_size-1 
        
        jj_begin_para(i)=j
        jj_end_para(i)=j+jj_Nb_para(i)-1
        j=j+jj_Nb_para(i)
      
      enddo
      
      jj_begin = jj_begin_para(mpi_rank)
      jj_end   = jj_end_para(mpi_rank)
      jj_nb    = jj_nb_para(mpi_rank)
      
      ij_begin=(jj_begin-1)*iip1+1
      ij_end=jj_end*iip1
      
      if (mpi_rank.eq.0) then
        pole_nord=.TRUE.
      else 
        pole_nord=.FALSE.
      endif
      
      if (mpi_rank.eq.mpi_size-1) then
        pole_sud=.TRUE.
      else 
        pole_sud=.FALSE.
      endif
        
      write(lunout,*)"init_parallel: jj_begin",jj_begin
      write(lunout,*)"init_parallel: jj_end",jj_end
      write(lunout,*)"init_parallel: ij_begin",ij_begin
      write(lunout,*)"init_parallel: ij_end",ij_end
      jjb_u=MAX(jj_begin-halo_max,1)
      jje_u=MIN(jj_end+halo_max,jjp1)
      jjnb_u=jje_u-jjb_u+1

      jjb_v=MAX(jj_begin-halo_max,1)
      jje_v=MIN(jj_end+halo_max,jjm)
      jjnb_v=jje_v-jjb_v+1

      ijb_u=MAX(ij_begin-halo_max*iip1,1)
      ije_u=MIN(ij_end+halo_max*iip1,ip1jmp1)
      ijnb_u=ije_u-ijb_u+1

      ijb_v=MAX(ij_begin-halo_max*iip1,1)
      ije_v=MIN(ij_end+halo_max*iip1,ip1jm)
      ijnb_v=ije_v-ijb_v+1
      
!$OMP PARALLEL

#ifdef CPP_OMP
!$OMP MASTER
        omp_size=OMP_GET_NUM_THREADS()
!$OMP END MASTER
!$OMP BARRIER
        omp_rank=OMP_GET_THREAD_NUM()  

!Config  Key  = omp_chunk
!Config  Desc = taille des blocs openmp
!Config  Def  = 1
!Config  Help = defini la taille des packets d'it�ration openmp
!Config         distribue a chaque tache lors de l'entree dans une
!Config         boucle parallelisee

!$OMP MASTER
      omp_chunk=(llm+1)/omp_size
      IF (MOD(llm+1,omp_size)/=0) omp_chunk=omp_chunk+1
      CALL getin('omp_chunk',omp_chunk)
!$OMP END MASTER
!$OMP BARRIER        
#else    
        omp_size=1
        omp_rank=0
#endif
!$OMP END PARALLEL         
      CALL create_distrib(jj_nb_para,current_dist)
      
      IF ((mpi_rank==0).and.(omp_rank==0)) THEN
        is_master=.true.
      ELSE
        is_master=.false.
      ENDIF
      
    end subroutine init_parallel

    SUBROUTINE create_distrib(jj_nb_new,d)
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE


      
      INTEGER,INTENT(IN) :: jj_Nb_New(0:MPI_Size-1)
      TYPE(distrib),INTENT(INOUT) :: d
      INTEGER :: i  
  
      IF (.NOT. ASSOCIATED(d%jj_nb_para)) ALLOCATE(d%jj_nb_para(0:MPI_Size-1))
      IF (.NOT. ASSOCIATED(d%jj_begin_para)) ALLOCATE(d%jj_begin_para(0:MPI_Size-1))
      IF (.NOT. ASSOCIATED(d%jj_end_para)) ALLOCATE(d%jj_end_para(0:MPI_Size-1))
      
      d%jj_Nb_Para=jj_Nb_New
      
      d%jj_begin_para(0)=1
      d%jj_end_para(0)=d%jj_Nb_Para(0)
      
      do i=1,mpi_size-1 
        
        d%jj_begin_para(i)=d%jj_end_para(i-1)+1
        d%jj_end_para(i)=d%jj_begin_para(i)+d%jj_Nb_para(i)-1
      
      enddo
      
      d%jj_begin = d%jj_begin_para(mpi_rank)
      d%jj_end   = d%jj_end_para(mpi_rank)
      d%jj_nb    = d%jj_nb_para(mpi_rank)
      
      d%ij_begin=(d%jj_begin-1)*iip1+1
      d%ij_end=d%jj_end*iip1

      d%jjb_u=MAX(d%jj_begin-halo_max,1)
      d%jje_u=MIN(d%jj_end+halo_max,jjp1)
      d%jjnb_u=d%jje_u-d%jjb_u+1

      d%jjb_v=MAX(d%jj_begin-halo_max,1)
      d%jje_v=MIN(d%jj_end+halo_max,jjm)
      d%jjnb_v=d%jje_v-d%jjb_v+1

      d%ijb_u=MAX(d%ij_begin-halo_max*iip1,1)
      d%ije_u=MIN(d%ij_end+halo_max*iip1,ip1jmp1)
      d%ijnb_u=d%ije_u-d%ijb_u+1

      d%ijb_v=MAX(d%ij_begin-halo_max*iip1,1)
      d%ije_v=MIN(d%ij_end+halo_max*iip1,ip1jm)
      d%ijnb_v=d%ije_v-d%ijb_v+1      

    END SUBROUTINE create_distrib

      
    SUBROUTINE Set_Distrib(d)
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE



    TYPE(distrib),INTENT(IN) :: d

      jj_begin = d%jj_begin
      jj_end = d%jj_end
      jj_nb = d%jj_nb
      ij_begin = d%ij_begin
      ij_end = d%ij_end

      jjb_u = d%jjb_u
      jje_u = d%jje_u
      jjnb_u = d%jjnb_u
      jjb_v = d%jjb_v
      jje_v = d%jje_v
      jjnb_v = d%jjnb_v
 
      ijb_u = d%ijb_u
      ije_u = d%ije_u
      ijnb_u = d%ijnb_u
    
      ijb_v = d%ijb_v
      ije_v = d%ije_v
      ijnb_v = d%ijnb_v
     
    
      jj_begin_para(:) = d%jj_begin_para(:)
      jj_end_para(:) = d%jj_end_para(:)
      jj_nb_para(:) = d%jj_nb_para(:)
      current_dist=d 

    END SUBROUTINE Set_Distrib

    SUBROUTINE copy_distrib(dist,new_dist)
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE



    TYPE(distrib),INTENT(INOUT) :: dist
    TYPE(distrib),INTENT(IN) :: new_dist

     dist%jj_begin = new_dist%jj_begin
     dist%jj_end = new_dist%jj_end
     dist%jj_nb = new_dist%jj_nb
     dist%ij_begin = new_dist%ij_begin
     dist%ij_end = new_dist%ij_end

     dist%jjb_u = new_dist%jjb_u
     dist%jje_u = new_dist%jje_u
     dist%jjnb_u = new_dist%jjnb_u
     dist%jjb_v = new_dist%jjb_v
     dist%jje_v = new_dist%jje_v
     dist%jjnb_v = new_dist%jjnb_v
   
     dist%ijb_u = new_dist%ijb_u
     dist%ije_u = new_dist%ije_u
     dist%ijnb_u = new_dist%ijnb_u
     
     dist%ijb_v = new_dist%ijb_v
     dist%ije_v = new_dist%ije_v
     dist%ijnb_v = new_dist%ijnb_v
         
     
     dist%jj_begin_para(:) = new_dist%jj_begin_para(:)
     dist%jj_end_para(:) = new_dist%jj_end_para(:)
     dist%jj_nb_para(:) = new_dist%jj_nb_para(:)
 
    END SUBROUTINE copy_distrib
    
    
    SUBROUTINE get_current_distrib(d)
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
USE paramet_mod_h
IMPLICIT NONE



    TYPE(distrib),INTENT(OUT) :: d

     d=current_dist

    END SUBROUTINE get_current_distrib
    
    subroutine Finalize_parallel
    USE lmdz_mpi
    ! ug Pour les sorties XIOS
  	use wxios_mod
    USE control_mod, only : ok_dyn_xios
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
    USE pycpl

#ifdef CPP_COUPLE
! Use of Oasis-MCT coupler 
#if defined CPP_OMCT
    use mod_prism
#else
    use mod_prism_proto
#endif
! Ehouarn: surface_data module is in 'phylmd' ...
      use surface_data, only : type_ocean
      implicit none
#else
      implicit none
! without the surface_data module, we declare (and set) a dummy 'type_ocean'
      character(len=6),parameter :: type_ocean="dummy"
#endif
! #endif of #ifdef CPP_EARTH
      integer :: ierr
      integer :: i

      if (allocated(jj_begin_para)) deallocate(jj_begin_para)
      if (allocated(jj_end_para))   deallocate(jj_end_para)
      if (allocated(jj_nb_para))    deallocate(jj_nb_para)

      if (type_ocean == 'couple' .OR. lk_pycpl) then
#ifdef CPP_COUPLE
        IF (lk_pycpl) THEN
           CALL finalize_python_coupling()
        ENDIF
        IF (using_xios) THEN
          !Fermeture propre de XIOS
          ! close xios dynamic context if is defined (call LMDZDYN)   
          IF (ok_dyn_xios) THEN
             CALL xios_context_finalize()
          ENDIF
          CALL wxios_close()
          CALL prism_terminate_proto(ierr)
          IF (ierr .ne. PRISM_Ok) THEN
            CALL abort_gcm('Finalize_parallel',' Probleme dans prism_terminate_proto ',1)
          ENDIF
        ELSE
           call prism_terminate_proto(ierr)
           IF (ierr .ne. PRISM_Ok) THEN
              call abort_gcm('Finalize_parallel',' Probleme dans prism_terminate_proto ',1)
           endif
        ENDIF
#else
        call abort_gcm('Finalize_parallel','type_ocean = couple but CPP_COUPLE not defined',1)
#endif 
      else
        IF (using_xios) THEN
          !Fermeture propre de XIOS
          IF (ok_dyn_xios) THEN
             CALL xios_context_finalize()
          ENDIF
          CALL wxios_close()
        ENDIF
        IF (using_mpi) call MPI_FINALIZE(ierr)
      end if
      
    end subroutine Finalize_parallel
        
    subroutine Pack_Data(Field,ij,ll,row,Buffer)
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
implicit none




      integer, intent(in) :: ij,ll,row
      real,dimension(ij,ll),intent(in) ::Field
      real,dimension(ll*iip1*row), intent(out) :: Buffer 
            
      integer :: Pos
      integer :: i,l
      
      Pos=0
      do l=1,ll
        do i=1,row*iip1
          Pos=Pos+1
          Buffer(Pos)=Field(i,l)
        enddo
      enddo
      
    end subroutine Pack_data 
     
    subroutine Unpack_Data(Field,ij,ll,row,Buffer)
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
implicit none




      integer, intent(in) :: ij,ll,row
      real,dimension(ij,ll),intent(out) ::Field
      real,dimension(ll*iip1*row), intent(in) :: Buffer 
            
      integer :: Pos
      integer :: i,l
      
      Pos=0
      
      do l=1,ll
        do i=1,row*iip1
          Pos=Pos+1
          Field(i,l)=Buffer(Pos)
        enddo
      enddo
      
    end subroutine UnPack_data

    
    SUBROUTINE barrier
    USE lmdz_mpi
    IMPLICIT NONE
    INTEGER :: ierr
    
!$OMP CRITICAL (MPI)      
      IF (using_mpi) CALL MPI_Barrier(COMM_LMDZ,ierr)
!$OMP END CRITICAL (MPI)
    
    END SUBROUTINE barrier
       
      
    subroutine exchange_hallo(Field,ij,ll,up,down)
    USE lmdz_mpi
    USE Vampir
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
implicit none


      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: up,down
      
      INTEGER :: ierr
      LOGICAL :: SendUp,SendDown
      LOGICAL :: RecvUp,RecvDown
      INTEGER, DIMENSION(4) :: Request
      INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: Status

      INTEGER :: NbRequest
      REAL, dimension(:),allocatable :: Buffer_Send_up,Buffer_Send_down
      REAL, dimension(:),allocatable :: Buffer_Recv_up,Buffer_Recv_down
      INTEGER :: Buffer_size      

      IF (using_mpi) THEN

        CALL barrier
      
        call VTb(VThallo)
      
        SendUp=.TRUE.
        SendDown=.TRUE.
        RecvUp=.TRUE.
        RecvDown=.TRUE.
          
        IF (pole_nord) THEN
          SendUp=.FALSE.
          RecvUp=.FALSE.
        ENDIF
    
        IF (pole_sud) THEN
          SendDown=.FALSE.
          RecvDown=.FALSE.
        ENDIF
        
        if (up.eq.0) then
          SendDown=.FALSE.
          RecvUp=.FALSE.
        endif
      
        if (down.eq.0) then
          SendUp=.FALSE.
          RecvDown=.FALSE.
        endif
      
        NbRequest=0
  
        IF (SendUp) THEN
          NbRequest=NbRequest+1
          buffer_size=down*iip1*ll
          allocate(Buffer_Send_up(Buffer_size))
          call PACK_Data(Field(ij_begin,1),ij,ll,down,Buffer_Send_up)
!$OMP CRITICAL (MPI)
          call MPI_ISEND(Buffer_send_up,Buffer_Size,MPI_REAL8,MPI_Rank-1,1,     &
                          COMM_LMDZ,Request(NbRequest),ierr)
!$OMP END CRITICAL (MPI)
        ENDIF
  
        IF (SendDown) THEN
          NbRequest=NbRequest+1
           
          buffer_size=up*iip1*ll
          allocate(Buffer_Send_down(Buffer_size))
          call PACK_Data(Field(ij_end+1-up*iip1,1),ij,ll,up,Buffer_send_down)
        
!$OMP CRITICAL (MPI)
          call MPI_ISEND(Buffer_send_down,Buffer_Size,MPI_REAL8,MPI_Rank+1,1,     &
                          COMM_LMDZ,Request(NbRequest),ierr)
!$OMP END CRITICAL (MPI)
        ENDIF
    
  
        IF (RecvUp) THEN
          NbRequest=NbRequest+1
          buffer_size=up*iip1*ll
          allocate(Buffer_recv_up(Buffer_size))
              
!$OMP CRITICAL (MPI)
          call MPI_IRECV(Buffer_recv_up,Buffer_size,MPI_REAL8,MPI_Rank-1,1,  &
                          COMM_LMDZ,Request(NbRequest),ierr)
!$OMP END CRITICAL (MPI)
     
       
        ENDIF
  
        IF (RecvDown) THEN
          NbRequest=NbRequest+1
          buffer_size=down*iip1*ll
          allocate(Buffer_recv_down(Buffer_size))
        
!$OMP CRITICAL (MPI)
          call MPI_IRECV(Buffer_recv_down,Buffer_size,MPI_REAL8,MPI_Rank+1,1,     &
                          COMM_LMDZ,Request(NbRequest),ierr)
!$OMP END CRITICAL (MPI)
        
        ENDIF
  
        if (NbRequest > 0) call MPI_WAITALL(NbRequest,Request,Status,ierr)
        IF (RecvUp)  call Unpack_Data(Field(ij_begin-up*iip1,1),ij,ll,up,Buffer_Recv_up)
        IF (RecvDown) call Unpack_Data(Field(ij_end+1,1),ij,ll,down,Buffer_Recv_down)  

        call VTe(VThallo)
        call barrier
      
      ENDIF  ! using_mpi
      
      RETURN
      
    end subroutine exchange_Hallo
    

    subroutine Gather_Field(Field,ij,ll,rank)
    USE iniprint_mod_h
    USE lmdz_mpi
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
implicit none


      INTEGER :: ij,ll,rank
      REAL, dimension(ij,ll) :: Field
      REAL, dimension(:),allocatable :: Buffer_send   
      REAL, dimension(:),allocatable :: Buffer_Recv
      INTEGER, dimension(0:MPI_Size-1) :: Recv_count, displ
      INTEGER :: ierr
      INTEGER ::i
      
      IF (using_mpi) THEN

        if (ij==ip1jmp1) then 
           allocate(Buffer_send(iip1*ll*(jj_end-jj_begin+1)))
           call Pack_Data(Field(ij_begin,1),ij,ll,jj_end-jj_begin+1,Buffer_send)
        else if (ij==ip1jm) then
           allocate(Buffer_send(iip1*ll*(min(jj_end,jjm)-jj_begin+1)))
           call Pack_Data(Field(ij_begin,1),ij,ll,min(jj_end,jjm)-jj_begin+1,Buffer_send)
        else
           write(lunout,*)ij  
        CALL abort_gcm("parallel_lmdz","erreur dans Gather_Field",1)
        endif
        
        if (MPI_Rank==rank) then
          allocate(Buffer_Recv(ij*ll))

!CDIR NOVECTOR
          do i=0,MPI_Size-1
             
            if (ij==ip1jmp1) then 
              Recv_count(i)=(jj_end_para(i)-jj_begin_para(i)+1)*ll*iip1
            else if (ij==ip1jm) then
              Recv_count(i)=(min(jj_end_para(i),jjm)-jj_begin_para(i)+1)*ll*iip1
            else
              CALL abort_gcm("parallel_lmdz","erreur dans Gather_Field",1)
            endif
                   
            if (i==0) then 
              displ(i)=0 
            else
              displ(i)=displ(i-1)+Recv_count(i-1)
            endif
            
          enddo
          
        else
          ! Ehouarn: When in debug mode, ifort complains (for call MPI_GATHERV
          !          below) about Buffer_Recv() being not allocated.
          !          So make a dummy allocation.
          allocate(Buffer_Recv(1))
        endif ! of if (MPI_Rank==rank)
  
!$OMP CRITICAL (MPI)
        call MPI_GATHERV(Buffer_send,(min(ij_end,ij)-ij_begin+1)*ll,MPI_REAL8,   &
                          Buffer_Recv,Recv_count,displ,MPI_REAL8,rank,COMM_LMDZ,ierr)
!$OMP END CRITICAL (MPI)
      
        if (MPI_Rank==rank) then                  
      
          if (ij==ip1jmp1) then 
            do i=0,MPI_Size-1
              call Unpack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                 &
                               jj_end_para(i)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
            enddo
          else if (ij==ip1jm) then
            do i=0,MPI_Size-1
               call Unpack_Data(Field((jj_begin_para(i)-1)*iip1+1,1),ij,ll,                       &
                               min(jj_end_para(i),jjm)-jj_begin_para(i)+1,Buffer_Recv(displ(i)+1))
            enddo
          endif
        endif 
      ENDIF ! using_mpi
      
    end subroutine Gather_Field


    subroutine AllGather_Field(Field,ij,ll)
    USE lmdz_mpi
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
implicit none


      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: ierr
      
      IF (using_mpi) THEN
        call Gather_Field(Field,ij,ll,0)
!$OMP CRITICAL (MPI)
      call MPI_BCAST(Field,ij*ll,MPI_REAL8,0,COMM_LMDZ,ierr)
!$OMP END CRITICAL (MPI)
      ENDIF
      
    end subroutine AllGather_Field
    
   subroutine Broadcast_Field(Field,ij,ll,rank)
    USE lmdz_mpi
    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
    USE paramet_mod_h
implicit none


      INTEGER :: ij,ll
      REAL, dimension(ij,ll) :: Field
      INTEGER :: rank
      INTEGER :: ierr
      
      IF (using_mpi) THEN
      
!$OMP CRITICAL (MPI)
      call MPI_BCAST(Field,ij*ll,MPI_REAL8,rank,COMM_LMDZ,ierr)
!$OMP END CRITICAL (MPI)
      
      ENDIF
    end subroutine Broadcast_Field
        
   
!  Subroutine verif_hallo(Field,ij,ll,up,down)
!    USE lmdz_mpi
!    USE dimensions_mod, ONLY: iim, jjm, llm, ndm
!    USE paramet_mod_h, ONLY: iip1, iip2, iip3, jjp1, llmp1, llmp2, llmm1, kftd, ip1jm, ip1jmp1, &
!          ip1jmi1, ijp1llm, ijmllm, mvar, jcfil, jcfllm
!implicit none
!
!
!    
!      INTEGER :: ij,ll
!      REAL, dimension(ij,ll) :: Field
!      INTEGER :: up,down 
!      
!      REAL,dimension(ij,ll): NewField
!      
!      NewField=0
!      
!      ijb=ij_begin
!      ije=ij_end
!      if (pole_nord) 
!      NewField(ij_be       

  end MODULE parallel_lmdz
