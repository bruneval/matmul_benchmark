!===================================================
! TEST symmetric diago
!
!===================================================
program diago
#if defined(SCALAPACK)
 use mpi
#endif
 use,intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
 use,intrinsic :: iso_c_binding
 implicit none

 integer,parameter   :: nmat=10000
 integer,parameter   :: nb=64

#ifdef SINGLE_PREC
 real(4),allocatable :: matrix_a(:,:),matrix_b(:,:),matrix_c(:,:)
#else
 real(8),allocatable :: matrix_a(:,:),matrix_b(:,:),matrix_c(:,:)
#endif
 real(8) :: time
 integer :: imat,jmat
 integer :: info
!SCALAPACK
 integer            :: count_rate,count_max,count_start,count_end
 integer            :: mlocal,nlocal
 integer            :: ilocal,jlocal
 integer            :: cntxt
 integer            :: iproc,nproc
 integer            :: nprow,npcol
 integer            :: iprow,ipcol
 integer            :: desc(9)
 integer            :: neigval,neigvec
 integer,allocatable :: iwork(:)
 integer             :: liwork,lwork
#if defined(SCALAPACK)
 integer,external :: NUMROC,INDXL2G
#endif
!========================

 call system_clock(COUNT_RATE=count_rate,COUNT_MAX=count_max)

#ifdef SCALAPACK
 call MPI_INIT(info)
 call BLACS_PINFO( iproc, nproc )

 write(*,*) 'Hello world',iproc

 nprow = FLOOR(SQRT(DBLE(nproc)))
 npcol = nproc / nprow

 do while( nprow * npcol /= nproc )
   write(*,*) 'test',nprow, ' x ' , npcol,' != ',nproc
   nprow = nprow - 1
   npcol = nproc / nprow
 enddo

 if( nprow * npcol /= nproc ) stop 'change the number of procs'

 if( iproc /= 0 ) then
   close(OUTPUT_UNIT)
   open(unit=OUTPUT_UNIT,file='/dev/null')
 endif

 call BLACS_GET( -1, 0, cntxt )
 call BLACS_GRIDINIT( cntxt, 'R', nprow, npcol )
 call BLACS_GRIDINFO( cntxt, nprow, npcol, iprow, ipcol )

 write(*,*) 'Number of MPI threads:',nproc
 write(*,*) 'SCALAPACK GRID',nprow,' x ',npcol

 mlocal = NUMROC(nmat,nb,iprow,0,nprow)
 nlocal = NUMROC(nmat,nb,ipcol,0,npcol)

 call DESCINIT(desc,nmat,nmat,nb,nb,0,0,cntxt,MAX(1,mlocal),info)

#else
 write(*,*) 'Hello world no SCALAPACK'

 mlocal = nmat
 nlocal = nmat
#endif

 write(*,*) 'Global matrix:',  nmat,' x ',nmat
 write(*,*) ' Local matrix:',mlocal,' x ',nlocal


 allocate(matrix_a(mlocal,nlocal))
 allocate(matrix_b(mlocal,nlocal))
 allocate(matrix_c(mlocal,nlocal))

 !
 ! Setup the matrix
 !
#ifndef SCALAPACK
 do jmat=1,nmat
   do imat=1,nmat
     matrix_a(imat,jmat) = matrix_element1(imat,jmat)
     matrix_b(imat,jmat) = matrix_element2(imat,jmat)
   enddo
 enddo
#else
 do jlocal=1,nlocal
   do ilocal=1,mlocal
     imat = INDXL2G(ilocal,nb,iprow,0,nprow)
     jmat = INDXL2G(jlocal,nb,ipcol,0,npcol)

     matrix_a(ilocal,jlocal) = matrix_element1(imat,jmat)
     matrix_b(ilocal,jlocal) = matrix_element2(imat,jmat)
   enddo
 enddo
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Matrix timer starts here !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call system_clock(COUNT=count_start)
#ifndef SCALAPACK
#ifdef SINGLE_PREC
 call SGEMM('T','N',nmat,nmat,nmat,1.0,matrix_a,nmat,matrix_b,nmat,0.0,matrix_c,nmat)
#else
 call DGEMM('T','N',nmat,nmat,nmat,1.0d0,matrix_a,nmat,matrix_b,nmat,0.0d0,matrix_c,nmat)
#endif
#else
 call PDGEMM('T','N',nmat,nmat,nmat,1.0d0,matrix_a,1,1,desc,matrix_b,1,1,desc,0.0d0,matrix_c,1,1,desc)
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! Matrix timer stops here !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call system_clock(COUNT=count_end)

 write(*,*) 'resultat',matrix_c(1,1),matrix_c(2,1)
 deallocate(matrix_a,matrix_b,matrix_c)


 time = MODULO( count_end - count_start , count_max) / DBLE(count_rate)
 write(*,*)
 write(*,*) 'Estimated number of Floating Point operations',2.0d0*DBLE(nmat)**3
 write(*,*) 'Time (s) =',time
 write(*,*) 'GFLOPS',2.0d0*DBLE(nmat)**3 / time / 1.0d9
 write(*,*) 'End'

#ifdef SCALAPACK
 call BLACS_GRIDEXIT(cntxt)
 call BLACS_EXIT(0)
#endif

contains

real(8) function matrix_element1(i,j) result(m_ij)
 integer,intent(in) :: i,j
!===== 
 
 m_ij = SQRT(DBLE(i))/DBLE(j)

end function matrix_element1

real(8) function matrix_element2(i,j) result(m_ij)
 integer,intent(in) :: i,j
!===== 

 m_ij = 3.14159d0 * DBLE(i)**0.11 *  SQRT(DBLE(j))

end function matrix_element2


end program
