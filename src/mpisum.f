!     this subrotine summarizes the values 
!     this is only used for mpi 
      
      subroutine  mpisum(tflx, bflx, dflx, 
     &     tpfd, bpfd, dpfd, scmpf, scmpp, rflx, rbflx, rdflx)

      implicit none

      include 'mpif.h'
      include 'common.inc'
      include 'math.inc'

      integer ierr
      integer spsize,brfsize,refsize,parsize
      real tflx, bflx, dflx
      real tpfd, bpfd, dpfd
      real scmpf(3,100), scmpp(3100)
      real rflx, rbflx, rdflx

!     flux 
      call MPI_Reduce(tflx,tflx,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(bflx,bflx,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dflx,dflx,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tpfd,tpfd,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(bpfd,bpfd,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(dpfd,dpfd,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(rflx,rflx,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(rbflx,rbflx,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(rdflx,rdflx,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
     
!     flux spectrum
      spsize = 3 * 100
      call MPI_Reduce(scmpf,scmpf,spsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(scmpp,scmpp,spsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 


!     brf 
      brfsize = 2 * 700
      call MPI_Reduce(brf,brf,brfsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(brfc,brfc,brfsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(brfs,brfs,brfsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)      
      call MPI_Reduce(brff,brff,brfsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 

!     ref images
      refsize = 2 * size * size
      call MPI_Reduce(refl,refl,refsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(irefl,irefl,refsize,MPI_INT,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 

!     par 

      call MPI_Reduce(cfpr,cfpr,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(bfpr,bfpr,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(ffpr,ffpr,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(sfpr,sfpr,1,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 

      parsize = size * size
      call MPI_Reduce(apf,apf,parsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr) 
      call MPI_Reduce(aps,aps,parsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(apb,apb,parsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(apfd,apfd,parsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)

      parsize = parsize * 100
      call MPI_Reduce(ap,ap,parsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)
      call MPI_Reduce(apd,apd,parsize,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)

      call MPI_Reduce(apnp,apnp,100,MPI_REAL,MPI_SUM,
     &     0,MPI_COMM_WORLD,ierr)      

      return
      end
