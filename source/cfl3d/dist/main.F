      program main
c
c     $Id$
c
c***********************************************************************
c     CFL3D :  Characteristic-based scheme for steady/unsteady solutions
c              to the Euler/Navier-Stokes equations.
c
c     VERSION 6.0 :  Computational Fluids Laboratory, Mail Stop 128,
c                    NASA Langley Research Center, Hampton, VA
c
c
c***********************************************************************
c
c     The code CFL3D has evolved from the contributions of a number of
c     researchers at NASA Langley, principally from the Computational 
c     Fluids Lab, including AAMB, UAB, and ICASE. Points of contact for
c     particular expertise in certain features of the code are listed 
c     below.
c     
c     
c     POINTS OF CONTACT AND AREAS OF PARTICULAR EXPERTISE:
c
c     Kyle Anderson        NASA Langley AAMB            (757)864-2164
c       ** multigrid algorithm development
c          and unsteady flows
c
c     Bob Biedron          NASA Langley AAMB            (757)864-2156
c       ** multiblock applications; patched-grid
c          and chimera-grid preprocessors;
c          unsteady flows/dynamic grid applications;
c          parallel implementation; complex-variable
c          implementation; aeroelastic implementation
c
c     Farhad Ghaffari      NASA Langley TAB             (757)864-2856
c       ** transonic applications
c
c     Chris Rumsey         NASA Langley AAMB            (757)864-2165
c       ** viscous steady/unsteady flows, turbulence
c          models, and CGNS implementation
c
c     Jim Thomas           NASA Langley AAMB            (757)864-2163
c       ** algorithm development and applications
c
c     Veer Vatsa           NASA Langley AAMB            (757)864-2236
c       ** transonic viscous flows
c
c     Bob Walters          VPI&SU                       (703)961-6748
c       ** algorithm development and applications,
c          related PNS and equilibrium air codes,
c          and generalized chemistry codes (GASP)
c
c     Oktay Baysal         ODU                          (757)864-5570
c       ** overlapped grids/chimera scheme
c
c     Jerry Mall           CSC                          (757)766-8223
c       ** parallel implementation
c
c     Jack Edwards, NCSU                                (919)515-5264
c       ** low mach number preconditioning
c
c     Bob Bartels          NASA Langley AB              (757)864-2813
c       ** aeroelastic applications, spring analogy
c          mesh deformation scheme
c
c***********************************************************************
c
c      Brief description of major code features/attributes : 
c
c       (1) Finite volume discretization with q= (rho,u,v,w,p) at cell 
c           centers.
c
c       (2) Upwind-biased convective/pressure term differencing using
c           either flux-vector-splitting of Van Leer or flux-difference-
c           splitting of Roe.
c
c       (3) Central diferencing of thin-layer Navier-Stokes terms with
c           turbulence accounted for via zero-equation, one-equation or 
c           two-equation models.
c
c       (4) Spatially-factored 3-factor implicit algorithm with either
c           5x5 block inversion or diagonal inversions applicable to
c           both steady and unsteady flows.
c
c       (5) Subiterations for reduction of linearization and 
c           factorization errors on the left hand side
c
c       (6) FAS multigrid acceleration for both steady and unsteady
c           flows.
c
c       (7) Local solution refinement via embedded meshes
c
c       (8) Domain decomposition using  C-0 continuous interfaces, 
c           patched interfaces, and overlapped/embedded grids through
c           the chimera scheme.
c    
c       (9) Boundary conditions can be set over subsets of block faces 
c
c      (10) Grid motion via translation or rotation; blocks may move in
c           relative motion past one another (patched grids only)
c
c      (11) Parallel processing via MPI protocals; each processor may
c           have one or more grid blocks assigned to it.
c
c*******IMPORTANT*******IMPORTANT*******IMPORTANT*******IMPORTANT*******
c
c  If you discover any bugs, please let us know so that they may be
c  corrected
c
c***********************************************************************
c
c     Description of the files read/written.
c
c     unit  1 - Input binary grid file.
c     unit  2 - Input/output file containing binary restart data.
c     unit  3 - Output binary grid file for PLOT3D.
c     unit  4 - Output binary dependent variable file for PLOT3D.
c     unit  5 - Input case data file.
c     unit  8 - Parameter files
c     unit 11 - Primary output file.
c     unit 12 - Output file for convergence history (residual, forces,
c               moments)
c     unit 15 - Secondary output file (programmer's output).
c     unit 17 - Output file for flowfield printing.
c     unit 20 - Output file of unsteady pressure coeffs (dynamic
c               grids only).
c     unit 21 - Input file for overlapped grid information,
c     unit 22 - Input file for patched grid information.
c
c     Note: Units 1, and 5 must be available at beginning of program
c           execution (i.e. the files attached to them must exist and
c           contain appropriate data). Unit 2 is the binary restart
c           file and is overwritten during program execution; it needs
c           to be available only when restarting. Unit 21 must be
c           available if using the overlapped grid option (iovrlp > 0).
c           Unit 22 must be available if using the (static) patched grid
c           option (ninter < 0). Dynamic (i.e. vary with time) pathced
c           interfaces are handled via the main input deck.
c
c           File names attached to the following unit numbers are not 
c           set by the user, but are opened for output with hardwired
c           names:
c           23......cfl3d.subit_res
c           24......cfl3d.subit_turres
c           25......cfl3d.dynamic_patch
c           27......cfl3d.alpha
c           28-31...(opened in subroutine yplusout only if ifunc is
c                   set > 0 in that subroutine)
c           32......tempz3y2x1.inp (copy of input file, deleted at normal 
c                   termination)
c           33......aesurf.dat
c           34......genforce.dat
c           35......cfl3d.sd_res
c           39......newsurf.p3d
c           66......precfl3d.out
c           96......cfl3d_avgg.p3d
c           97......cfl3d_avgq.p3d
c           99......precfl3d.error / cfl3d.error
c
c           Unit 26 is reserved for reading in bc data via the
c           2000 series of boundary conditions.
c
c***********************************************************************
c
c     Principle variables required for distributed computing.
c
c     myhost  - host process ID
c     myid    - local process ID
c     mycomm  - MPI's global group communicator
c     nnodes  - number of node processes (myhost not included)
c     mblk2nd - assigment of blocks to nodes
c
c***********************************************************************
c
c     DEVELOPER NOTES:
c
c     1) Subroutine findmin (module lbcx, used to find the minimum
c        distance from a surface for use with 1 and 2 eq. turb.
c        models), key off boundary condition type 2004 (viscous
c        surface). At the current time, this is the only BC type in
c        the standard suite of BC types that makes sense to include
c        in findmin_new. However, in future, if other viscous surface 
c        BC's are added to the suite, then those BC types will also
c        have to be included in findmin_new.
c        NOTE: the old, slow method for finding the minimum distance
c              has been eliminated from this version, and findmin_new
c              has been extended to include the extra data needed for 
c              then Baldwin-Barth turb. model.
c
c     2) Subroutine metric (module lbcx) NO LONGER keys off various
c        singular axis BC's - detection of collapsed metrics is automatic
c        Detection of collapsed metrics keys off a parameter atol, which
c        may need to be adjusted in some cases (but to date has proven
c        to be robust). "Trust but verify"
c
c     3) Storage space has been provided to store the acceleration terms
c        for points on block boundaries. However, at this time, the 
c        acceleration terms are set to zero in subroutines trans and 
c        rotate. For some high rotation rate problems, it may in be
c        necessary to add the acceleration terms to the normal 
c        pressure gradient at solid walls.
c
c     4) the subroutine getdhdr must be modified if any new bc's of the 
c        2000 series are added (bc's that need user-specified data). 
c        getdhdr sets up character variables that are used for headers
c        in the output file for the 2000 series bc's
c
c     5) 10.**(-iexp) is machine zero; iexp is available in common 
c        block /zero/ and may be used for setting tolerences
c
c     6) Beware the expanding MY_MPI_REAL: this variable is used in
c        MPI calls that pass real variables, but if double precision 
c        compilation is used, MY_MPI_REAL will expand to the longer
c        MPI_DOUBLE_PRECISION. Thus, the MPI calls may appear to 
c        fit within 72 characters in the fortan statement as coded,
c        but during compilation will expand to > 72 columns. The effect
c        will typically be to hang the MPI process. 
c
c     7) A dummy module, development.F is provided for testing code
c        modifications/enhancements. This module is always compiled
c        by the makefile, but as issued, just contains a short dummy
c        subroutine. Subroutines placed in development.F ((including
c        routines involving MPI), will supercede any other version of
c        the routine in the code, with the exception of the driver
c        routine (main.F).  NOTE: for code transfer, the packitup
c        script temporarily replaces the current contents of 
c        development.F with a clean version, so that the packaged
c        version has the short dummy routine. Upon completion of the
c        packitup script, the prior contents of development.F are
c        restored.
c
c     8) To help coordinate printout of the various processors, 
c        an internal write buffer (array) is used, which is periodically
c        flushed. The array, bou, is dimensioned as (ibufdim,nbuf),
c        where nbuf is the number of output files that can be supported
c        by the buffer, and ibufdim is the max number of 120 character 
c        lines that can be stored in the file. Currently supported 
c        output files are:
c
c        ibuf    unit/file
c           1    11  / main output (e.g. "cfl3d.out") 
c           2    09  / fort.9, aux. file for dynamic patching
c           3    14  / baldwin-lomax output (e.g. "cfl3d.blomax")
c           4    25  / dynamic patch output (e.g. "cfl3d.dynamic_patch")
c
c     9) The code is now "self sizing" - precfl3d is built in, along
c        with memory allocation routines employing Cray pointers. These
c        pointers work in f77 and f90 for all but character variables.
c        Therefore, the following parameters pertaining to character
c        variable arrays are hardwired below to reasonably large values
c        that should cover most cases:
c
c        ibufdim...see 8) above
c        nbuf......see 8) above
c        mxbcfil...maximum number of bc file names for 2000 series bc's
c                  i.e. up to mxbcfil bc files may be specified in the 
c                  cfl3d input deck.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
c
#     include "mpif.h"
#   ifdef DBLE_PRECSN
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_DOUBLE_COMPLEX
#      else
#        define MY_MPI_REAL MPI_DOUBLE_PRECISION
#      endif
#   else
#      ifdef CMPLX
#        define MY_MPI_REAL MPI_COMPLEX
#      else
#        define MY_MPI_REAL MPI_REAL
#      endif
#   endif
      character errmsg(MPI_MAX_ERROR_STRING)
      dimension istat(MPI_STATUS_SIZE)
#endif
c
      parameter (ibufdim=2000,nbuf=4,mxbcfil=100)
c
      external usrint
c
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      character*80 inpstring
      character*80 bcfiles(mxbcfil)
      character*120 bou(ibufdim,nbuf)
      character*32 basedesired
c
      dimension itrnsfr(27),titlw(20),nou(nbuf)
c
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
      common /proces/ numprocs
      common /mydist2/ nnodes,myhost,myid,mycomm
      common /unit5/ iunit5
      common /key/ nkey
      common /zero/ iexp
c
#if defined DIST_MPI
c***********************************************************************
c     MPI initialization
c***********************************************************************
c
      call MPI_Init (ierr)
      if (ierr.ne.MPI_SUCCESS) then
         write (11,*)'MPI_Init not successful...stopping'
         call MPI_Error_String(ierr, errmsg, len, ier)
c        write (11,a1) (errmsg(l),l=1,len)
      endif
      call MPI_Comm_Rank (MPI_COMM_WORLD, myid, ierr)
      call MPI_Comm_Size (MPI_COMM_WORLD, numprocs, ierr)
c
      mycomm     = MPI_COMM_WORLD
      myhost     = 0
      istat_size = MPI_STATUS_SIZE
      nnodes = numprocs -1
c
c     check for inappropriate number of processors
c
      if (numprocs.lt.2) then
         write (11,66)
   66    format('stopping...only running a host process; ',
     .          'no node process',/,'you must run the ',
     .          'mpi version with more than 1 process')
         call MPI_ABORT(MPI_COMM_WORLD, myid, mpierror)
         call MPI_Finalize (ierr)
      end if
c
#else
c     initialization of sequential version
c
      myhost     = 0
      myid       = 0
      mycomm     = 0
      istat_size = 1
      nnodes     = 1
#endif
c
      do ii=1,nbuf
         nou(ii) = 0
         do mm=1,ibufdim
            bou(mm,ii) = ' '
         end do
      end do 
c
c     determine machine zero for use in setting tolerances
c     (10.**(-iexp) is machine zero) 
c
      icount = 0
      compare = 1.0
      do i = 1,20
         icount = icount + 1
         add = 1.
         do n=1,i
            add = add*.1
         enddo
         x11 = compare + add
         if (x11.eq.compare)then
            iexp = i-1
            goto 4010
         end if
      end do
 4010 continue
c
c***********************************************************************
c     set some signals for receiving a SIGTERM
c***********************************************************************
c
c      8...arithmetic / floating point exception
c      9...killed
#  if defined LINUX
c      7...bus error
#  else
c     10...bus error
#  endif
c     11...segmentation fault
c
#if defined IBM
      call signal(8,usrint)
      call signal(9,usrint)
      call signal(10,usrint)
      call signal(11,usrint)
#else
      call signal(8,usrint,-1)
      call signal(9,usrint,-1)
#  if defined LINUX
      call signal(7,usrint,-1)
#  else
      call signal(10,usrint,-1)
#  endif
      call signal(11,usrint,-1)
#endif
c
c***********************************************************************
c     open files
c***********************************************************************
c
#if defined DIST_MPI
c
      if (myid.eq.myhost) then
#endif
c
      iunit5 = 32
c
#if defined NOREDIRECT
      open(iunit5,file='cfl3d.inp',form='formatted',status='old')
# else
      open(iunit5,file='tempz3y2x1.inp',form='formatted',
     .     status='unknown')
c
c     first copy the command-line input file to a temporary input file
c     this is necessary because the input file must be read more 
c     more than once to accomodate dynamic memory allocation, and
c     some operating systems do not permit a rewind of a file
c     that has not been explicitly opened
c
c     NOTE: this assumes the cfl3d input file contains lines of no
c           more than 80 characters!!!
c 
      do n=1,9999
         read(5,'(a80)',end=999) inpstring
         write(iunit5,'(a80)') inpstring
      end do
  999 continue
      rewind(iunit5)
# endif
c
      read(iunit5,*)
      read(iunit5,'(a60)')grid
      read(iunit5,'(a60)')plt3dg
      read(iunit5,'(a60)')plt3dq
      read(iunit5,'(a60)')output
      read(iunit5,'(a60)')residual
      read(iunit5,'(a60)')turbres
      read(iunit5,'(a60)')blomx
      read(iunit5,'(a60)')output2
      read(iunit5,'(a60)')printout
      read(iunit5,'(a60)')pplunge
      read(iunit5,'(a60)')ovrlap
      read(iunit5,'(a60)')patch
      read(iunit5,'(a60)')restrt
c
c     asn file calls used on cray if single-precision unformated
c     files are used/desired. Such files are directly readable
c     on workstations (e.g. in FAST, PLOT3D, etc)
c
#if defined ASN_GRD
      call asnfile(grid, '-F f77 -N ieee', IER)
#endif
#if defined ASN_P3D
      call asnfile(plt3dg, '-F f77 -N ieee', IER)
      call asnfile(plt3dq, '-F f77 -N ieee', IER)
#endif
c
      if (icgns .ne. 1) then
      open(unit=1,file=grid,form='unformatted',status='old')
      open(unit=2,file=restrt,form='unformatted',status='unknown')
      else
      write(iunit5,'(/,'' icgns=1.  Therefore, restart filename:'',/,
     . 6x,a60,/,'' is not used.  All restart info is written to cgns'',
     . '' file:'',/,6x,a60,/,'' (where grid is also located)'')') 
     . restrt,grid
      end if
      open(unit=11,file=output,form='formatted',status='unknown')
      open(unit=12,file=residual,form='formatted',status='unknown')
      open(unit=13,file=turbres,form='formatted',status='unknown')
      open(unit=14,file=blomx,form='formatted',status='unknown')
      open(unit=15,file=output2,form='formatted',status='unknown')
      open(unit=20,file=pplunge,form='formatted',status='unknown')
      open(unit=21,file=ovrlap,form='unformatted',status='unknown')
      open(unit=22,file=patch,form='unformatted',status='unknown')
c
c     the following file names are not set from the input deck
c     NOTE: unit 26 is reserved for boundary condition data files, 
c     the name(s) of which are read in via the input deck
c
      subres    = 'cfl3d.subit_res'
      subtur    = 'cfl3d.subit_turres'
      grdmov    = 'cfl3d.dynamic_patch'
      alphahist = 'cfl3d.alpha'
      preout    = 'precfl3d.out'
      errfile   = 'precfl3d.error'
      aeinp     = 'aesurf.dat'
      aeout     = 'genforce.dat'
      sdhist    = 'cfl3d.sd_res'
      avgg      = 'cfl3d_avgg.p3d'
      avgq      = 'cfl3d_avgq.p3d'
c
      open(unit=23,file=subres,form='formatted',status='unknown')
      open(unit=24,file=subtur,form='formatted',status='unknown')
      open(unit=25,file=grdmov,form='formatted',status='unknown')
      open(unit=27,file=alphahist,form='formatted',status='unknown')
      open(unit=33,file=aeinp,form='formatted',status='unknown')
      open(unit=34,file=aeout,form='formatted',status='unknown')
#   ifdef CMPLX
      open(unit=35,file=sdhist,form='formatted',status='unknown')
#   endif
      open(unit=66,file=preout,form='formatted',status='unknown')
      open(unit=96,file=avgg,form='unformatted',status='unknown')
      open(unit=97,file=avgq,form='unformatted',status='unknown')
      open(unit=99,file=errfile,form='formatted',status='unknown')
c
      rewind( 1)
      rewind( 2)
      rewind( 3)
      rewind( 4)
      rewind(11)
      rewind(12)
      rewind(13)
      rewind(14)
      rewind(15)
      rewind(17)
      rewind(20)
      rewind(21)
      rewind(22)
      rewind(23)
      rewind(24)
      rewind(25)
      rewind(27)
      rewind(33)
      rewind(34)
      rewind(66)
      rewind(99)
c
c***********************************************************************
c     determine array size requirements
c***********************************************************************
c
c     read input file to get the array dimensions needed by precfl3d
c
      rewind(iunit5)
c
      ibufdim0 = ibufdim
      nbuf0    = nbuf
      mxbcfil0 = mxbcfil
c
c     global0 sets the parameters needed for precfl3d (sizer)
c
      iunit11 = 99
      call global0(nplots0,maxnode0,mxbli0,lbcprd0,lbcemb0,
     .             lbcrad0,maxbl0,maxgr0,maxseg0,maxcs0,ncycmax0,
     .             intmax0,nsub10,intmx0,mxxe0,mptch0,msub10,
     .             ibufdim0,nbuf0,mxbcfil0,nmds0,maxaes0,
     .             maxsegdg0,ntr,nnodes,nou,bou,iunit11,myid,
     .             idm0,jdm0,kdm0)
c
c     rewind the input file so it can be read again
c
      rewind(iunit5)
c
c     sizer sets the parameters needed for cfl3d
c
      imode = 1
      mwork = 1
      call sizer(mwork,mworki,nplots,minnode,iitot,intmax,
     .           maxxe,mxbli,nsub1,lbcprd,lbcemb,lbcrad,
     .           maxbl,maxgr,maxseg,maxcs,ncycmax,intmx,mxxe,
     .           mptch,msub1,nmds,maxaes,maxsegdg,nnodes,nslave,
     .           nmaster,myhost,myid,mycomm,nplots0,maxnode0,mxbli0,
     .           lbcprd0,lbcemb0,lbcrad0,maxbl0,maxgr0,maxseg0,
     .           maxcs0,ncycmax0,intmax0,nsub10,intmx0,mxxe0,
     .           mptch0,msub10,ibufdim0,nbuf0,mxbcfil0,nmds0,
     .           maxaes0,maxsegdg0,imode,ntr,bcfiles,bou,nou)
c
#if defined DIST_MPI
c     send the array dimensions to all processors
c
      itrnsfr(1 ) = mwork
      itrnsfr(2 ) = mworki
      itrnsfr(3 ) = maxgr
      itrnsfr(4 ) = maxbl
      itrnsfr(5 ) = minnode
      itrnsfr(6 ) = maxseg
      itrnsfr(7 ) = nplots
      itrnsfr(8 ) = ncycmax
      itrnsfr(9 ) = mxbli
      itrnsfr(10) = intmax
      itrnsfr(11) = intmx
      itrnsfr(12) = maxcs
      itrnsfr(13) = mxxe
      itrnsfr(14) = mptch
      itrnsfr(15) = nsub1
      itrnsfr(16) = msub1
      itrnsfr(17) = lbcprd
      itrnsfr(18) = lbcemb
      itrnsfr(19) = lbcrad
      itrnsfr(20) = maxxe
      itrnsfr(21) = iitot
      itrnsfr(22) = nmds
      itrnsfr(23) = maxaes
      itrnsfr(24) = nslave
      itrnsfr(25) = maxsegdg
      itrnsfr(26) = nkey
      itrnsfr(27) = nmaster
c
      do inode = 1,nnodes
         itag = inode
         call MPI_Send (itrnsfr, 27, MPI_INTEGER, inode,
     .                 itag, mycomm, ierr)
      end do
c
      else
c
      itag = myid
      call MPI_Recv (itrnsfr, 27, MPI_INTEGER, myhost,
     .               itag, mycomm, istat, ierr)
c
      mwork    = itrnsfr(1 )
      mworki   = itrnsfr(2 )
      maxgr    = itrnsfr(3 )
      maxbl    = itrnsfr(4 )
      minnode  = itrnsfr(5 )
      maxseg   = itrnsfr(6 )
      nplots   = itrnsfr(7 )
      ncycmax  = itrnsfr(8 )
      mxbli    = itrnsfr(9 )
      intmax   = itrnsfr(10)
      intmx    = itrnsfr(11)
      maxcs    = itrnsfr(12)
      mxxe     = itrnsfr(13)
      mptch    = itrnsfr(14)
      nsub1    = itrnsfr(15)
      msub1    = itrnsfr(16)
      lbcprd   = itrnsfr(17)
      lbcemb   = itrnsfr(18)
      lbcrad   = itrnsfr(19)
      maxxe    = itrnsfr(20)
      iitot    = itrnsfr(21)
      nmds     = itrnsfr(22)
      maxaes   = itrnsfr(23)
      nslave   = itrnsfr(24)
      maxsegdg = itrnsfr(25)
      nkey     = itrnsfr(26)
      nmaster  = itrnsfr(27)
c
      end if
#endif
c
c     close/delete precfl3d error file, and reopen as cfl3d error file
c
      close(99,status='delete')
      errfile   = 'cfl3d.error'
      open(unit=99,file=errfile,form='formatted',status='unknown')
c
c***********************************************************************
c     allocate memory and run flow solver
c***********************************************************************
c
      call cfl3d(mwork,mworki,nplots,minnode,iitot,intmax,
     .           maxxe,mxbli,nsub1,lbcprd,lbcemb,lbcrad,maxbl,
     .           maxgr,maxseg,maxcs,ncycmax,intmx,mxxe,mptch,
     .           msub1,ibufdim,nbuf,mxbcfil,istat_size,ntr,
     .           nmds,maxaes,nslave,maxsegdg,nmaster,bcfiles,bou,nou)
c
c***********************************************************************
c     normal program termnination
c***********************************************************************
c   Write input file used to CGNS file
#if defined CGNS
      if (icgns .eq. 1 .and. myid .eq. myhost) then
        basedesired='Base'
        idimdesired=3
        call wopencgns(grid,basedesired,idimdesired,iccg,
     .              ibase,nzones)
        call writeinput(iccg,ibase,iunit5)
c   Close data base
        call cg_close_f(iccg, ier)
        if (ier .ne. 0) then
          call cg_error_print_f
          call termn8(myid,-1,ibufdim,nbuf,bou,nou)
        end if
      end if
#endif
c
c     remove the temporary input file, if used
c
#if defined NOREDIRECT
#else
      if (myid.eq.myhost) then
         close(iunit5,status='delete')
      end if
#endif
c
      call termn8(myid,0,ibufdim,nbuf,bou,nou)
c
      stop
      end
