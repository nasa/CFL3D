      subroutine newalpha(ncycmax,rms,clw,myid,myhost,mycomm)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  determine a new estimate for alpha that will yield
c     the specified value of Cl
c
c     original coding courtesy Neal Frink, cfl3d mods by R. Biedron
c     more mods by C. Rumsey, 3/01
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
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
      dimension istat(MPI_STATUS_SIZE)
c
#endif
      dimension rms(ncycmax),clw(ncycmax)
c
      common /alphait/ cltarg,resupdt,cltol,dalim,dcldal,alphalast,
     .                 ialphit,nttlast,icycupdt
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /conversion/ radtodeg
c
c     cltol....stop alpha iteration when abs(delcl) < cltol
c     dalim....limit of alpha change (deg.) per update
c     ioalph...if > 0, output alpha history to cfl3d.alpha
c
      cltol  = 1.e-6
      ioalph = 1
c
      if (myid.eq.myhost) then
c
         clntt = clw(ntt)
         delcl = cltarg - clntt
c     compute d(cl)/d(alpha)
         dclprv = clw(ntt) - clw(nttlast)
         daprev = alpha - alphalast
         if (daprev .ne. 0.) then
           claprv = dclprv/daprev
         else
           claprv = dcldal
         end if
         dcldal0 = 0.
c      if computed d(cl)/d(alpha) of opposite sign from
c      input value, then convergence probably not good
c      enough to use it... default to input value
         if ((real(claprv))*(real(dcldal)) .lt. 0.) claprv = dcldal
c
         if (real(rms(ntt)).gt.0) then
            rms0 = rms(nttlast)
            if(real(rms0) .eq. 0.) rms0 = 1.0
            alogres = cclog10(rms(ntt)/rms0)
         else
            alogres = 1.0
         end if
c
         if (real(cltarg).ne.0. .and.
     .       abs(real(delcl)) .gt. real(cltol)) then
c      check if residual has dropped enough
         if(real(alogres).le.real(resupdt)) then
           iupdt=1
         else
           iupdt=0
         end if
c      check if cycle number for update has been reached
c      (this takes precedence over residual update method)
         if(icycupdt .gt. 0) then
           if(mod(ntt,icycupdt) .eq. 0) then
             iupdt=1
           else
             iupdt=0
           end if
         end if
         if (iupdt .eq. 1) then
c      use computed d(cl)/d(alpha) in certain circumstances
c      (this technique copied from OVERFLOW):
             if (real(claprv) .ne. 0. .and. real(dcldal) .ne. 0.) then
                if (max(real(dcldal)/(real(claprv)),
     .                  real(claprv)/(real(dcldal))) .le. 2.) then
                  dcldal0 = claprv
                else
                  dcldal0 = dcldal
                end if
             end if
             if (real(dcldal0) .eq. 0.) then
                dal = 0.
             else
                dal = radtodeg*(cltarg-clntt)/dcldal0
             end if
c      update alpha (limited to be no greater than dalim)
             alphalast = alpha
             nttlast = ntt + 1
             if (abs(real(dal)).gt.real(dalim)) then
                dal = dalim*dal/ccabs(dal)
             end if
             alpha = alpha + dal/radtodeg
         end if
         end if
c
      end if
c
#if defined DIST_MPI
      call mpi_Bcast(alpha,1,MY_MPI_REAL,myhost,mycomm,ierr)
c
#endif
c     reset freestream for new alpha
      call init_mast
c
      if (myid.eq.myhost) then
         if (ioalph .gt. 0) then
            write(27,'(3x,i6,5e16.8)') ntt,log10(real(rms(ntt))),
     .                                 real(clw(ntt)),
     .                                 real(alpha)*real(radtodeg),
     .                                 real(dcldal0),real(claprv)
         end if
      end if
c
      return
      end
