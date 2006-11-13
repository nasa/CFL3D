      subroutine rb_pred(nt)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Update the rigid body displacements and velocities via
c              a predictor step of the rigid body equations of motion.
c
c              
c        Reference: Vinh, Lam-Son, Edwards, J.W., Seidel, D. A., Batina,
c                   J. T., "Transonic Stability and Control of Aircraft
c                   Using CFD Methods," AIAA Paper 88-4374.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /rbstmt1/ bmatrb(4,4),gforcnrb(4),gforcnmrb(4),gforcsrb(4),
     .                 stmrb(4,4),stmirb(4,4),xsrb(4),xxnrb(4),xnmrb(4),
     .                 x0rb(4)
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat 
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4),zlfct,epstr,relax,ittrst
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
c
c     rigid body displacement and velocity prediction via rigid body
c     equations of motion
c
      do j = 1,4
         gforcnmrb(j) = gforcnrb(j)
         gforcnrb(j)  = gforcsrb(j)
      end do
c
      do n=1,4
         xxnrb(n) = xsrb(n)
      end do
c
      do n=1,4
         xsrb(n) = 0.
         do j=1,4
            xsrb(n) = xsrb(n) + stmrb(n,j)*xxnrb(j)
     .              + .5*stmirb(n,j)*(3.*gforcnrb(j) - gforcnmrb(j))
         end do
      end do
c
      zrg1 = zrg0 - xsrb(1)/greflrb
      dzdt = - xsrb(2)/greflrb
      alf1 =  xsrb(3)
      thtd1=  xsrb(4)
      write(299,31029) dt*real(nt),zrg0-xsrb(1),xsrb(3),gforcnrb(2),
     .                 gforcnrb(4)
31029 format(5(1x,e16.8))
c
      return
      end
