      subroutine init_trim
c
c     $Id$
c
c***********************************************************************
c     Purpose: initialize trim data.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /trim/ dmtrmn,dmtrmnm,dlcln,dlclnm,trtol,cmy,cnw,alf0,
     .              alf1,dzdt,thtd0,thtd1,zrg0,zrg1,dtrmsmx,dtrmsmn,
     .              dalfmx,ddtmx,ddtrm0,ddtrm1,itrmt,itrminc,fp(4,4),
     .              tp(4,4)
c
c     The following relaxation coefficients and limits are hardwired:
c
      itrminc = 5
      dclda   = 6.
      dcldd   = 1.4
      dcmda   =-.2
      dcmdd   =-.88
      dtrmsmx =  0.50
      dtrmsmn = -0.50
      dmtrmn  =  0.
      dmtrmnm =  0.
      dlcln   =  0.
      dlclnm  =  0.
      ddtmx   =  0.02
      dalfmx  =  0.02
      trtol   =  0.0002
      itrmt   =  0
c
      tp(1,1) = dclda
      tp(1,2) = dcldd
      tp(2,1) = dcmda
      tp(2,2) = dcmdd
      dtr     = tp(1,1)*tp(2,2)-tp(1,2)*tp(2,1)
      fp(1,1) = tp(2,2)/dtr
      fp(1,2) =-tp(1,2)/dtr
      fp(2,1) =-tp(2,1)/dtr
      fp(2,2) = tp(1,1)/dtr
c
      return
      end
