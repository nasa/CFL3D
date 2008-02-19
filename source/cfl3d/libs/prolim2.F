      subroutine prolim2(n,x1,x2,xc,leq)
c
c     $Id$
c
c***********************************************************************
c     Purpose: to limit static pressure and static density to be
c              greater than or equal to a minimum value set via the
c              data statement. For points that have the pressure/density
c              so limited, the difference stencil is also make locally
c              first order accurate.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension x1(n),x2(n),xc(n)
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10(7)
c
      data pmin,rhomin /0.01,0.01/
c
      if (leq.eq.1) then
         do 100 izz = 1,n
         xc(izz) = ccmax(xc(izz),rhomin)
         if (real(xc(izz)).eq.real(rhomin)) then
            x1(izz) = 0.
            x2(izz) = 0.
c           write(6,*)'in PROLIM2: rho,izz = ',real(xc(izz)),izz
         end if
100      continue
      else if (leq.eq.5) then
         do 200 izz = 1,n
         xc(izz) = ccmax(xc(izz),pmin)
         if (real(xc(izz)).eq.real(pmin)) then
            x1(izz) = 0.
            x2(izz) = 0.
c           write(6,*)'in PROLIM2: p,izz = ',real(xc(izz)),izz
         end if
200      continue
      end if
c
      return
      end
