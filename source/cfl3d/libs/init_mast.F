      subroutine init_mast         
c***********************************************************************
c     Purpose:  Initialize the initial conditions for use on the master
c     (host) node
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10,tur20
c
      rho0   = 1.e0
      c0     = 1.e0
      p0     = rho0*c0*c0/gamma
c
c   The wind axis system follows NASA SP-3070 (1972), with the exception that
c   positive beta is in the opposite direction
c
      u0     = xmach*cos(alpha)*cos(beta)
      w0     = xmach*sin(alpha)*cos(beta)
      v0     = -xmach*sin(beta)
      ei0    = p0/((gamma-1.e0)*rho0)
      et0    = rho0*(ei0+.5e0*(u0*u0+v0*v0+w0*w0))
      h0     = (et0+p0)/rho0
      pt0    = p0*(1.e0+.5e0*gm1*xmach*xmach)**(gamma/gm1)
      qiv(1) = rho0
      qiv(2) = u0
      qiv(3) = v0
      qiv(4) = w0
      qiv(5) = p0
c
      return
      end
