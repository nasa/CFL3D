      subroutine lesdiag(myid,jdim,kdim,idim,q,ux,vist3d,vol,si,sj,sk,
     +                   nou,bou,nbuf,ibufdim)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute diagnostics for LES-type runs.  If a 
c     Smagorinsky-type SGS model is being used, the appropriate eddy
c     viscosity is computed and set here.
c     Output goes to files fort.50x, depending on myid (i.e., if 8
c     processors are used, files fort.501-508 will be written).
c
c     Currently coded: standard Smagorinsky model, no wall damping on Cs
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /reyue/ reue,tinf,ivisc(3)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /sklton/ isklton
      common /lesinfo/ cs_smagorinsky
c
      dimension q(jdim,kdim,idim,5)
      dimension ux(jdim-1,kdim-1,idim-1,9)
      dimension vist3d(jdim,kdim,idim)
      dimension sj(jdim,kdim,idim-1,5),
     + sk(jdim,kdim,idim-1,5),si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
      imodel=1
      cs=cs_smagorinsky
      if(isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Computing LES-type turbulent'',
     +'' viscosity, block='',i5)') nbl
         if (imodel .eq. 0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''        no model'')')
         else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''        Smag model, Cs='',f10.5)') cs
         end if
      end if
c
c   acquire mu_t if employing an SGS model
c
      if (imodel .eq. 0) then
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            vist3d(j,k,i)=0.
          enddo
        enddo
      enddo
      else
      do i=1,idim-1
        do j=1,jdim-1
          do k=1,kdim-1
            deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
            deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
            deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
            if( i2d .ne. 1 ) then
               delta = (deltai*deltaj*deltak)**0.333333
            else
               delta = sqrt(deltaj*deltak)
            end if
            s11 = ux(j,k,i,1)
            s22 = ux(j,k,i,5)
            s33 = ux(j,k,i,9)
            s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
            s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
            s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
            xis = s11*s11 + s22*s22 + s33*s33 +
     +            2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
c           Standard Smagorinsky model:
            vist3d(j,k,i)=q(j,k,i,1)*(cs*delta)**2*sqrt(2.*xis)*
     +            reue/xmach
          enddo
        enddo
      enddo
      end if
c
c   obtain LES diagnostic info
c
      if (icyc .eq. ncyc) then
      c2b=cbar/tinf
      c2bp=c2b+1.0
      eddytime=1./xmach
      sumvel2=0.
      sumux2=0.
      sumvy2=0.
      sumwz2=0.
      sumux3=0.
      sumvy3=0.
      sumwz3=0.
      sumu2=0.
      sumv2=0.
      sumw2=0.
      sumuv=0.
      sumuw=0.
      sumvw=0.
      sumenergy=0.
      fnu=0.
      sumxis=0.
      sumwis=0.
      rho=0.
      n=0
      do i=1,idim-1
      do j=1,jdim-1
      do k=1,kdim-1
        n=n+1
        sumvel2=sumvel2+q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2
        sumux2=sumux2+(ux(j,k,i,1)**2)
        sumvy2=sumvy2+(ux(j,k,i,5)**2)
        sumwz2=sumwz2+(ux(j,k,i,9)**2)
        sumux3=sumux3+(ux(j,k,i,1)**3)
        sumvy3=sumvy3+(ux(j,k,i,5)**3)
        sumwz3=sumwz3+(ux(j,k,i,9)**3)
        sumu2=sumu2+q(j,k,i,2)**2
        sumv2=sumv2+q(j,k,i,3)**2
        sumw2=sumw2+q(j,k,i,4)**2
        sumuv=sumuv+q(j,k,i,2)*q(j,k,i,3)
        sumuw=sumuw+q(j,k,i,2)*q(j,k,i,4)
        sumvw=sumvw+q(j,k,i,3)*q(j,k,i,4)
        sumenergy=sumenergy+q(j,k,i,5)
        s11 = ux(j,k,i,1)
        s22 = ux(j,k,i,5)
        s33 = ux(j,k,i,9)
        s12 = 0.5*(ux(j,k,i,2) + ux(j,k,i,4))
        s13 = 0.5*(ux(j,k,i,3) + ux(j,k,i,7))
        s23 = 0.5*(ux(j,k,i,6) + ux(j,k,i,8))
        w12 = 0.5*(ux(j,k,i,2) - ux(j,k,i,4))
        w13 = 0.5*(ux(j,k,i,3) - ux(j,k,i,7))
        w23 = 0.5*(ux(j,k,i,6) - ux(j,k,i,8))
        xis = s11*s11 + s22*s22 + s33*s33 +
     +        2.*s12*s12 + 2.*s13*s13 + 2.*s23*s23
        wis = 2.*w12*w12 + 2.*w13*w13 + 2.*w23*w23
        tt=gamma*q(j,k,i,5)/q(j,k,i,1)
        fnu=fnu+(c2bp*tt*sqrt(tt)/(c2b+tt))
        sumxis=sumxis+xis
        sumwis=sumwis+wis
        rho=rho+q(j,k,i,1)
      enddo
      enddo
      enddo
      avgux2=sumux2/float(n)
      avgvy2=sumvy2/float(n)
      avgwz2=sumwz2/float(n)
      avgux3=sumux3/float(n)
      avgvy3=sumvy3/float(n)
      avgwz3=sumwz3/float(n)
      avgvel2=sumvel2/float(n)
      avgu2=sumu2/float(n)
      avgv2=sumv2/float(n)
      avgw2=sumw2/float(n)
      avguv=sumuv/float(n)
      avguw=sumuw/float(n)
      avgvw=sumvw/float(n)
      avgfnu=fnu/float(n)
      avgrho=rho/float(n)
      avgwis=sumwis/float(n)
      xtime=time/eddytime
c   for each processor, write out:
c   ntt,icyc,xtime,sumvel2,sumxis,avgux2,avgvy2,avgwz2,
c   avgux3,avgvy3,avgwz3,avgvel2,avgu2,avgv2,avgw2,avgfnu,avgrho,avgwis
      write(500+myid,'(2i5,20e15.5)') ntt,icyc,xtime,sumvel2,sumxis,
     +   avgux2,avgvy2,avgwz2,avgux3,avgvy3,avgwz3,avgvel2,
     +   avgu2,avgv2,avgw2,avgfnu,avgrho,avgwis,sumenergy,
     +   avguv,avguw,avgvw
c
      end if
      return
      end
