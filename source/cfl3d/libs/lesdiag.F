      subroutine lesdiag(myid,jdim,kdim,idim,q,ux,vist3d,vol,si,sj,sk,
     +                   vor,smin,xjb,tursav,xkb,blnum,
     +                   nou,bou,nbuf,ibufdim,nbl,nummem,x,y,z)
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
c     Currently coded: standard Smagorinsky model
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
      dimension vor(jdim-1,kdim-1,idim-1),smin(jdim-1,kdim-1,idim-1)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
c     the screwy method for keeping track of nearest wall location (for
c     getting y+ value associated with a field point) is taken from method
c     used for Baldwin-Barth... it's kind of a memory nightmare:
c     xjb=j index, tursav(2)=i index, xkb=k index, blnum=corresp block number
c     Actually, the Baldwin-Barth index tracking appears to be screwed up -
c     so for now, do not use it!!
c     i.e., xjb, xkb, tursav(2), and blnum are not used
      dimension xjb(jdim-1,kdim-1,idim-1),xkb(jdim-1,kdim-1,idim-1),
     + tursav(jdim,kdim,idim,nummem),blnum(jdim-1,kdim-1,idim-1)
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
c
c     imodel controls whether or not a model is used (0=no model)
c     (for Smagorinsky, this can also be achieved by setting cs_smagorinsky=0)
      imodel=1
c     iwallscale controls whether the Delta is scaled by van Driest near walls
      iwallscale=1
c
      cs=cs_smagorinsky
      c2b=cbar/tinf
      c2bp=c2b+1.0
      re=reue/xmach
c
      if(isklton .gt. 0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     Computing LES-type turbulent'',
     +'' viscosity, block='',i5)') nbl
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'(''     WARNING: ivisc=25 still under'',
     .    '' development... use at your own risk!'')')
         if (imodel .eq. 0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''        no model'')')
         else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''        Smag model, Cs='',f10.5)') cs
         end if
         if (iwallscale .eq. 0) then
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''        no wall damping'')')
         else
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'(''        van Driest wall damping'')')
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
c   Set ypls = big number if nearest body is NOT in current block
c   ***NOTE: This is a crude approximation which should be ok as long
c            as block boundaries are not VERY near bodies.  There is
c            really no EASY way to get correct information from across
c            blocks for ypls
c           nblb=int(blnum(j,k,i)+0.1)
c           if (nblb .ne. nbl) then
c             ypls=1000.
c           else
c             jbb=int(xjb(j,k,i)+0.1)
c             kbb=int(xkb(j,k,i)+0.1)
c             ibb=int(tursav(j,k,i,2)+0.1)
c             if (jbb .eq. jdim) jbb=jbb-1
c             if (kbb .eq. kdim) kbb=kbb-1
c             if (ibb .eq. idim) ibb=ibb-1
c             ibb=max(ibb,1)
c   The B-B type finding of nearest wall appears to be screwed up big time
c   So, assume k-walls only, and simply get location to nearest wall
c   location at same j and i (note- wall could be at k=1 or kdim)
              jbb=j
              ibb=i
              dist1=(x(j,k,i)-x(j,1,i))**2 + 
     +              (y(j,k,i)-y(j,1,i))**2 +
     +              (z(j,k,i)-z(j,1,i))**2
              dist2=(x(j,k,i)-x(j,kdim,i))**2 + 
     +              (y(j,k,i)-y(j,kdim,i))**2 +
     +              (z(j,k,i)-z(j,kdim,i))**2
              if (dist1 .le. dist2) then
                kbb=1
              else
                kbb=kdim-1
              end if
              tt=gamma*q(jbb,kbb,ibb,5)/q(jbb,kbb,ibb,1)
              wnu=c2bp*tt*sqrt(tt)/(c2b+tt)/q(jbb,kbb,ibb,1)
              utau=sqrt(wnu*vor(jbb,kbb,ibb)/(q(jbb,kbb,ibb,1)*re))
              ypls=re*q(jbb,kbb,ibb,1)*utau*ccabs(smin(j,k,i))/wnu
c           end if
            deltaj = 2.*vol(j,k,i)/(sj(j,k,i,4)+sj(j+1,k,i,4))
            deltak = 2.*vol(j,k,i)/(sk(j,k,i,4)+sk(j,k+1,i,4))
            deltai = 2.*vol(j,k,i)/(si(j,k,i,4)+si(j,k,i+1,4))
            if( i2d .ne. 1 ) then
               delta = (deltai*deltaj*deltak)**0.333333
            else
               delta = sqrt(deltaj*deltak)
            end if
c           Scale delta by van Driest function
            if (iwallscale .eq. 1) then
              vandriest=1.0-exp(-ypls/25.)
              delta=delta*vandriest
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
c     if(icyc.eq.ncyc1(1) .or. icyc.eq.ncyc1(2) .or. icyc.eq.ncyc1(3)
c    +.or. icyc.eq.ncyc1(4) .or. icyc.eq.ncyc1(5)) then
c       write(51,*) idim-1,jdim-1,kdim-1
c       write(51,*) (((x(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +              (((y(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +              (((z(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1)
c       write(52,*) idim-1,jdim-1,kdim-1
c       time=0.
c       write(52,*) xmach,alpha,reue,time
c       write(52,*) (((xjb(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +              (((xkb(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +         (((tursav(j,k,i,2),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +         (((tursav(j,k,i,1),i=1,idim-1),j=1,jdim-1),k=1,kdim-1),
c    +         (((vist3d(j,k,i),i=1,idim-1),j=1,jdim-1),k=1,kdim-1)
c     end if
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
