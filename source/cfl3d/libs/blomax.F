c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
      subroutine blomax(jdim,kdim,idim,q,qi0,qj0,qk0,vor,snj0,snk0,sni0,
     .                  snjm,snkm,snim,
     .                  vist3d,eoms,iprint,inmx,eomui,fbl,rhon,amun,
     .                  vortn,disn,utot,eomuo,damp,nblt,x,y,z,blank,
     .                  iover,bci,bcj,bck,nou,bou,nbuf,ibufdim,
     .                  myid,mblk2nd,maxbl)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Compute turbulent viscosity distributions using
c     Baldwin-Lomax two-layer eddy-viscosity model.
c***********************************************************************
c 
c     coding history  :  V.N.Vatsa of Langley (June-1984) 
c                        J.L. Thomas          (Dec -1987)
c 
c      input arrays   : 
c       q             : primitive variables at cell centers
c                       (rho,u,v,w,p)
c       vor           : vorticity magnitude at cell centers
c       snj0,snjm     : directed distance from j=0,max wall 
c                        of cell-centers
c       snk0,snkm     : directed distance from k=0,max wall 
c                        of cell-centers
c       sni0,snim     : directed distance from i=0,max wall 
c                       of cell-centers
c      output arrays  :
c       vist3d        : total eddy viscosity at cell centers
c      scratch arrays :
c       eoms          : sidewall eddy viscosity at cell centers
c
c**********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,5), qj0(kdim,idim-1,5,4),
     .          qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
      dimension vist3d(jdim,kdim,idim),eoms(jdim-1,kdim-1,idim-1)
      dimension vor(jdim-1,kdim-1,idim-1),blank(jdim,kdim,idim)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension snj0(jdim-1,kdim-1,idim-1),snk0(jdim-1,kdim-1,idim-1),
     .          sni0(jdim-1,kdim-1,idim-1)
      dimension snjm(jdim-1,kdim-1,idim-1),snkm(jdim-1,kdim-1,idim-1),
     .          snim(jdim-1,kdim-1,idim-1)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
      dimension mblk2nd(maxbl)
c
c      one-dimensional scratch arrays
c
      dimension eomui(inmx),fbl(inmx),rhon(inmx),amun(inmx),vortn(inmx),
     .           disn(inmx),utot(inmx),eomuo(inmx),damp(inmx)
c 
      data aplus/26.e0/,ccp/1.6e0/,ckleb/.3e0/,cwk/1.0e0/,vk/.4e0/, 
     .          clauser/.0168e0/
c
      common /degshf/ ideg(3)
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
      common /fluid2/ pr,prt,cbar
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /wallfun/ iwf(3)
      common /mgrd/ levt,kode,mode,ncyc,mtt,icyc,level,lglobal
      common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi,
     .        i_lam_forcezero
      common /reyue/ reue,tinf,ivisc(3)
      common /sklton/ isklton
      common /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze,
     .                  iwarneddy,itime2read,itaturb,tur1cut,tur2cut,
     .                  iturbord,tur1cutlev,tur2cutlev
c
      if (icyc .le. nfreeze) then
        if (isklton .gt. 0) then
           nss=min(ncyc,nfreeze)
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),*)
           nou(1) = min(nou(1)+1,ibufdim)
           write(bou(nou(1),1),'('' B-L turbulence model is frozen '',
     +     ''for '',i5,'' iterations or subits'')') nss
        end if
        return
      end if
c
c      iprint             : printing flag
c                           =0 none
c                           =1 print length scale parameters
c                           =2 print above and profile data
c
c      ilfreq1, ilfreq2   : frequency at which length scales are output 
c                           in each of the two tangential directions
c
c      ipfreq1, ipfreq2   : frequency at which profile data are output 
c                           in each of the two tangential directions
c
c                           where 1,2 denote: j,i  for k=0 viscous layer
c                          
c                                             k,i  for j=0 viscous layer
c             
c                                             j,k  for i=0 viscous layer
c
c                           e.q. for k=0 viscous layer, ilfreq1=1,ilfreq2=5,
c                           ipfreq1=10 ipfreq2=20 will print length scales 
c                           for all j stations at i=5,10,15,20..., and velocity
c                           profiles for every tenth j station at i=20,40,60...
c
       ilfreq1 = 1
       ilfreq2 = 5
       ipfreq1 = 10
       ipfreq2 = 10
       if (idim.eq.2 .or. jdim.eq.2 .or. kdim.eq.2) then
          ilfreq2 = 1
          ipfreq2 = 1
       end if
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
      aplusi = 1.e0/aplus 
      xmi    = 1.e0/xmach 
      ckout  = reue*xmi*clauser*ccp
      c2b    = cbar/tinf
      c2bp   = c2b+1.e0
c
c     iwall = 0 for yplus damping scaled with local values
c     iwall = anything else for yplus damping scaled with wall values
      iwall  = 0
c
c     ifit  = 1 for curve fit of f near infmax to improve fmax location
c     ifit  = 0 take fmax as largest descrete value of f(in) (no curve fit)
      ifit  = 0
c 
c   determine which walls to use, based on bcj,bck,bci
      ibcjmin=0
      ibcjmax=0
      do i=1,idim1
      do k=1,kdim
        if (real(bcj(k,i,1)) .gt. 0.) ibcjmin=1
        if (real(bcj(k,i,2)) .gt. 0.) ibcjmax=1
      enddo
      enddo
      ibckmin=0
      ibckmax=0
      do i=1,idim1
      do j=1,jdim
        if (real(bck(j,i,1)) .gt. 0.) ibckmin=1
        if (real(bck(j,i,2)) .gt. 0.) ibckmax=1
      enddo
      enddo
      ibcimin=0
      ibcimax=0
      do k=1,kdim
      do j=1,jdim
        if (real(bci(j,k,1)) .gt. 0.) ibcimin=1
        if (real(bci(j,k,2)) .gt. 0.) ibcimax=1
      enddo
      enddo
      if (ivisc(3) .gt. 1) then
c*************************************************************** 
c      evaluate turbulent viscosity along normal to k=0 wall 
c***************************************************************  
c 
c   defaults to using min face for distance if neither min nor max contains wall
      if (ibckmin.eq.1 .or. ibckmax.eq.0) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' computing turb viscosity using '',
     . ''Baldwin-Lomax, K=1   , block='',i4)') nblt
      end if
      if (isklton.gt.0 .and. ideg(3).eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1525)
      end if
 1525 format(33h   Degani-Schiff option turned on)
      if (isklton.gt.0 .and. ideg(3).ne.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1526)
      end if
 1526 format(34h   Degani-Schiff option turned off)
c
c      loop through i stations
c 
      if (inmx.lt.kdim) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),9)inmx,kdim
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
    9 format(19h error: inmx, kdim=,2i5)
c
      ihead = 0
      do 1000 i=1,idim1
c
c      loop through j stations
c 
      do 1000 j=1,jdim1
c 
c     create working arrays in K-direction
c 
c      bound turbulent region
c
      kloop  = .80*kdim
      kloop1 = kloop-1
      inmax  = kloop
      inmax1 = inmax-1
c
c      bound search region for fmax
c
      in1    = inmax1*0.20+1
      if(iwf(3) .gt. 0) then
        in1    = 3
      end if
      in2    = inmax1*0.80+1
c
c     Chimera scheme modification....limit upper bound of search region to
c     the lower edge of a hole.  If hole extends to lower search bound, set
c     eddy viscosity to zero at all points this station.  
c
      if (iover.eq.1) then
         if (real(blank(j,in1,i)) .eq. 0.) then
            do 1519 k=1,kdim1
            vist3d(j,k,i) = 0.e0
 1519       continue
            if (iprint.ge.1) then
            if (j.eq.j/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  nou(3) = min(nou(3)+1,ibufdim)
                  write(bou(nou(3),3),1409)
               end if   
               nou(3) = min(nou(3)+1,ibufdim)
               write(bou(nou(3),3),1408) j,i
 1408          format(i9,i5,5x,34hskipping eddy viscosity evaluation,
     .                41h at this station due to hole at-near wall)
               ihead = ihead+1
    	    end if
            end if
  	    go to 1000
         end if
         in21 = in2
         do 1523 in=in1,in2
         if (real(blank(j,in,i)) .eq. 0.) go to 1524
         in21 = in
 1523    continue
 1524    continue
         in2 = in21    
      end if 
c
c      viscosity thru sutherland relation at cell centers 
c
cdir$ ivdep
      do 9108 k=1,kloop1
      t5         =  gamma*q(j,k,i,5)/q(j,k,i,1)
      t6         =  sqrt(t5) 
      amun (k+1) =  c2bp*t5*t6/(c2b+t5)
 9108 continue
cdir$ ivdep
      do 9109 k=1,kloop1
      rhon (k+1) =  q(j,k,i,1)
      vortn(k+1) =  vor(j,k,i)
      disn (k+1) =  ccabs( snk0(j,k,i) )
      utot (k+1) =  sqrt( q(j,k,i,2)*q(j,k,i,2)
     .            +       q(j,k,i,3)*q(j,k,i,3)
     .            +       q(j,k,i,4)*q(j,k,i,4) )
 9109 continue
c
      rhon(1)   =  rhon(2)
      amun(1)   =  amun(2)
      disn(1)   =  0.
      utot(1)   =  utot(2)
      vortn(1)  =  vortn(2)
c
      if (real(bck(j,i,1)) .eq. 1.) then
c
c     modify wall values in no-slip region
c
      t5 = gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
      t6 = sqrt(t5)
c
      amun(1)  = c2bp*t5*t6/(c2b+t5)
      rhon(1)  = qk0(j,i,1,1)
      utot(1)  = sqrt(qk0(j,i,2,1)*qk0(j,i,2,1)+
     .               qk0(j,i,3,1)*qk0(j,i,3,1)+
     .               qk0(j,i,4,1)*qk0(j,i,4,1))
c
      vortn(1) = (utot(2)-utot(1))/disn(2)
      end if
c 
c     convenient groupings and initial values at k=1
c 
      yplusoy = ccabs(rhon(1)*vortn(1)*reue*xmi/amun(1) ) 
      yplusoy = sqrt(yplusoy)*aplusi 
      yplus   = yplusoy*disn(2)*aplus*2.0
c 
      if (real(bck(j,i,1)) .eq. 0.) then
c        wake damping
         term     = 1.-exp(-50.e0)
         do 9111 in=2,inmax
         damp(in) = term
 9111    continue
      else
c        wall damping 
         do 9112 in=2,inmax
         if (iwall .eq. 0) then
            term = yplusoy*disn(in)*sqrt(rhon(in)/rhon(1))*
     .      (amun(1)/amun(in))
         else
            term = yplusoy*disn(in)
         end if
         damp(in) = 1.0e0 -exp(-term) 
 9112    continue
      end if
c
      fbl(1)    = 0.0e0
      eomui(1)  = 0.0e0
      do 9113 in=2,inmax
      fbl(in)   = vortn(in)*disn(in)*damp(in)
      amixl     = vk*disn(in)*damp(in)
      eomui(in) = reue*xmi*rhon(in)*vortn(in)* 
     .            amixl*amixl
 9113 continue
c
      utmax = q8smax(inmax,utot(1))
      utmin = q8smin(inmax,utot(1))
      if (real(bck(j,i,1)) .eq. 1.) utmin = utot(1)
c 
c     locate max value and location of fbl
c
      fblmax = 1.e-10
      infmax = in1
      do 9115 in=in1,in2
      if (real(fbl(in)).gt.real(fblmax)) then
         fblmax = fbl(in)
         infmax = in
      end if
c
c      first maximum (degani-schiff)
c
      if (ideg(3).eq.1.and.real(fbl(in)).lt.0.9e0*real(fblmax))goto 9125
 9115 continue
 9125 continue
c
c     curve fit of f near infmax to improve fmax location
c
      if (ifit.gt.0) then
         if (infmax.gt.in1 .and. infmax.lt.in2) then
            dfm = fbl(infmax) - fbl(infmax-1)
            dfp = fbl(infmax) - fbl(infmax+1)
            if (real(fbl(infmax-1)).lt.real(fbl(infmax+1))) then
               if (real(dfm).gt.0.) then
                  ymax = disn(infmax) + 0.5*(disn(infmax+1) - 
     .                                 disn(infmax))*(1.-dfp/dfm)
               else
                  ymax = disn(infmax)
               end  if
            else
               if (real(dfp).gt.0.) then
                  ymax = disn(infmax) - 0.5*(disn(infmax) - 
     .                                 disn(infmax-1))*(1.-dfm/dfp)
               else
                  ymax = disn(infmax)
               end if
            end if
         end if
      else
         ymax   = disn(infmax) 
      end if
c
      ymax   = ccmaxcr(ymax,1.e-20)
      udif   = utmax - utmin 
      fwake  = ymax*fblmax 
      fwake2 = cwk*udif*udif*ymax/fblmax 
      fwake  = ccmin(fwake2,fwake)
c
      do 9116 in=2,inmax 
      eomuo(in) = 1.e0/(1.e0+5.5e0*(ckleb*disn(in)/ymax)**6)
      eomuo(in) = ckout*eomuo(in)*rhon(in)*fwake
 9116 continue
      eomuo(1)  = eomuo(2)
c
      do 9117 in=2,in2
      if (real(eomui(in)).gt.real(eomuo(in))) go to 9127
 9117 continue
 9127 inswtch   = in-1
      do 9118 in=inswtch,inmax
      eomui(in) = eomuo(in)
 9118 continue
c
c      fill eomu array corresponding to k=0 wall 
c
      do 900 k=1,kloop1
      vist3d(j,k,i) = eomui(k+1)
  900 continue
      if (kloop1.lt.kdim1) then
         do 1907 k=kloop,kdim1
         vist3d(j,k,i) = 0.
 1907    continue
      end if
c
      if (iprint.ge.1) then
      if (j.eq.(j/ilfreq1)*ilfreq1 .and. i.eq.(i/ilfreq2)*ilfreq2) then
c
c      length scale parameters
c
      if (ihead.eq.0) then
         nou(3) = min(nou(3)+1,ibufdim)
         write(bou(nou(3),3),1409)
      end if
 1409 format(1x,3hk=0,4x,1hj,4x,1hi,5x,5hyplus,5x,5hutmax,5x,5hut(2),
     .6x,4hfmax,6x,4hymax,4x,6hemumax,1x,5hswtch,
     .2x,4hfmax,1x,3hin1,1x,3hin2,1x,3hmax,1x,3h1/2)
      ihead = ihead+1
      eomax = q8smax(inmax1,eomui(2))
      utmax = utmax/xmach
      utmin = utmin/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1410) j,i,real(yplus),real(utmax),
     .      real(utmin),real(fblmax),real(ymax),real(eomax),
     .      inswtch,infmax,in1,in2,inmax
 1410 format(i9,i5,3e10.3,3e10.3,i6,i6,3i4)
      end if
      end if
c
c      profile data
c
      if (iprint.gt.1) then
      if (j.eq.(j/ipfreq1)*ipfreq1 .and. i.eq.(i/ipfreq2)*ipfreq2) then
      ihead = 0
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1413)
 1413 format(4x,1hj,4x,1hi,4x,1hk,6x,4hdisn,6x,4hutot,6x,4heomu,
     .7x,3hfbl,6x,4hamun,4x,6heomu-o,6x,4hvort)
      do 1415 in=1,inmax
      utw = utot(in)/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1414) j,i,in,real(disn(in)),real(utw),
     .               real(eomui(in)),real(fbl(in)),
     .               real(amun(in)),real(eomuo(in)),real(vortn(in))
 1414 format(3i5,7e10.3)
 1415 continue
      end if
      end if
 1000 continue
c   end if on k=0 face:
      end if
c***************************************************************
c      evaluate turbulent viscosity along normal to k=kdim wall
c***************************************************************
c
      if (ibckmax.eq.1) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' computing turb viscosity using '',
     . ''Baldwin-Lomax, K=kdim, block='',i4)') nblt
      end if
      if (isklton.gt.0 .and. ideg(3).eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1525)
      end if
      if (isklton.gt.0 .and. ideg(3).ne.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1526)
      end if
c
c      loop through i stations
c
      if (inmx.lt.kdim) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),9)inmx,kdim
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      ihead = 0
      do 10 i=1,idim1
c
c      loop through j stations
c
      do 10 j=1,jdim1
c
c     create working arrays in K-direction
c
c      bound turbulent region
c
      kloop  = .80*kdim
      kloop1 = kloop-1
      inmaxt = kloop
      inmaxt1 = inmaxt-1
c
      kloop  = kdim1 - kloop + 1
      kloop  = max(1,kloop)
      kloop1 = kloop + 1
      inmax  = kloop
      inmax1 = inmax + 1
c
c      bound search region for fmax
c
      in1    = inmaxt1*0.20+1
      if(iwf(3) .gt. 0) then
        in1    = 3
      end if
      in2    = inmaxt1*0.80+1
      in1    = kdim1 - in1 + 1
      in2    = kdim1 - in2 + 1
c
c     Chimera scheme modification....limit upper bound of search region to
c     the lower edge of a hole.  If hole extends to lower search bound, set
c     eddy viscosity to zero at all points this station.
c
      if (iover.eq.1) then
         if (real(blank(j,in1,i)) .eq. 0.) then
            do 15 k=1,kdim1
            vist3d(j,k,i) = 0.e0
   15       continue
            if (iprint.ge.1) then
            if (j.eq.j/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  nou(3) = min(nou(3)+1,ibufdim)
                  write(bou(nou(3),3),1409)
               end if
               nou(3) = min(nou(3)+1,ibufdim)
               write(bou(nou(3),3),1408) j,i
               ihead = ihead+1
            end if
            end if
            go to 10
         end if
         in21 = in2
         do 23 in=in1,in2,-1
         if (real(blank(j,in,i)) .eq. 0.) go to 24
         in21 = in
   23    continue
   24    continue
         in2 = in21
      end if
c
c      viscosity thru sutherland relation at cell centers
c
cdir$ ivdep
      do 3478 k=kdim1,kloop1,-1
      t5         =  gamma*q(j,k,i,5)/q(j,k,i,1)
      t6         =  sqrt(t5)
      amun (k) =  c2bp*t5*t6/(c2b+t5)
 3478 continue
cdir$ ivdep
c     ierr = 0
      do 29 k=kdim1,kloop1,-1
      rhon (k) =  q(j,k,i,1)
      vortn(k) =  vor(j,k,i)
      disn (k) =  ccabs( snkm(j,k,i) )
      utot (k) =  sqrt(   q(j,k,i,2)*q(j,k,i,2)
     .              +       q(j,k,i,3)*q(j,k,i,3)
     .              +       q(j,k,i,4)*q(j,k,i,4) )
   29 continue
c
      rhon(kdim)   =  rhon(kdim1)
      amun(kdim)   =  amun(kdim1)
      disn(kdim)   =  0.
      utot(kdim)   =  utot(kdim1)
      vortn(kdim)  =  vortn(kdim1)
c
      if (real(bck(j,i,2)) .eq. 1.) then
c
c     modify wall values in no-slip region
c
      t5 = gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
      t6 = sqrt(t5)
c
      amun(kdim)  = c2bp*t5*t6/(c2b+t5)
      rhon(kdim)  = qk0(j,i,1,3)
      utot(kdim)  = sqrt(qk0(j,i,2,3)*qk0(j,i,2,3)+
     .                   qk0(j,i,3,3)*qk0(j,i,3,3)+
     .                   qk0(j,i,4,3)*qk0(j,i,4,3))
c
      vortn(kdim) = (utot(kdim1)-utot(kdim))/disn(kdim1)
      end if
c
c     convenient groupings and initial values at k=kdim
c
      yplusoy = ccabs(rhon(kdim)*vortn(kdim)*reue*xmi/amun(kdim) )
      yplusoy = sqrt(yplusoy)*aplusi
      yplus   = yplusoy*disn(kdim1)*aplus*2.0
c
      if (real(bck(j,i,2)) .eq. 0.) then
c        wake damping
         term     = 1.-exp(-50.e0)
         do 11 in=kdim1,inmax1,-1
         damp(in) = term
   11    continue
      else
c        wall damping
         do   12 in=kdim1,inmax1,-1
         if (iwall .eq. 0) then
           term = yplusoy*disn(in)*sqrt(rhon(in)/rhon(kdim))*
     .      (amun(kdim)/amun(in))
         else
            term = yplusoy*disn(in)
         end if
         damp(in) = 1.0e0 -exp(-term)
   12    continue
      end if
c
      fbl(kdim)    = 0.0e0
      eomui(kdim)  = 0.0e0
      do   13 in=kdim1,inmax1,-1
      fbl(in)   = vortn(in)*disn(in)*damp(in)
      amixl     = vk*disn(in)*damp(in)
      eomui(in) = reue*xmi*rhon(in)*vortn(in)*
     .            amixl*amixl
   13 continue
c
      utmax = q8smax(inmaxt1,utot(inmax1))
      utmin = q8smin(inmaxt1,utot(inmax1))
      if (real(bck(j,i,2)) .eq. 1.) utmin = utot(kdim)
c
c     locate max value and location of fbl
c
      fblmax = 1.e-10
      infmax = in1+1
      do 115 in=in1+1,in2+1,-1
      if (real(fbl(in)).gt.real(fblmax)) then
         fblmax = fbl(in)
         infmax = in
      end if
c
c      first maximum (degani-schiff)
c
      if (ideg(3).eq.1.and.real(fbl(in)).lt.0.9e0*real(fblmax))goto 125
  115 continue
  125 continue
c
c     curve fit of f near infmax to improve fmax location
c
      if (ifit.gt.0) then
         if (infmax.lt.in1 .and. infmax.gt.in2) then
            dfm = fbl(infmax) - fbl(infmax+1)
            dfp = fbl(infmax) - fbl(infmax-1)
            if (real(fbl(infmax+1)).lt.real(fbl(infmax-1))) then
               if (real(dfm).gt.0.) then
                  ymax = disn(infmax) + 0.5*(disn(infmax-1) -
     .                                 disn(infmax))*(1.-dfp/dfm)
               else
                  ymax = disn(infmax)
               end  if
            else
               if (real(dfp).gt.0.) then
                  ymax = disn(infmax) - 0.5*(disn(infmax) -
     .                                 disn(infmax+1))*(1.-dfm/dfp)
               else
                  ymax = disn(infmax)
               end if
            end if
         end if
      else
         ymax   = disn(infmax)
      end if
c
      ymax   = ccmaxcr(ymax,1.e-20)
      udif   = utmax - utmin
      fwake  = ymax*fblmax
      fwake2 = cwk*udif*udif*ymax/fblmax
      fwake  = ccmin(fwake2,fwake)
c
      do 116 in=kdim1,inmax1,-1
      eomuo(in) = 1.e0/(1.e0+5.5e0*(ckleb*disn(in)/ymax)**6)
      eomuo(in) = ckout*eomuo(in)*rhon(in)*fwake
  116 continue
      eomuo(kdim)  = eomuo(kdim1)
c
      do 117 in=kdim1,in2,-1
      if (real(eomui(in)).gt.real(eomuo(in))) go to 127
  117 continue
  127 inswtch   = in+1
      do 118 in=inswtch,inmax1,-1
      eomui(in) = eomuo(in)
  118 continue
c
c      fill eomu array corresponding to k=kdim wall
c
      if (ibckmin.eq.1) then
c
c     for k=1 & k=kdim walls
      do 90 k=kdim1,kloop1,-1
      s1 = snk0(j,k,i)*snk0(j,k,i)
      s2 = snkm(j,k,i)*snkm(j,k,i)
      vist3d(j,k,i) = ( s1*eomui(k) + s2*vist3d(j,k,i) )/( s1 + s2 )
   90 continue
c
      else
c  
c     for k=kdim wall only
      do k=kdim1,kloop1,-1
         vist3d(j,k,i) = eomui(k)
      enddo
      if (kloop1 .gt. 1) then
         do k = kloop,1,-1
            vist3d(j,k,i) = 0.
         enddo
      end if
c
      end if
c
c
      if (iprint.ge.1) then
      if (j.eq.(j/ilfreq1)*ilfreq1 .and. i.eq.(i/ilfreq2)*ilfreq2) then
c
c      length scale parameters
c
      if (ihead.eq.0) then
         nou(3) = min(nou(3)+1,ibufdim)
         write(bou(nou(3),3),1409)
      end if
      ihead = ihead+1
      eomax = q8smax(inmax1,eomui(kdim1))
      utmax = utmax/xmach
      utmin = utmin/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1410) j,i,real(yplus),real(utmax),
     .      real(utmin),real(fblmax),real(ymax),real(eomax),
     .      inswtch,infmax,in1,in2,inmax
      end if
      end if
c
c      profile data
c
      if (iprint.gt.1) then
      if (j.eq.(j/ipfreq1)*ipfreq1 .and. i.eq.(i/ipfreq2)*ipfreq2) then
      ihead = 0
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1413)
      do 95 in=kdim1,inmax1,-1
      utw = utot(in)/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1414) j,i,in,real(disn(in)),real(utw),
     .               real(eomui(in)),real(fbl(in)),
     .               real(amun(in)),real(eomuo(in)),real(vortn(in))
   95 continue
      end if
      end if
   10 continue
c
c   end if on k=kdim face:
      end if
c
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-L turb model'',
     .   '' has no laminar regions'')') nblt
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-L turb model -'',
     .   '' laminar region is:'')') nblt
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   i='',i5,'' to'',i5,'', j='',i5,
     .   '' to'',i5,'', k='',i5,'' to'',i5)') ilamlo,ilamhi,jlamlo,
     .   jlamhi,klamlo,klamhi
        end if
      end if
      do 61 i=ilamlo,ilamhi-1
      do 61 k=klamlo,klamhi-1
      do 61 j=jlamlo,jlamhi-1
        vist3d(j,k,i)=0.
   61 continue
c  end if on (ivisc(3) .gt. 1):
      end if
c
      if (ivisc(2) .gt. 1) then
c*************************************************************** 
c      evaluate turbulent viscosity along normal to j=0 wall  
c*************************************************************** 
c 
c   defaults to using min face for distance if neither min nor max contains wall
      if (ibcjmin.eq.1 .or. ibcjmax.eq.0) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' computing turb viscosity using '',
     . ''Baldwin-Lomax, J=1   , block='',i4)') nblt
      end if
      if (isklton.gt.0 .and. ideg(2).eq.1) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),1525)
      end if
      if (isklton.gt.0 .and. ideg(2).ne.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1526)
      end if
c 
c      loop through i stations
c 
      if (inmx.lt.jdim) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),39)inmx,jdim
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
   39 format(19h error: inmx, jdim=,2i5)
c
      ihead = 0
      do 7000 i=1,idim1
c
c      loop through k stations
c 
      do 7000 k=1,kdim1
c 
c     create working arrays in J-direction (normal to j=0 wall)
c 
c      bound turbulent region
c
      jloop  = .80*jdim
      jloop1 = jloop-1
      inmax  = jloop
      inmax1 = inmax-1
c
c      bound search region for fmax
c
      in1    = inmax1*0.20+1
      if(iwf(2) .gt. 0) then
        in1    = 3
      end if
      in2    = inmax1*0.80+1
c
c     Chimera scheme modification....limit upper bound of search region to
c     the lower edge of a hole.  If hole extends to lower search bound, set
c     eddy viscosity to zero at all points this station.  
c
      if (iover.eq.1) then
         if (real(blank(in1,k,i)) .eq. 0.) then
            do 1619 j=1,jdim1
            eoms(j,k,i) = 0.e0
 1619       continue                                                  
            if (iprint.ge.1) then
            if (k.eq.k/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  nou(3) = min(nou(3)+1,ibufdim)
                  write(bou(nou(3),3),1409)                           
               end if
               nou(3) = min(nou(3)+1,ibufdim)
               write(bou(nou(3),3),1408) k,i
               ihead = ihead+1
	    end if
            end if
	    go to 7000
         end if
         in21 = in2
         do 1623 in=in1,in2
         if (real(blank(in,k,i)) .eq. 0.) go to 1624
         in21 = in
 1623    continue
 1624    continue    
         in2 = in21
      end if 
c
c      viscosity thru sutherland relation at cell centers 
c
cdir$ ivdep
      do 6108 j=1,jloop1
      t5         = gamma*q(j,k,i,5)/q(j,k,i,1)
      t6         = sqrt(t5) 
      amun (j+1) = c2bp*t5*t6/(c2b+t5)
 6108 continue
cdir$ ivdep
      do 6109 j=1,jloop1
      rhon (j+1) = q(j,k,i,1)
      vortn(j+1) = vor(j,k,i)
      disn (j+1) = ccabs( snj0(j,k,i) )
      utot (j+1) = sqrt( q(j,k,i,2)*q(j,k,i,2)
     .            +      q(j,k,i,3)*q(j,k,i,3)
     .            +      q(j,k,i,4)*q(j,k,i,4) )
 6109 continue
c
      rhon(1)  =  rhon(2)
      amun(1)  =  amun(2)
      disn(1)  =  0.
      utot(1)  =  utot(2)
      vortn(1) =  vortn(2)
c
      if (real(bcj(k,i,1)) .eq. 1.) then
c
c     modify wall values in no-slip region
c
      t5 = gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
      t6 = sqrt(t5)
c
      amun(1) = c2bp*t5*t6/(c2b+t5)
      rhon(1) = qj0(k,i,1,1)
      utot(1) = sqrt(qj0(k,i,2,1)*qj0(k,i,2,1)+
     .               qj0(k,i,3,1)*qj0(k,i,3,1)+
     .               qj0(k,i,4,1)*qj0(k,i,4,1))
c
      vortn(1) = (utot(2)-utot(1))/disn(2)
      end if
c 
c     convenient groupings and initial values at j=1
c 
      yplusoy = ccabs(rhon(1)*vortn(1)*reue*xmi/amun(1) ) 
      yplusoy = sqrt(yplusoy)*aplusi 
      yplus   = yplusoy*disn(2)*aplus*2.0
c
      if (real(bcj(k,i,1)) .eq. 0.) then
c        wake damping
         term     = 1.-exp(-50.e0)
         do 6111 in=2,inmax
         damp(in) = term
 6111    continue
      else
c        wall damping 
         do 6112 in=2,inmax
         if (iwall .eq. 0) then
            term = yplusoy*disn(in)*sqrt(rhon(in)/rhon(1))*
     .      (amun(1)/amun(in))
         else
            term = yplusoy*disn(in)
         end if
         damp(in) = 1.0e0 -exp(-term) 
 6112    continue
      end if
c
      fbl(1)    = 0.0e0
      eomui(1)  = 0.0e0
      do 6113 in=2,inmax
      fbl(in)   = vortn(in)*disn(in)*damp(in)
      amixl     = vk*disn(in)*damp(in)
      eomui(in) = reue*xmi*rhon(in)*vortn(in)* 
     .            amixl*amixl
 6113 continue
c
      utmax = q8smax(inmax,utot(1))
      utmin = q8smin(inmax,utot(1))
      if (real(bcj(k,i,1)) .eq. 1.) utmin = utot(1)
c 
c     locate max value and location of fbl
c
      fblmax = 1.e-10
      infmax = in1
      do 6115 in=in1,in2
      if (real(fbl(in)).gt.real(fblmax)) then
         fblmax = fbl(in)
         infmax = in
      end if
c
c      first maximum (degani-schiff)
c
      if (ideg(2).eq.1.and.real(fbl(in)).lt.0.9e0*real(fblmax))goto 6125
 6115 continue
 6125 continue
c
c     curve fit of f near infmax to improve fmax location
c
      if (ifit.gt.0) then
         if (infmax.gt.in1 .and. infmax.lt.in2) then
            dfm = fbl(infmax) - fbl(infmax-1)
            dfp = fbl(infmax) - fbl(infmax+1)
            if (real(fbl(infmax-1)).lt.real(fbl(infmax+1))) then
               if (real(dfm).gt.0.) then
                  ymax = disn(infmax) + 0.5*(disn(infmax+1) -
     .                                 disn(infmax))*(1.-dfp/dfm)
               else
                  ymax = disn(infmax)
               end  if
            else
               if (real(dfp).gt.0.) then
                  ymax = disn(infmax) - 0.5*(disn(infmax) -
     .                                 disn(infmax-1))*(1.-dfm/dfp)
               else
                  ymax = disn(infmax)
               end if
            end if
         end if
      else
         ymax   = disn(infmax)
      end if
c
      ymax  = ccmaxcr(ymax,1.e-20)
      udif  = utmax - utmin 
      fwake = ymax*fblmax 
      fwake2= cwk*udif*udif*ymax/fblmax 
      fwake = ccmin(fwake2,fwake)
c
      do 6116 in=2,inmax 
      eomuo(in) = 1.e0/(1.e0+5.5e0*(ckleb*disn(in)/ymax)**6)
      eomuo(in) = ckout*eomuo(in)*rhon(in)*fwake
 6116 continue
      eomuo(1)  =  eomuo(2)
c
      do 6117 in=2,in2
      if (real(eomui(in)).gt.real(eomuo(in))) go to 6127
 6117 continue
 6127 inswtch   = in-1
      do 6118 in=inswtch,inmax
      eomui(in) = eomuo(in)
 6118 continue
c
c      fill eomu array corresponding to j=0 wall 
c
      do 600 j=1,jloop1
      eoms(j,k,i) = eomui(j+1)
  600 continue
      if (jloop1.lt.jdim1) then
         do 1607 j=jloop,jdim1
         eoms(j,k,i) = 0.
 1607    continue
      end if
c 
      if (iprint.ge.1) then
      if (k.eq.(k/ilfreq1)*ilfreq1 .and. i.eq.(i/ilfreq2)*ilfreq2) then
c
c      length scale parameters
c
      if (ihead.eq.0) then
         nou(3) = min(nou(3)+1,ibufdim)
         write(bou(nou(3),3),1406)
      end if
 1406 format(1x,3hj=0,4x,1hk,4x,1hi,5x,5hyplus,5x,5hutmax,5x,5hut(2),
     .6x,4hfmax,6x,4hymax,4x,6hemumax,1x,5hswtch,
     .2x,4hfmax,1x,3hin1,1x,3hin2,1x,3hmax,1x,3h1/2)
      ihead = ihead+1
      eomax = q8smax(inmax1,eomui(2))
      utmax = utmax/xmach
      utmin = utmin/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1610) k,i,real(yplus),real(utmax),
     .      real(utmin),real(fblmax),real(ymax),real(eomax),
     .      inswtch,infmax,in1,in2,inmax
 1610 format(i9,i5,3e10.3,3e10.3,i6,i6,4i4)
      end if
      end if
c
c      profile data
c
      if (iprint.gt.1) then
      if (k.eq.(k/ipfreq1)*ipfreq1 .and. i.eq.(i/ipfreq2)*ipfreq2) then
      ihead = 0
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1613)
 1613 format(4x,1hk,4x,1hi,4x,1hj,6x,4hdisn,6x,4hutot,6x,4heomu,
     .7x,3hfbl,6x,4hamun,4x,6heomu-o,6x,4hvort)
      do 1615 in=1,inmax
      utw = utot(in)/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1414) k,i,in,real(disn(in)),real(utw),
     .               real(eomui(in)),real(fbl(in)),
     .               real(amun(in)),real(eomuo(in)),real(vortn(in))
 1614 format(3i5,7e10.3)
 1615 continue
      end if
      end if
 7000 continue
c
c   end if on j=0 face:
      end if
c***************************************************************
c      evaluate turbulent viscosity along normal to j=jdim wall
c***************************************************************
c
      if (ibcjmax.eq.1) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' computing turb viscosity using '',
     . ''Baldwin-Lomax, J=jdim, block='',i4)') nblt
      end if
      if (isklton.gt.0 .and. ideg(2).eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1525)
      end if
      if (isklton.gt.0 .and. ideg(2).ne.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1526)
      end if
c
c      loop through i stations
c
      if (inmx.lt.jdim) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),39)inmx,jdim
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
c
      ihead = 0
      do 7610 i=1,idim1
c
c      loop through k stations
c
      do 7610 k=1,kdim1
c
c     create working arrays in J-direction
c
c      bound turbulent region
c
      jloop  = .80*jdim
      jloop1 = jloop-1
      inmaxt = jloop
      inmaxt1 = inmaxt-1
c
      jloop  = jdim1 - jloop + 1
      jloop  = max(1,jloop)
      jloop1 = jloop + 1
      inmax  = jloop
      inmax1 = inmax + 1
c
c      bound search region for fmax
c
      in1    = inmaxt1*0.20+1
      if(iwf(2) .gt. 0) then
        in1    = 3
      end if
      in2    = inmaxt1*0.80+1
      in1    = jdim1 - in1 + 1
      in2    = jdim1 - in2 + 1
c
c     Chimera scheme modification....limit upper bound of search region to
c     the lower edge of a hole.  If hole extends to lower search bound, set
c     eddy viscosity to zero at all points this station.
c
      if (iover.eq.1) then
         if (real(blank(in1,k,i)) .eq. 0.) then
            do 615 j=1,jdim1
            eoms(j,k,i) = 0.e0
  615       continue
            if (iprint.ge.1) then
            if (k.eq.k/ilfreq1*ilfreq1.and.i.eq.i/ilfreq2*ilfreq2) then
               if (ihead.eq.0) then
                  nou(3) = min(nou(3)+1,ibufdim)
                  write(bou(nou(3),3),1409)
               end if
               nou(3) = min(nou(3)+1,ibufdim)
               write(bou(nou(3),3),1408) k,i
               ihead = ihead+1
            end if
            end if
            go to 7610
         end if
         in21 = in2
         do 623 in=in1,in2,-1
         if (real(blank(in,k,i)) .eq. 0.) go to 624
         in21 = in
  623    continue
  624    continue
         in2 = in21
      end if
c
c      viscosity thru sutherland relation at cell centers
c
cdir$ ivdep
      do 8 j=jdim1,jloop1,-1
      t5         =  gamma*q(j,k,i,5)/q(j,k,i,1)
      t6         =  sqrt(t5)
      amun (j) =  c2bp*t5*t6/(c2b+t5)
    8 continue
cdir$ ivdep
c     ierr = 0
      do 629 j=jdim1,jloop1,-1
      rhon (j) =  q(j,k,i,1)
      vortn(j) =  vor(j,k,i)
      disn (j) =  ccabs( snjm(j,k,i) )
      utot (j) =  sqrt(   q(j,k,i,2)*q(j,k,i,2)
     .              +       q(j,k,i,3)*q(j,k,i,3)
     .              +       q(j,k,i,4)*q(j,k,i,4) )
  629 continue
c
      rhon(jdim)   =  rhon(jdim1)
      amun(jdim)   =  amun(jdim1)
      disn(jdim)   =  0.
      utot(jdim)   =  utot(jdim1)
      vortn(jdim)  =  vortn(jdim1)
c
      if (real(bcj(k,i,2)) .eq. 1.) then
c
c     modify wall values in no-slip region
c
      t5 = gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
      t6 = sqrt(t5)
c
      amun(jdim)  = c2bp*t5*t6/(c2b+t5)
      rhon(jdim)  = qj0(k,i,1,3)
      utot(jdim)  = sqrt(qj0(k,i,2,3)*qj0(k,i,2,3)+
     .                   qj0(k,i,3,3)*qj0(k,i,3,3)+
     .                   qj0(k,i,4,3)*qj0(k,i,4,3))
c
      vortn(jdim) = (utot(jdim1)-utot(jdim))/disn(jdim1)
      end if
c
c     convenient groupings and initial values at j=jdim
c
      yplusoy = ccabs(rhon(jdim)*vortn(jdim)*reue*xmi/amun(jdim) )
      yplusoy = sqrt(yplusoy)*aplusi
      yplus   = yplusoy*disn(jdim1)*aplus*2.0
c
      if (real(bcj(k,i,2)) .eq. 0.) then
c        wake damping
         term     = 1.-exp(-50.e0)
         do 511 in=jdim1,inmax1,-1
         damp(in) = term
  511    continue
      else
c        wall damping
         do  512 in=jdim1,inmax1,-1
         if (iwall .eq. 0) then
           term = yplusoy*disn(in)*sqrt(rhon(in)/rhon(jdim))*
     .      (amun(jdim)/amun(in))
         else
            term = yplusoy*disn(in)
         end if
         damp(in) = 1.0e0 -exp(-term)
  512    continue
      end if
c
      fbl(jdim)    = 0.0e0
      eomui(jdim)  = 0.0e0
      do  513 in=jdim1,inmax1,-1
      fbl(in)   = vortn(in)*disn(in)*damp(in)
      amixl     = vk*disn(in)*damp(in)
      eomui(in) = reue*xmi*rhon(in)*vortn(in)*
     .            amixl*amixl
  513 continue
c
      utmax = q8smax(inmaxt1,utot(inmax1))
      utmin = q8smin(inmaxt1,utot(inmax1))
      if (real(bcj(k,i,2)) .eq. 1.) utmin = utot(jdim)
c
c     locate max value and location of fbl
c
      fblmax = 1.e-10
      infmax = in1+1
      do 915 in=in1+1,in2+1,-1
      if (real(fbl(in)).gt.real(fblmax)) then
         fblmax = fbl(in)
         infmax = in
      end if
c
c      first maximum (degani-schiff)
c
      if (ideg(2).eq.1.and.real(fbl(in)).lt.0.9e0*real(fblmax))goto  925
  915 continue
  925 continue
c
c     curve fit of f near infmax to improve fmax location
c
      if (ifit.gt.0) then
         if (infmax.lt.in1 .and. infmax.gt.in2) then
            dfm = fbl(infmax) - fbl(infmax+1)
            dfp = fbl(infmax) - fbl(infmax-1)
            if (real(fbl(infmax+1)).lt.real(fbl(infmax-1))) then
               if (real(dfm).gt.0.) then
                  ymax = disn(infmax) + 0.5*(disn(infmax-1) -
     .                                 disn(infmax))*(1.-dfp/dfm)
               else
                  ymax = disn(infmax)
               end  if
            else
               if (real(dfp).gt.0.) then
                  ymax = disn(infmax) - 0.5*(disn(infmax) -
     .                                 disn(infmax+1))*(1.-dfm/dfp)
               else
                  ymax = disn(infmax)
               end if
            end if
         end if
      else
         ymax   = disn(infmax)
      end if
c
      ymax   = ccmaxcr(ymax,1.e-20)
      udif   = utmax - utmin
      fwake  = ymax*fblmax
      fwake2 = cwk*udif*udif*ymax/fblmax
      fwake  = ccmin(fwake2,fwake)
c
      do 186 in=jdim1,inmax1,-1
      eomuo(in) = 1.e0/(1.e0+5.5e0*(ckleb*disn(in)/ymax)**6)
      eomuo(in) = ckout*eomuo(in)*rhon(in)*fwake
  186 continue
      eomuo(jdim)  = eomuo(jdim1)
c
      do 187 in=jdim1,in2,-1
      if (real(eomui(in)).gt.real(eomuo(in))) go to 137
  187 continue
  137 inswtch   = in+1
      do 178 in=inswtch,inmax1,-1
      eomui(in) = eomuo(in)
  178 continue
c
c      fill eomu array corresponding to j=jdim wall
c
      if (ibcjmin.eq.1) then
c
c     for j=1 & j=jdim walls
      do 790 j=jdim1,jloop1,-1
      s1 = snj0(j,k,i)*snj0(j,k,i)
      s2 = snjm(j,k,i)*snjm(j,k,i)
      eoms(j,k,i) = ( s1*eomui(j) + s2*eoms(j,k,i) )/( s1 + s2 )
  790 continue
c
      else
c  
c     for j=jdim wall only
      do j=jdim1,jloop1,-1
         eoms(j,k,i) = eomui(j)
      enddo
      if (jloop1 .gt. 1) then
         do j = jloop,1,-1
            eoms(j,k,i) = 0.
         enddo
      end if
c
      end if
c
c
      if (iprint.ge.1) then
      if (k.eq.(k/ilfreq1)*ilfreq1 .and. i.eq.(i/ilfreq2)*ilfreq2) then
c
c      length scale parameters
c
      if (ihead.eq.0) then
         nou(3) = min(nou(3)+1,ibufdim)
         write(bou(nou(3),3),1409)
      end if
      ihead = ihead+1
      eomax = q8smax(inmax1,eomui(jdim1))
      utmax = utmax/xmach
      utmin = utmin/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1410) k,i,real(yplus),real(utmax),
     .      real(utmin),real(fblmax),real(ymax),real(eomax),
     .      inswtch,infmax,in1,in2,inmax
      end if
      end if
c
c      profile data
c
      if (iprint.gt.1) then
      if (k.eq.(k/ipfreq1)*ipfreq1 .and. i.eq.(i/ipfreq2)*ipfreq2) then
      ihead = 0
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1413)
      do 495 in=jdim1,inmax1,-1
      utw = utot(in)/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1414) k,i,in,real(disn(in)),real(utw),
     .               real(eomui(in)),real(fbl(in)),
     .               real(amun(in)),real(eomuo(in)),real(vortn(in))
  495 continue
      end if
      end if
 7610 continue
c
c   end if on j=jdim face:
      end if
c
      if (ivisc(3).gt.1) then
c     form the composite eddy-viscosity from k & j walls
c 
      do 150 i=1,idim1
      do 150 k=1,kdim1
      do 150 j=1,jdim1
      if(ibckmin.eq.0 .and. ibckmax.eq.1) then
        s1            = snkm(j,k,i)**2
      else if(ibckmin.eq.1 .and. ibckmax.eq.1) then
        s1            = ccmin(snk0(j,k,i),snkm(j,k,i))
        s1            = s1**2
      else
        s1            = snk0(j,k,i)**2
      end if
      if(ibcjmin.eq.0 .and. ibcjmax.eq.1) then
        s2            = snjm(j,k,i)**2
      else if(ibcjmin.eq.1 .and. ibcjmax.eq.1) then
        s2            = ccmin(snj0(j,k,i),snjm(j,k,i))
        s2            = s2**2
      else
        s2            = snj0(j,k,i)**2
      end if
      vist3d(j,k,i) = ( vist3d(j,k,i)*s2  + eoms(j,k,i)*s1 )/( s1 + s2 )
  150 continue
c
c      smooth interior values
c
      nsmurf = 3 
      if (nsmurf.eq.0) return
c
      do 185 nnn=1,nsmurf
c
      do 170 i=1,idim1
      do 170 k=1,kdim1
      do 170 j=1,jdim1
      eoms(j,k,i) = vist3d(j,k,i)
  170 continue
c
      do 176 i=2,idim1-1
      do 176 k=2,kdim1-1
      do 176 j=2,jdim1-1
      vist3d(j,k,i) = .125e0* ( eoms(j+1,k-1,i-1)+ eoms(j-1,k-1,i-1)
     .               +          eoms(j+1,k+1,i-1)+ eoms(j-1,k+1,i-1)
     .               +          eoms(j+1,k-1,i+1)+ eoms(j-1,k-1,i+1)
     .               +          eoms(j+1,k+1,i+1)+ eoms(j-1,k+1,i+1)) 
  176 continue
  185 continue
c
      else
c
c      single turbulent wall - install vist3d array as eoms array
c
      do 173 i=1,idim1
      do 173 k=1,kdim1
      do 173 j=1,jdim1
      vist3d(j,k,i) = eoms(j,k,i)
 173  continue
c
c   end if on k & j composite:
      end if
c
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-L turb model'',
     .   '' has no laminar regions'')') nblt
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-L turb model -'',
     .   '' laminar region is:'')') nblt
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   i='',i5,'' to'',i5,'', j='',i5,
     .   '' to'',i5,'', k='',i5,'' to'',i5)') ilamlo,ilamhi,jlamlo,
     .   jlamhi,klamlo,klamhi
        end if
      end if
      do 7865 i=ilamlo,ilamhi-1
      do 7865 k=klamlo,klamhi-1
      do 7865 j=jlamlo,jlamhi-1
        vist3d(j,k,i)=0.
 7865 continue
c
      return
c  end if on (ivisc(2) .gt. 1):
      end if      
c
      if (ivisc(1) .le. 1) return
c***************************************************************  
c      evaluate turbulent viscosity along normal to i=0 wall 
c***************************************************************  
c
c   defaults to using min face for distance if neither min nor max contains wall
      if (ibcimin.eq.1 .or. ibcimax.eq.0) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' computing turb viscosity using '',
     . ''Baldwin-Lomax, I=1   , block='',i4)') nblt
      end if
      if (isklton.gt.0 .and. ideg(1).eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1525)
      end if
      if (isklton.gt.0 .and. ideg(1).ne.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1526)
      end if
c
c      loop through k stations
c 
      ihead = 0
c
      if (inmx.lt.idim) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1219)inmx,idim
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
 1219 format(19h error: inmx, idim=,2i5)
c
      do 2000 k=1,kdim1
c
c      loop through j stations
c 
      do 2000 j=1,jdim1
c 
c     create working arrays in I-direction 
c 
c      bound turbulent region
c
      iloop  = .80*idim
      iloop1 = iloop-1
      inmax  = iloop
      inmax1 = inmax-1
c
c      bound search region for fmax
c
      in1    = inmax1*0.20+1
      if(iwf(1) .gt. 0) then
        in1    = 3
      end if
      in2    = inmax1*0.80+1
c
c     chimera scheme modification....limit upper bound of search region to
c     the lower edge of a hole.  if hole extends to lower search bound, set
c     eddy viscosity to zero at all points this station.  
c
      if (iover.eq.1) then
        if(real(blank(j,k,in1)) .eq. 0.) then
          do 1719 i=1,idim1
          eoms(j,k,i) = 0.e0
1719      continue                                                  
          if (iprint.ge.1) then
            if (j.eq.j/ilfreq1*ilfreq1.and.k.eq.k/ilfreq2*ilfreq2) then
              if (ihead.eq.0) then
                 nou(3) = min(nou(3)+1,ibufdim)
                 write(bou(nou(3),3),1409)                           
              end if
              nou(3) = min(nou(3)+1,ibufdim)
              write(bou(nou(3),3),1408) j,k
              ihead = ihead+1
	    end if
          end if
	  go to 2000
        end if
        in21=in2
        do 1723 in=in1,in2
        if(real(blank(j,k,in)) .eq. 0.) go to 1724
        in21 = in
1723    continue
1724    continue    
        in2  = in21
      end if 
c
c      viscosity thru sutherland relation at cell centers 
c
cdir$ ivdep
      do 8108 i=1,iloop1
      t5         =  gamma*q(j,k,i,5)/q(j,k,i,1)
      t6         =  sqrt(t5) 
      amun (i+1) =  c2bp*t5*t6/(c2b+t5)
 8108 continue
cdir$ ivdep
      do 8109 i=1,iloop1
      rhon (i+1) =  q(j,k,i,1)
      vortn(i+1) =  vor(j,k,i)
      disn (i+1) =  ccabs( sni0(j,k,i) )
      utot (i+1) =  sqrt( q(j,k,i,2)*q(j,k,i,2)
     .            +       q(j,k,i,3)*q(j,k,i,3)
     .            +       q(j,k,i,4)*q(j,k,i,4) )
 8109 continue
c
      rhon(1)  =  rhon(2)
      amun(1)  =  amun(2)
      disn(1)  =  0.
      utot(1)  =  utot(2)
      vortn(1) =  vortn(2)
c
      if (real(bci(j,k,1)) .eq. 1.) then
c
c     modify wall values in no-slip region
c
      t5 = gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
      t6 = sqrt(t5)
c
      amun(1) = c2bp*t5*t6/(c2b+t5)
      rhon(1) = qi0(j,k,1,1)
      utot(1) = sqrt(qi0(j,k,2,1)*qi0(j,k,2,1)+
     .               qi0(j,k,3,1)*qi0(j,k,3,1)+
     .               qi0(j,k,4,1)*qi0(j,k,4,1))
c
      vortn(1) = (utot(2)-utot(1))/disn(2)
      end if
c 
c     convenient groupings and initial values at i=1
c 
      yplusoy = ccabs(rhon(1)*vortn(1)*reue*xmi/amun(1) ) 
      yplusoy = sqrt(yplusoy)*aplusi 
      yplus   = yplusoy*disn(2)*aplus*2.0
c
      if (real(bci(j,k,1)) .eq. 0.) then
c        wake damping
         term     = 1.-exp(-50.e0)
         do 8111 in=2,inmax
         damp(in) = term
 8111    continue
      else
c        wall damping 
         do 8112 in=2,inmax
         if (iwall .eq. 0) then
            term = yplusoy*disn(in)*sqrt(rhon(in)/rhon(1))*
     .      (amun(1)/amun(in))
         else
            term = yplusoy*disn(in)
         end if
         damp(in) = 1.0e0 -exp(-term) 
 8112    continue
      end if
c
      fbl(1)    = 0.0e0
      eomui(1)  = 0.0e0
      do 8113 in=2,inmax
      fbl(in)   = vortn(in)*disn(in)*damp(in)
      amixl     = vk*disn(in)*damp(in)
      eomui(in) = reue*xmi*rhon(in)*vortn(in)* 
     .            amixl*amixl
 8113 continue
c
      utmax = q8smax(inmax,utot(1))
      utmin = q8smin(inmax,utot(1))
      if (real(bci(j,k,1)) .eq. 1.) utmin = utot(1)
c 
c     locate max value and location of fbl
c
      fblmax = 1.e-10
      infmax = in1
      do 8115 in=in1,in2
      if (real(fbl(in)).gt.real(fblmax)) then
         fblmax = fbl(in)
         infmax = in
      end if
c
c      first maximum (degani-schiff)
c
      if (ideg(1).eq.1.and.real(fbl(in)).lt.0.9e0*real(fblmax))goto 8125
 8115 continue
 8125 continue
c
c     curve fit of f near infmax to improve fmax location
c
      if (ifit.gt.0) then
         if (infmax.gt.in1 .and. infmax.lt.in2) then
            dfm = fbl(infmax) - fbl(infmax-1)
            dfp = fbl(infmax) - fbl(infmax+1)
            if (real(fbl(infmax-1)).lt.real(fbl(infmax+1))) then
               if (real(dfm).gt.0.) then
                  ymax = disn(infmax) + 0.5*(disn(infmax+1) -
     .                                 disn(infmax))*(1.-dfp/dfm)
               else
                  ymax = disn(infmax)
               end  if
            else
               if (real(dfp).gt.0.) then
                  ymax = disn(infmax) - 0.5*(disn(infmax) -
     .                                 disn(infmax-1))*(1.-dfm/dfp)
               else
                  ymax = disn(infmax)
               end if
            end if
         end if
      else
         ymax   = disn(infmax)
      end if
c
      ymax   = ccmaxcr(ymax,1.e-20)
      udif   = utmax - utmin 
      fwake  = ymax*fblmax 
      fwake2 = cwk*udif*udif*ymax/fblmax 
      fwake  = ccmin(fwake2,fwake)
c
      do 8116 in=2,inmax 
      eomuo(in) = 1.e0/(1.e0+5.5e0*(ckleb*disn(in)/ymax)**6)
      eomuo(in) = ckout*eomuo(in)*rhon(in)*fwake
 8116 continue
      eomuo(1)  =  eomuo(2)
c
      do 8117 in=2,in2
      if (real(eomui(in)).gt.real(eomuo(in))) go to 8127
 8117 continue
 8127 inswtch   = in-1
      do 8118 in=inswtch,inmax
      eomui(in) =  eomuo(in)
 8118 continue
c
c      fill eomu array corresponding to i=0 wall 
c
      do 1900 i=1,iloop1
      eoms(j,k,i) = eomui(i+1)
 1900 continue
      if (iloop1.lt.idim1) then
         do 1807 i=iloop,idim1
         eoms(j,k,i) = 0.
 1807    continue
      end if
c
      if (iprint.ge.1) then
      if (j.eq.(j/ilfreq1)*ilfreq1 .and. k.eq.(k/ilfreq2)*ilfreq2) then
      if (ihead.eq.0) then
         nou(3) = min(nou(3)+1,ibufdim)
         write(bou(nou(3),3),1509)
      end if
 1509 format(1x,3hi=0,4x,1hj,4x,1hk,5x,5hyplus,5x,5hutmax,5x,5hut(2),
     .6x,4hfmax,6x,4hymax,4x,6hemumax,1x,5hswtch,
     .2x,4hfmax,1x,3hin1,1x,3hin2,1x,3hmax)
      ihead = ihead+1
      eomax = q8smax(inmax1,eomui(2))
      utmax = utmax/xmach
      utmin = utmin/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1510) j,k,real(yplus),real(utmax),
     .      real(utmin),real(fblmax),real(ymax),real(eomax),
     .      inswtch,infmax,in1,in2,inmax
 1510 format(i9,i5,3e10.3,3e10.3,2i6,3i4)
      end if
      end if
c
c      profile data
c
      if (iprint.gt.1) then
      if (j.eq.(j/ipfreq1)*ipfreq1 .and. k.eq.(k/ipfreq2)*ipfreq2) then
      ihead = 0
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1513)
 1513 format(4x,1hj,4x,1hk,4x,1hi,6x,4hdisn,6x,4hutot,6x,4heomu,
     .7x,3hfbl,6x,4hamun,4x,6heomu-o)
      do 1515 in=1,inmax
      utw = utot(in)/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),1514) j,k,in,real(disn(in)),real(utw),
     .               real(eomui(in)),real(fbl(in)),
     .               real(amun(in)),real(eomuo(in))
 1514 format(3i5,6e10.3)
 1515 continue
      end if
      end if
 2000 continue
c
c   end if on i=0 face
      end if
c***************************************************************  
c      evaluate turbulent viscosity along normal to i=idim wall 
c***************************************************************  
c
      if (ibcimax .eq. 1) then
c
      if (isklton.gt.0) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),'('' computing turb viscosity using '',
     . ''Baldwin-Lomax, I=idim, block='',i4)') nblt
      end if
      if (isklton.gt.0 .and. ideg(1).eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1525)
      end if
      if (isklton.gt.0 .and. ideg(1).ne.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),1526)
      end if
c
c      loop through k stations
c 
      if (inmx.lt.idim) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),19)inmx,idim
         call termn8(myid,-1,ibufdim,nbuf,bou,nou)
      end if
   19 format(19h error: inmx, idim=,2i5)
c
      ihead = 0
      do 2007 k=1,kdim1
c
c      loop through j stations
c 
      do 2007 j=1,jdim1
c 
c     create working arrays in I-direction 
c 
c      bound turbulent region
c
      iloop  = .80*idim
      iloop1 = iloop-1
      inmaxt = iloop
      inmaxt1= inmaxt-1
c
      iloop  = idim1-iloop+1
      iloop  = max(1,iloop)
      iloop1 = iloop+1
      inmax  = iloop
      inmax1 = inmax+1
c
c      bound search region for fmax
c
      in1    = inmaxt1*0.20+1
      if(iwf(1) .gt. 0) then
        in1    = 3
      end if
      in2    = inmaxt1*0.80+1
      in1    = idim1-in1+1
      in2    = idim1-in2+1
c
c     chimera scheme modification....limit upper bound of search region to
c     the lower edge of a hole.  if hole extends to lower search bound, set
c     eddy viscosity to zero at all points this station.  
c
      if (iover.eq.1) then
        if(real(blank(j,k,in1)) .eq. 0.) then
          do 6719 i=1,idim1
          eoms(j,k,i) = 0.e0
6719      continue                                                  
          if (iprint.ge.1) then
            if (j.eq.j/ilfreq1*ilfreq1.and.k.eq.k/ilfreq2*ilfreq2) then
              if (ihead.eq.0) then
                 nou(3) = min(nou(3)+1,ibufdim)
                 write(bou(nou(3),3),1409)                           
              end if
              nou(3) = min(nou(3)+1,ibufdim)
              write(bou(nou(3),3),1408) j,k
              ihead = ihead+1
	    end if
          end if
	  go to 2007
        end if
        in21=in2
        do 6723 in=in1,in2,-1
        if(real(blank(j,k,in)) .eq. 0.) go to 6724
        in21 = in
6723    continue
6724    continue    
        in2=in21
      end if 
c
c      viscosity thru sutherland relation at cell centers 
c
cdir$ ivdep
      do 8158 i=idim1,iloop1,-1
      t5         =  gamma*q(j,k,i,5)/q(j,k,i,1)
      t6         =  sqrt(t5) 
      amun (i) =  c2bp*t5*t6/(c2b+t5)
 8158 continue
cdir$ ivdep
      do 8159 i=idim1,iloop1,-1
      rhon (i) =  q(j,k,i,1)
      vortn(i) =  vor(j,k,i)
      disn (i) =  ccabs( snim(j,k,i) )
      utot (i) =  sqrt( q(j,k,i,2)*q(j,k,i,2)
     .          +       q(j,k,i,3)*q(j,k,i,3)
     .          +       q(j,k,i,4)*q(j,k,i,4) )
 8159 continue
c
      rhon(idim)  =  rhon(idim1)
      amun(idim)  =  amun(idim1)
      disn(idim)  =  0.
      utot(idim)  =  utot(idim1)
      vortn(idim) =  vortn(idim1)
c
      if (real(bci(j,k,2)) .eq. 1.) then
c
c     modify wall values in no-slip region
c
      t5 = gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
      t6 = sqrt(t5)
c
      amun(idim) = c2bp*t5*t6/(c2b+t5)
      rhon(idim) = qi0(j,k,1,1)
      utot(idim) = sqrt(qi0(j,k,2,3)*qi0(j,k,2,3)+
     .                  qi0(j,k,3,3)*qi0(j,k,3,3)+
     .                  qi0(j,k,4,3)*qi0(j,k,4,3))
c
      vortn(idim) = (utot(idim1)-utot(idim))/disn(idim1)
      end if
c 
c     convenient groupings and initial values at i=1
c 
      yplusoy = ccabs(rhon(idim)*vortn(idim)*reue*xmi/amun(idim) ) 
      yplusoy = sqrt(yplusoy)*aplusi 
      yplus   = yplusoy*disn(idim1)*aplus*2.0
c
      if (real(bci(j,k,2)) .eq. 0.) then
c        wake damping
         term     = 1.-exp(-50.e0)
         do 8171 in=idim1,inmax1,-1
         damp(in) = term
 8171    continue
      else
c        wall damping 
         do 8172 in=idim1,inmax1,-1
         if (iwall .eq. 0) then
            term = yplusoy*disn(in)*sqrt(rhon(in)/rhon(idim))*
     .      (amun(idim)/amun(in))
         else
            term = yplusoy*disn(in)
         end if
         damp(in) = 1.0e0 -exp(-term) 
 8172    continue
      end if
c
      fbl(idim)    = 0.0e0
      eomui(idim)  = 0.0e0
      do 8173 in=idim1,inmax1,-1
      fbl(in)   = vortn(in)*disn(in)*damp(in)
      amixl     = vk*disn(in)*damp(in)
      eomui(in) = reue*xmi*rhon(in)*vortn(in)* 
     .            amixl*amixl
 8173 continue
c
      utmax = q8smax(inmaxt1,utot(inmax1))
      utmin = q8smin(inmaxt1,utot(inmax1))
      if (real(bci(j,k,2)) .eq. 1.) utmin = utot(idim)
c 
c     locate max value and location of fbl
c
      fblmax = 1.e-10
      infmax = in1+1
      do 8975 in=in1+1,in2+1,-1
      if (real(fbl(in)).gt.real(fblmax)) then
         fblmax = fbl(in)
         infmax = in
      end if
c
c      first maximum (degani-schiff)
c
      if (ideg(1).eq.1.and.real(fbl(in)).lt.0.9e0*real(fblmax))goto 8175
 8975 continue
 8175 continue
c
c     curve fit of f near infmax to improve fmax location
c
      if (ifit.gt.0) then
         if (infmax.lt.in1 .and. infmax.gt.in2) then
            dfm = fbl(infmax) - fbl(infmax+1)
            dfp = fbl(infmax) - fbl(infmax-1)
            if (real(fbl(infmax+1)).lt.real(fbl(infmax-1))) then
               if (real(dfm).gt.0.) then
                  ymax = disn(infmax) + 0.5*(disn(infmax-1) -
     .                                 disn(infmax))*(1.-dfp/dfm)
               else
                  ymax = disn(infmax)
               end  if
            else
               if (real(dfp).gt.0.) then
                  ymax = disn(infmax) - 0.5*(disn(infmax) -
     .                                 disn(infmax+1))*(1.-dfm/dfp)
               else
                  ymax = disn(infmax)
               end if
            end if
         end if
      else
         ymax   = disn(infmax)
      end if
c
      ymax   = ccmaxcr(ymax,1.e-20)
      udif   = utmax - utmin 
      fwake  = ymax*fblmax 
      fwake2 = cwk*udif*udif*ymax/fblmax 
      fwake  = ccmin(fwake2,fwake)
c
      do 8176 in=idim1,inmax1,-1
      eomuo(in) = 1.e0/(1.e0+5.5e0*(ckleb*disn(in)/ymax)**6)
      eomuo(in) = ckout*eomuo(in)*rhon(in)*fwake
 8176 continue
      eomuo(idim)  =  eomuo(idim1)
c
      do 8177 in=idim1,in2,-1
      if (real(eomui(in)).gt.real(eomuo(in))) go to 8197
 8177 continue
 8197 inswtch   = in+1
      do 8178 in=inswtch,inmax1,-1
      eomui(in) =  eomuo(in)
 8178 continue
c
c      fill eomu array corresponding to i=idim wall 
c
      if (ibcimin.eq.1) then
c
c     for i=1 & i=idim walls
      do 177 i=idim1,iloop1,-1
      s1 = sni0(j,k,i)*sni0(j,k,i)
      s2 = snim(j,k,i)*snim(j,k,i)
      eoms(j,k,i) = ( s1*eomui(k) + s2*eoms(j,k,i) )/( s1 + s2 )
  177 continue
c
      else
c
c     for i=idim wall only
      do 6900 i=idim1,iloop1,-1
      eoms(j,k,i) = eomui(i)
 6900 continue
      if (iloop1.gt.1) then
         do 6807 i=iloop,1,-1
         eoms(j,k,i) = 0.
 6807    continue
      end if
c
      end if
c
      if (iprint.ge.1) then
      if (j.eq.(j/ilfreq1)*ilfreq1 .and. k.eq.(k/ilfreq2)*ilfreq2) then
      if (ihead.eq.0) then
         nou(3) = min(nou(3)+1,ibufdim)
         write(bou(nou(3),3),6509)
      end if
 6509 format(1x,3hi=0,4x,1hj,4x,1hk,5x,5hyplus,5x,5hutmax,5x,5hut(2),
     .6x,4hfmax,6x,4hymax,4x,6hemumax,1x,5hswtch,
     .2x,4hfmax,1x,3hin1,1x,3hin2,1x,3hmax)
      ihead = ihead+1
      eomax = q8smax(inmax1,eomui(idim1))
      utmax = utmax/xmach
      utmin = utmin/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),6510) j,k,real(yplus),real(utmax),
     .      real(utmin),real(fblmax),real(ymax),real(eomax),
     .      inswtch,infmax,in1,in2,inmax
 6510 format(i9,i5,3e10.3,3e10.3,2i6,3i4)
      end if
      end if
c
c      profile data
c
      if (iprint.gt.1) then
      if (j.eq.(j/ipfreq1)*ipfreq1 .and. k.eq.(k/ipfreq2)*ipfreq2) then
      ihead = 0
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),6513)
 6513 format(4x,1hj,4x,1hk,4x,1hi,6x,4hdisn,6x,4hutot,6x,4heomu,
     .7x,3hfbl,6x,4hamun,4x,6heomu-o)
      do 6515 in=idim1,inmax1,-1
      utw = utot(in)/xmach
      nou(3) = min(nou(3)+1,ibufdim)
      write(bou(nou(3),3),6514) j,k,in,real(disn(in)),real(utw),
     .               real(eomui(in)),real(fbl(in)),
     .               real(amun(in)),real(eomuo(in))
 6514 format(3i5,6e10.3)
 6515 continue
      end if
      end if
 2007 continue
c
c   end if on i=idim face
      end if
c
      if (ivisc(3).gt.1) then
c 
c     form the composite eddy-viscosity from k & i walls
c 
      do 6250 i=1,idim1
      do 6250 k=1,kdim1
      do 6250 j=1,jdim1
      if(ibckmin.eq.0 .and. ibckmax.eq.1) then
        s1            = snkm(j,k,i)**2
      else if(ibckmin.eq.1 .and. ibckmax.eq.1) then
        s1            = ccmin(snk0(j,k,i),snkm(j,k,i))
        s1            = s1**2
      else
        s1            = snk0(j,k,i)**2
      end if
      if(ibcimin.eq.0 .and. ibcimax.eq.1) then
        s2            = snim(j,k,i)**2
      else if(ibcimin.eq.1 .and. ibcimax.eq.1) then
        s2            = ccmin(sni0(j,k,i),snim(j,k,i))
        s2            = s2**2
      else
        s2            = sni0(j,k,i)**2
      end if
      vist3d(j,k,i) =( vist3d(j,k,i)*s2  + eoms(j,k,i)*s1 )/( s1 + s2 )
6250  continue
c
c      smooth interior values
c
      nsmurf = 3 
      if (nsmurf.eq.0) return
c
      do 6285 nnn=1,nsmurf
c
      do 6270 i=1,idim1
      do 6270 k=1,kdim1
      do 6270 j=1,jdim1
      eoms(j,k,i) = vist3d(j,k,i)
6270  continue
c
      do 6276 i=2,idim1-1
      do 6276 k=2,kdim1-1
      do 6276 j=2,jdim1-1
      vist3d(j,k,i) = .125e0* ( eoms(j+1,k-1,i-1)+ eoms(j-1,k-1,i-1)
     .               +          eoms(j+1,k+1,i-1)+ eoms(j-1,k+1,i-1)
     .               +          eoms(j+1,k-1,i+1)+ eoms(j-1,k-1,i+1)
     .               +          eoms(j+1,k+1,i+1)+ eoms(j-1,k+1,i+1)) 
6276  continue
6285  continue
c
      else
c
c      single turbulent wall - install vist3d array as eoms array
c
      do 6273 i=1,idim1
      do 6273 k=1,kdim1
      do 6273 j=1,jdim1
      vist3d(j,k,i) = eoms(j,k,i)
6273  continue
c
c   end if on composite
      end if
c
      if(isklton .gt. 0) then
        if(ilamlo.eq.0 .or. jlamlo.eq.0 .or. klamlo.eq.0) then
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-L turb model'',
     .   '' has no laminar regions'')') nblt
        else
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   block '',i4,'' in B-L turb model -'',
     .   '' laminar region is:'')') nblt
        nou(1) = min(nou(1)+1,ibufdim)
        write(bou(nou(1),1),'(''   i='',i5,'' to'',i5,'', j='',i5,
     .   '' to'',i5,'', k='',i5,'' to'',i5)') ilamlo,ilamhi,jlamlo,
     .   jlamhi,klamlo,klamhi
        end if
      end if
      do 7897 i=ilamlo,ilamhi-1
      do 7897 k=klamlo,klamhi-1
      do 7897 j=jlamlo,jlamhi-1
        vist3d(j,k,i)=0.
 7897 continue
      return
      end 
