      subroutine wrestg(nbl,jdim,kdim,idim,x,y,z,xnm2,ynm2,znm2,    
     .                  deltj,deltk,delti,qc0,nflagg,iuns,utrans,
     .                  vtrans,wtrans,omegax,omegay,omegaz,xorig,
     .                  yorig,zorig,dxmx,dymx,dzmx,dthxmx,dthymx,
     .                  dthzmx,thetax,thetay,thetaz,rfreqt,
     .                  rfreqr,xorig0,yorig0,zorig0,time2,
     .                  thetaxl,thetayl,thetazl,itrans,irotat,idefrm,
     .                  utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,omgzae,
     .                  xorgae,yorgae,zorgae,thtxae,thtyae,thtzae,
     .                  rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                  kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                  aesrfdat,perturb,myhost,myid,mycomm,mblk2nd,
     .                  maxbl,nsegdfrm,idfrmseg,iaesurf,maxsegdg,
c   ***CGNSstart
c    .                  wk,nwork)
     .                  wk,nwork,idima,jdima,kdima,igrid)
c   ***CGNSend
c***********************************************************************
c     Purpose:  Append to the end of the restart file the required info 
c     for a dynamic mesh (nflagg=1); write qc0 for 2nd order accurate
c     restart (in time) (nflagg=0).
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
#if defined DIST_MPI
#     include "mpif.h"
      dimension istat(MPI_STATUS_SIZE)
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
#endif
c
      dimension wk(nwork)
      dimension x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim)
      dimension xnm2(jdim,kdim,idim),ynm2(jdim,kdim,idim),
     .          znm2(jdim,kdim,idim)
      dimension qc0(jdim,kdim,idim-1,5)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension utrans(maxbl),vtrans(maxbl),wtrans(maxbl),
     .          omegax(maxbl),omegay(maxbl),omegaz(maxbl),
     .          xorig(maxbl),yorig(maxbl),zorig(maxbl),
     .          thetax(maxbl),thetay(maxbl),thetaz(maxbl),
     .          rfreqt(maxbl),rfreqr(maxbl),xorig0(maxbl),
     .          yorig0(maxbl),zorig0(maxbl),time2(maxbl),
     .          thetaxl(maxbl),thetayl(maxbl),thetazl(maxbl),
     .          itrans(maxbl),irotat(maxbl),idefrm(maxbl)
      dimension dxmx(maxbl),dymx(maxbl),dzmx(maxbl),dthxmx(maxbl),
     .          dthymx(maxbl),dthzmx(maxbl)
      dimension mblk2nd(maxbl)
      dimension utrnsae(maxbl,maxsegdg),vtrnsae(maxbl,maxsegdg),
     .          wtrnsae(maxbl,maxsegdg),omgxae(maxbl,maxsegdg),
     .          omgyae(maxbl,maxsegdg),omgzae(maxbl,maxsegdg),
     .          xorgae(maxbl,maxsegdg),yorgae(maxbl,maxsegdg),
     .          zorgae(maxbl,maxsegdg),thtxae(maxbl,maxsegdg),
     .          thtyae(maxbl,maxsegdg),thtzae(maxbl,maxsegdg),
     .          rfrqtae(maxbl,maxsegdg),rfrqrae(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension freq(nmds,maxaes),gmass(nmds,maxaes),x0(2*nmds,maxaes),
     .          gf0(2*nmds,maxaes),damp(nmds,maxaes),
     .          perturb(nmds,maxaes,4)
      dimension aesrfdat(5,maxaes)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
      common /motionmc/ xmc0,ymc0,zmc0,utransmc,vtransmc,wtransmc,
     .                  omegaxmc,omegaymc,omegazmc,xorigmc,yorigmc,
     .                  zorigmc,xorig0mc,yorig0mc,zorig0mc,thetaxmc,
     .                  thetaymc,thetazmc,dxmxmc,dymxmc,dzmxmc,
     .                  dthxmxmc,dthymxmc,dthzmxmc,rfreqtmc,
     .                  rfreqrmc,itransmc,irotatmc,time2mc
      common /unst/ time,cfltau,ntstep,ita,iunst
c   ***CGNSstart
      common /cgns/ icgns,iccg,ibase,nzones,nsoluse,irind,jrind,krind
      common /twod/ i2d
c   ***CGNSend
c
      jdim1 = jdim-1
      kdim1 = kdim-1
      idim1 = idim-1
c
#if defined DIST_MPI
      nd_srce = mblk2nd(nbl)
c
c     set baseline tag values
c
      ioffset    = maxbl
      itag_qc0   = 1
      itag_dmdat = itag_qc0   + ioffset
      itag_x     = itag_dmdat + ioffset
      itag_y     = itag_x     + ioffset
      itag_z     = itag_y     + ioffset
      itag_xnm2  = itag_z     + ioffset
      itag_ynm2  = itag_xnm2  + ioffset
      itag_znm2  = itag_ynm2  + ioffset
      itag_deltj = itag_znm2  + ioffset
      itag_deltk = itag_deltj + ioffset
      itag_delti = itag_deltk + ioffset
#endif
      if (nflagg .eq. 0) then
c
c ************************************
c        write out 2nd order time data 
c ************************************
c
#if defined DIST_MPI
c
c        send/receive second-order time data qc0 to/on
c        the appropriate node
c
         jkim5 = jdim*kdim*idim1*5
         mytag = itag_qc0 + nbl
         if (myid .eq. mblk2nd(nbl)) then
            call MPI_Send(qc0,jkim5,MY_MPI_REAL,
     .                    myhost,mytag,mycomm,ierr)
         else if (myid .eq. myhost) then
            call MPI_Recv(qc0,jkim5,MY_MPI_REAL,
     .                    nd_srce,mytag,mycomm,istat,ierr)
         end if
c
#endif
         if (myid.eq.myhost) then
c   ***CGNSstart
         if (icgns .ne. 1) then
c   ***CGNSend
            write(11,101) nbl
  101       format(45h writing 2nd order time data to restart file,,
     .             6h block,i4)
c
            write(2) jdim,kdim,idim
            write(2) ((((qc0(j,k,i,l),j=1,jdim1),k=1,kdim1),i=1,idim1),
     .                   l=1,5)
c   ***CGNSstart
         else
#if defined CGNS
         write(11,'('' writing 2nd order time data to CGNS file,'',
     +     '' block='',i4)') nbl
         call wsecord(iccg,ibase,igrid,idima,jdima,kdima,idim,
     +                jdim,kdim,wk,qc0,nsoluse,i2d)
#endif
         end if
c   ***CGNSend
         end if
c
      else if (nflagg .eq. 1) then
c
c ************************************
c        write out dynamic mesh data
c ************************************
c
         if (myid.eq.myhost) then
c   ***CGNSstart
         if (icgns .ne. 1) then
c   ***CGNSend
            write(2) iuns
c   ***CGNSstart
         end if
c   ***CGNSend
         end if
c
         if (iuns .ne. 0) then
c
            if (myid.eq.myhost) then
c   ***CGNSstart
            if (icgns .ne. 1) then
c   ***CGNSend
               write(11,102) nbl
  102          format(43h writing dynamic mesh data to restart file,,
     .             8h   block,i4)
c
               write(2) jdim,kdim,idim
c   ***CGNSstart
            end if
c   ***CGNSend
            end if
c
#if defined DIST_MPI
            if (myid.eq.mblk2nd(nbl)) then
               idir = 0
               call mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                    omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                    dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                    thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                    time2,thetaxl,thetayl,thetazl,itrans,irotat,
     .                    idefrm,utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                    omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                    thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                    maxsegdg,wk,nwork)
               mytag = itag_dmdat + nbl
               nval  = 58 + 20*nsegdfrm(nbl)
               call MPI_Send (wk, nval, MY_MPI_REAL,
     .                        myhost,mytag,mycomm,ierr)
            else if (myid.eq.myhost) then
               mytag = itag_dmdat + nbl
               nval  = 58 + 20*nsegdfrm(nbl)
               call MPI_Recv (wk, nval, MY_MPI_REAL,
     .                  nd_srce, mytag, mycomm, istat, ierr)
               idir = 1
               call mvdat(nbl,idir,maxbl,utrans,vtrans,wtrans,omegax,
     .                    omegay,omegaz,xorig,yorig,zorig,dxmx,dymx,
     .                    dzmx,dthxmx,dthymx,dthzmx,thetax,thetay,
     .                    thetaz,rfreqt,rfreqr,xorig0,yorig0,zorig0,
     .                    time2,thetaxl,thetayl,thetazl,itrans,irotat,
     .                    idefrm,utrnsae,vtrnsae,wtrnsae,omgxae,omgyae,
     .                    omgzae,xorgae,yorgae,zorgae,thtxae,thtyae,
     .                    thtzae,rfrqtae,rfrqrae,icsi,icsf,jcsi,jcsf,
     .                    kcsi,kcsf,freq,gmass,damp,x0,gf0,nmds,maxaes,
     .                    aesrfdat,perturb,nsegdfrm,idfrmseg,iaesurf,
     .                    maxsegdg,wk,nwork)
#endif
               if (myid.eq.myhost) then
c   ***CGNSstart
               if (icgns .ne. 1) then
c   ***CGNSend
                  write(2) itrans(nbl),rfreqt(nbl),xorig(nbl),
     .                  yorig(nbl),zorig(nbl),xorig0(nbl),yorig0(nbl),
     .                  zorig0(nbl),utrans(nbl),vtrans(nbl),wtrans(nbl),
     .                  dxmx(nbl),dymx(nbl),dzmx(nbl),itransmc,rfreqtmc,
     .                  xorigmc,yorigmc,zorigmc,xorig0mc,yorig0mc,
     .                  zorig0mc,utransmc,vtransmc,wtransmc,xmc,ymc,zmc,
     .                  dxmxmc,dymxmc,dzmxmc,irotat(nbl),rfreqr(nbl),
     .                  thetax(nbl),thetay(nbl),thetaz(nbl),omegax(nbl),
     .                  omegay(nbl),omegaz(nbl),dthxmx(nbl),dthymx(nbl),
     .                  dthzmx(nbl),irotatmc,rfreqrmc,thetaxmc,thetaymc,
     .                  thetazmc,omegaxmc,omegaymc,omegazmc,dthxmxmc,
     .                  dthymxmc,dthzmxmc,time2(nbl),time2mc,dt
c   ***CGNSstart
               else
#if defined CGNS
              write(11,'('' writing dynamic mesh data to CGNS file,'',
     +         '' block'',i4)') nbl
              call wgrdmov(iccg,ibase,igrid,iuns,
     .                  itrans(nbl),rfreqt(nbl),xorig(nbl),yorig(nbl),
     .                  zorig(nbl),xorig0(nbl),yorig0(nbl),
     .                  zorig0(nbl),utrans(nbl),vtrans(nbl),
     .                  wtrans(nbl),dxmx(nbl),dymx(nbl),dzmx(nbl),
     .                  itransmc,rfreqtmc,xorigmc,yorigmc,zorigmc,
     .                  xorig0mc,yorig0mc,zorig0mc,utransmc,vtransmc,
     .                  wtransmc,xmc,ymc,zmc,dxmxmc,dymxmc,dzmxmc,
     .                  irotat(nbl),rfreqr(nbl),thetax(nbl),
     .                  thetay(nbl),thetaz(nbl),omegax(nbl),omegay(nbl),
     .                  omegaz(nbl),dthxmx(nbl),dthymx(nbl),dthzmx(nbl),
     .                  irotatmc,rfreqrmc,thetaxmc,thetaymc,thetazmc,
     .                  omegaxmc,omegaymc,omegazmc,dthxmxmc,dthymxmc,
     .                  dthzmxmc,time2(nbl),time2mc,dt)
#endif
               end if
c   ***CGNSend
               end if
#if defined DIST_MPI
            end if
#endif 
            if (idefrm(nbl) .gt. 0) then
c
c
c              the restart file is written out at time step n; we need
c              to reconstruct the grid at time step n, and for second
c              order time accurate solutions, we also need the grid
c              at time n-1 (first order does not require the grid at
c              n-1). However, for the first order case the grid
c              at n is written out twice, so that the restart file
c              always contains the same amount of data.
c
c
#if defined DIST_MPI
               if (abs(ita).gt.1) then
                  jki = jdim*kdim*idim
                  if (myid.eq.mblk2nd(nbl)) then
                     mytag = itag_xnm2 + nbl
                     call MPI_Send (xnm2, jki, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                     mytag = itag_ynm2 + nbl
                     call MPI_Send (ynm2, jki, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                     mytag = itag_znm2 + nbl
                     call MPI_Send (znm2, jki, MY_MPI_REAL,
     .                              myhost, mytag, mycomm, ierr)
                  else if (myid.eq.myhost) then
                     mytag = itag_xnm2 + nbl
                     call MPI_Recv (xnm2, jki, MY_MPI_REAL,
     .                        nd_srce, mytag, mycomm, istat, ierr)
                     mytag = itag_ynm2 + nbl
                     call MPI_Recv (ynm2, jki, MY_MPI_REAL,
     .                        nd_srce, mytag, mycomm, istat, ierr)
                     mytag = itag_znm2 + nbl
                     call MPI_Recv (znm2, jki, MY_MPI_REAL,
     .                        nd_srce, mytag, mycomm, istat, ierr)
                  end if
               end if
               jki = jdim*kdim*idim
               if (myid.eq.mblk2nd(nbl)) then
                  mytag = itag_x + nbl
                  call MPI_Send (x, jki, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, ierr)
                  mytag = itag_y + nbl
                  call MPI_Send (y, jki, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, ierr)
                  mytag = itag_z + nbl
                  call MPI_Send (z, jki, MY_MPI_REAL,
     .                           myhost, mytag, mycomm, ierr)
               else if (myid.eq.myhost) then
                  mytag = itag_x + nbl
                  call MPI_Recv (x, jki, MY_MPI_REAL,
     .                     nd_srce, mytag, mycomm, istat, ierr)
                  mytag = itag_y + nbl
                  call MPI_Recv (y, jki, MY_MPI_REAL,
     .                     nd_srce, mytag, mycomm, istat, ierr)
                  mytag = itag_z + nbl
                  call MPI_Recv (z, jki, MY_MPI_REAL,
     .                     nd_srce, mytag, mycomm, istat, ierr)
               end if
#endif
               if (myid.eq.myhost) then
c   ***CGNSstart
               if (icgns .ne. 1) then
c   ***CGNSend
                  write(2) idefrm(nbl),nsegdfrm(nbl)
                  do is=1,nsegdfrm(nbl)
                     write(2) idfrmseg(nbl,is),utrnsae(nbl,is),
     .                        vtrnsae(nbl,is),wtrnsae(nbl,is),
     .                        omgxae(nbl,is),omgyae(nbl,is),
     .                        omgzae(nbl,is),xorgae(nbl,is),
     .                        yorgae(nbl,is),zorgae(nbl,is),
     .                        thtxae(nbl,is),thtyae(nbl,is),
     .                        thtzae(nbl,is),rfrqtae(nbl,is),
     .                        rfrqrae(nbl,is),icsi(nbl,is),
     .                        icsf(nbl,is),jcsi(nbl,is),jcsf(nbl,is),
     .                        kcsi(nbl,is),kcsf(nbl,is)
                  end do
c
c                 xnm2, etc. only needed (and stored) for 2nd order, so 
c                 write out x,y,z twice in 1st order case in order to 
c                 simplify restart logic
c                 
                  if (abs(ita).gt.1) then
                     write(2) 
     .                     (((xnm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                     (((ynm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                     (((znm2(j,k,i),j=1,jdim),k=1,kdim),i=1,idim)
                  else
                     write(2) 
     .                     (((x(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                     (((y(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                     (((z(j,k,i),j=1,jdim),k=1,kdim),i=1,idim)
                  end if
                  write(2) (((x(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                     (((y(j,k,i),j=1,jdim),k=1,kdim),i=1,idim),
     .                     (((z(j,k,i),j=1,jdim),k=1,kdim),i=1,idim)
c   ***CGNSstart
               else
#if defined CGNS
                write(11,'('' writing deforming mesh data to CGNS'',
     .           '' file, block'',i4)') nbl
                if (abs(ita).gt.1) then
                call wdeform(iccg,ibase,igrid,idefrm(nbl),
     .            nsegdfrm(nbl),idfrmseg(nbl,1),utrnsae(nbl,1),
     .            vtrnsae(nbl,1),wtrnsae(nbl,1),omgxae(nbl,1),
     .            omgyae(nbl,1),omgzae(nbl,1),xorgae(nbl,1),
     .            yorgae(nbl,1),zorgae(nbl,1),thtxae(nbl,1),
     .            thtyae(nbl,1),thtzae(nbl,1),rfrqtae(nbl,1),
     .            rfrqrae(nbl,1),icsi(nbl,1),icsf(nbl,1),jcsi(nbl,1),
     .            jcsf(nbl,1),kcsi(nbl,1),kcsf(nbl,1),jdim,kdim,idim,
     .            x,y,z,xnm2,ynm2,znm2,wk)
                else
                call wdeform(iccg,ibase,igrid,idefrm(nbl),
     .            nsegdfrm(nbl),idfrmseg(nbl,1),utrnsae(nbl,1),
     .            vtrnsae(nbl,1),wtrnsae(nbl,1),omgxae(nbl,1),
     .            omgyae(nbl,1),omgzae(nbl,1),xorgae(nbl,1),
     .            yorgae(nbl,1),zorgae(nbl,1),thtxae(nbl,1),
     .            thtyae(nbl,1),thtzae(nbl,1),rfrqtae(nbl,1),
     .            rfrqrae(nbl,1),icsi(nbl,1),icsf(nbl,1),jcsi(nbl,1),
     .            jcsf(nbl,1),kcsi(nbl,1),kcsf(nbl,1),jdim,kdim,idim,
     .            x,y,z,x,y,z,wk)
                end if
#endif
               end if
c   ***CGNSend
               end if
            end if
         end if
      end if
c
      return
      end
