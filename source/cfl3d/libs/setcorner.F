      subroutine setcorner(j,k,i,xnm1,ynm1,znm1,dx,dy,dz,jdim,kdim,
     .                     idim,wk,nsurf,iflag,ivert,slavept,nslave,
     .                     nou,bou,ibufdim,nbuf,myid,ibl,nmaster,iseq)
c***********************************************************************
c     Purpose:  Compute dispacements for subgrid vertex points to use
c     for transfinite interpolation. Points off the body are given 
c     displacements that are reduced via a decay function from the 
c     displacements of the closest point on the surface. The decay
c     function is the one proposed by Peter Hartwich, with minor
c     modifications to the parameters, and to allow for more than 
c     one master point to influence the slave point (originally, only
c     the nearest master point was used; now, the nmaster nearest master
c     points can be used)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension xnm1(jdim,kdim,idim),ynm1(jdim,kdim,idim),
     .          znm1(jdim,kdim,idim)
      dimension dx(jdim,kdim,idim),dy(jdim,kdim,idim),
     .          dz(jdim,kdim,idim)
      dimension wk(9*nsurf),slavept(nslave,3,nmaster,5),ibl(nsurf)
c
      common /zero/ iexp
      common /deformz/ beta1,ismooth,negvol
      common /fsum/ sref,cref,bref,xmc,ymc,zmc
c
      ivert = ivert + 1
c
c     first time through, determine closest surface point to the subface
c     corner point j,k,i and store for subsequent use
c
      if (iflag.eq.1) then
         if (ivert.gt.nslave) then
            nou(1) = min(nou(1)+1,ibufdim)
            write(bou(nou(1),1),'(''stopping...parameter nslave '',
     .                            ''too small'')')
            call termn8(myid,-1,ibufdim,nbuf,bou,nou)
         end if
         do n=1,nsurf
            ibl(n) = 1
         end do
         do nm=1,nmaster
            dsmin = 1.e20
            lsmin = 0
            ll    = 0
            do n=1,nsurf
               if (ibl(n) .eq. 1) then
                  dssq = (xnm1(j,k,i)-wk(ll+7))**2 
     .                 + (ynm1(j,k,i)-wk(ll+8))**2
     .                 + (znm1(j,k,i)-wk(ll+9))**2
                  if (real(dssq) .lt. real(dsmin)) then
                     dsmin = dssq
                     lsmin = ll
                     nsmin = n
                  end if
                  if (real(dssq) .eq. real(dsmin)) ibl(n) = 0
               end if
               ll = ll + 9
            end do
            slavept(ivert,1,nm,iseq) = float(lsmin)
            slavept(ivert,2,nm,iseq) = dsmin
            if (nm .eq. 1) then
               slavept(ivert,3,nm,iseq) = 1.
            else
               slavept(ivert,3,nm,iseq) = 
     .         slavept(ivert,2,1,iseq)/slavept(ivert,2,nm,iseq)
            end if
         end do
         wttot = 0.
         do nm=1,nmaster
            wttot = wttot + slavept(ivert,3,nm,iseq)
         end do
         do nm=1,nmaster
            slavept(ivert,3,nm,iseq) = slavept(ivert,3,nm,iseq)/wttot
         end do
         return
      end if
c
c     Hartwich damping function, modified to allow distance-weighted 
c     average of nmaster master points for each slave point
c
      eps     = 10.**(-iexp+1)
      dampmax = -log(10.**(-iexp))
      termx   = 0.
      termy   = 0.
      termz   = 0.
c
      do nm=1,nmaster
         lsmin  = int(slavept(ivert,1,nm,iseq))
         dsmin  = slavept(ivert,2,nm,iseq)
         wt     = slavept(ivert,3,nm,iseq)
         dv     = dsmin
c        dm     = wk(lsmin+4)**2 + wk(lsmin+5)**2 + wk(lsmin+6)**2
         dm     = sref
         damp1  = ccabs(beta1)*dv/(eps+dm)
         damp1  = ccmin(dampmax,damp1)
         factor = exp(-damp1)
         termx  = termx + wt*factor*wk(lsmin+4)
         termy  = termy + wt*factor*wk(lsmin+5)
         termz  = termz + wt*factor*wk(lsmin+6)
      end do
c
      dx(j,k,i) = termx
      dy(j,k,i) = termy
      dz(j,k,i) = termz
c
      return
      end
