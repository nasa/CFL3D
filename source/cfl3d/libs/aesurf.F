      subroutine aesurf(nbl,jdim,kdim,idim,deltj,deltk,delti,xmdj,
     .                  xmdk,xmdi,maxbl,maxseg,nmds,maxaes,aesrfdat,
     .                  xs,xxn,icsi,icsf,jcsi,jcsf,kcsi,kcsf,nsegdfrm,
     .                  idfrmseg,iaesurf,maxsegdg)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Update surface deflections due to aeroelastic motion
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension icsi(maxbl,maxsegdg),icsf(maxbl,maxsegdg),
     .          jcsi(maxbl,maxsegdg),jcsf(maxbl,maxsegdg),
     .          kcsi(maxbl,maxsegdg),kcsf(maxbl,maxsegdg)
      dimension nsegdfrm(maxbl),idfrmseg(maxbl,maxsegdg),
     .          iaesurf(maxbl,maxsegdg)
      dimension aesrfdat(5,maxaes)
      dimension deltj(kdim,idim,3,2),deltk(jdim,idim,3,2),
     .          delti(jdim,kdim,3,2)
      dimension xmdj(kdim,idim,6,nmds,maxaes),
     .          xmdk(jdim,idim,6,nmds,maxaes),
     .          xmdi(jdim,kdim,6,nmds,maxaes),
     .          xs(2*nmds,maxaes),xxn(2*nmds,maxaes)
      dimension nsrfflg(kdim*idim*4+jdim*idim*4+jdim*kdim*4,4)
c
      common /elastic/ ndefrm,naesrf
c
      do iaes=1,naesrf
c
         iskyhk = aesrfdat(1,iaes)
         grefl  = aesrfdat(2,iaes)
         uinf   = aesrfdat(3,iaes)
         qinf   = aesrfdat(4,iaes)
         nmodes = aesrfdat(5,iaes)
         rgrefl = 1./grefl

         is1 = 0
         do is=1,nsegdfrm(nbl)
         if(jcsi(nbl,is) .eq. jcsf(nbl,is)) then
           kinc = kcsf(nbl,is)-kcsi(nbl,is)
           do k = kcsi(nbl,is),kcsf(nbl,is),kinc
             do i = icsi(nbl,is),icsf(nbl,is)
               is1            = is1 + 1
               nsrfflg(is1,1) = i
               nsrfflg(is1,2) = jcsi(nbl,is)
               nsrfflg(is1,3) = k
               nsrfflg(is1,4) = 0
             enddo
           enddo
           iinc = icsf(nbl,is)-icsi(nbl,is)
           do i = icsi(nbl,is),icsf(nbl,is),iinc
             do k = kcsi(nbl,is),kcsf(nbl,is)
               is1            = is1 + 1
               nsrfflg(is1,1) = i
               nsrfflg(is1,2) = jcsi(nbl,is)
               nsrfflg(is1,3) = k
               nsrfflg(is1,4) = 0
             enddo
           enddo
         else
     .   if(icsi(nbl,is) .eq. icsf(nbl,is)) then
           kinc = kcsf(nbl,is)-kcsi(nbl,is)
           do k = kcsi(nbl,is),kcsf(nbl,is),kinc
             do j = jcsi(nbl,is),jcsf(nbl,is)
               is1            = is1 + 1
               nsrfflg(is1,1) = icsi(nbl,is)
               nsrfflg(is1,2) = j
               nsrfflg(is1,3) = k
               nsrfflg(is1,4) = 0
             enddo
           enddo
           jinc = jcsf(nbl,is)-jcsi(nbl,is)
           do j = jcsi(nbl,is),jcsf(nbl,is),jinc
             do k = kcsi(nbl,is),kcsf(nbl,is)
               is1            = is1 + 1
               nsrfflg(is1,1) = icsi(nbl,is)
               nsrfflg(is1,2) = j
               nsrfflg(is1,3) = k
               nsrfflg(is1,4) = 0
             enddo
           enddo
         else
     .   if(kcsi(nbl,is) .eq. kcsf(nbl,is)) then
           iinc = icsf(nbl,is)-icsi(nbl,is)
           do i = icsi(nbl,is),icsf(nbl,is),iinc
             do j = jcsi(nbl,is),jcsf(nbl,is)
               is1            = is1 + 1
               nsrfflg(is1,1) = i
               nsrfflg(is1,2) = j
               nsrfflg(is1,3) = kcsi(nbl,is)
               nsrfflg(is1,4) = 0
             enddo
           enddo
           jinc = jcsf(nbl,is)-jcsi(nbl,is)
           do j = jcsi(nbl,is),jcsf(nbl,is),jinc
             do i = icsi(nbl,is),icsf(nbl,is)
               is1            = is1 + 1
               nsrfflg(is1,1) = i
               nsrfflg(is1,2) = j
               nsrfflg(is1,3) = kcsi(nbl,is)
               nsrfflg(is1,4) = 0
             enddo
           enddo
         end if
         enddo
         nbndry = is1
c
         do is=1,nsegdfrm(nbl)
c
         if (iaesurf(nbl,is).eq.iaes .and. idfrmseg(nbl,is).eq.99) then
c
c           convert modal deflections to x,y,z - j=const surfaces
c
            if (jcsi(nbl,is) .eq. jcsf(nbl,is)) then
c
               j  = jcsi(nbl,is)
               ll = 0
               m  = 1
               if (jcsi(nbl,is) .eq.jdim) then
                  ll = 3
                  m  = 2
               end if
               ist = icsi(nbl,is)
               ifn = icsf(nbl,is)
               kst = kcsi(nbl,is)
               kfn = kcsf(nbl,is)
c
               do 3020 i = ist,ifn
               do 3020 k = kst,kfn
c
c              Flag point to indicate it has already been displaced, 
c              to eliminate doubling of displacement
c
               do ii = 1,nbndry
                   if(i.eq.nsrfflg(ii,1).and.j.eq.nsrfflg(ii,2).and.
     .                k.eq.nsrfflg(ii,3)) then
                      if(nsrfflg(ii,4).eq.1) goto 3015
                      nsrfflg(ii,4) = 1
                      goto 3009
                   end if
               enddo
3009           continue
c              save any existing displacements
               deltj1 = deltj(k,i,1,m)
               deltj2 = deltj(k,i,2,m)
               deltj3 = deltj(k,i,3,m)
c              initialize modal displacements
               deltj(k,i,1,m) = 0.
               deltj(k,i,2,m) = 0.
               deltj(k,i,3,m) = 0.
c              sum delta contributions from all modes
               do 3010 n = 1,nmodes
               deltj(k,i,1,m) = deltj(k,i,1,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdj(k,i,ll+1,n,iaes)*rgrefl
               deltj(k,i,2,m) = deltj(k,i,2,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdj(k,i,ll+2,n,iaes)*rgrefl
               deltj(k,i,3,m) = deltj(k,i,3,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdj(k,i,ll+3,n,iaes)*rgrefl
3010           continue
c              add modal displacements to existing displacements
               deltj(k,i,1,m) = deltj(k,i,1,m) + deltj1
               deltj(k,i,2,m) = deltj(k,i,2,m) + deltj2
               deltj(k,i,3,m) = deltj(k,i,3,m) + deltj3
3015           continue
3020           continue
3030           continue
c
            end if
c
c           convert modal deflections to x,y,z - k=const surfaces
c
            if (kcsi(nbl,is) .eq. kcsf(nbl,is)) then
c
               k  = kcsi(nbl,is)
               ll = 0
               m  = 1
               if (kcsi(nbl,is) .eq.kdim) then
                  ll = 3
                  m  = 2
               end if
               ist = icsi(nbl,is)
               ifn = icsf(nbl,is)
               jst = jcsi(nbl,is)
               jfn = jcsf(nbl,is)
c
               do 3050 i = ist,ifn
               do 3050 j = jst,jfn
c
c              Flag point to indicate it has already been displaced, 
c              to eliminate doubling of displacement
c
               do ii = 1,nbndry
                   if(i.eq.nsrfflg(ii,1).and.j.eq.nsrfflg(ii,2).and.
     .                k.eq.nsrfflg(ii,3)) then
                      if(nsrfflg(ii,4).eq.1) goto 3045
                      nsrfflg(ii,4) = 1
                      goto 3039
                   end if
               enddo
3039           continue
c              save any existing displacements
               deltk1 = deltk(j,i,1,m)
               deltk2 = deltk(j,i,2,m)
               deltk3 = deltk(j,i,3,m)
c              initialize modal displacements
               deltk(j,i,1,m) = 0.
               deltk(j,i,2,m) = 0.
               deltk(j,i,3,m) = 0.
c              sum delta contributions from all modes
               do 3040 n = 1,nmodes
               deltk(j,i,1,m) = deltk(j,i,1,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdk(j,i,ll+1,n,iaes)*rgrefl
               deltk(j,i,2,m) = deltk(j,i,2,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdk(j,i,ll+2,n,iaes)*rgrefl
               deltk(j,i,3,m) = deltk(j,i,3,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdk(j,i,ll+3,n,iaes)*rgrefl
3040           continue
c              add modal displacements to existing displacements
               deltk(j,i,1,m)=deltk(j,i,1,m) + deltk1
               deltk(j,i,2,m)=deltk(j,i,2,m) + deltk2
               deltk(j,i,3,m)=deltk(j,i,3,m) + deltk3
3045           continue
3050           continue
3060           continue
c
            end if
c
c           convert modal deflections to x,y,z - i=const surfaces
c
            if (icsi(nbl,is) .eq. icsf(nbl,is)) then
c
               i  = icsi(nbl,is)
               ll = 0
               m  = 1
               if (icsi(nbl,is) .eq.idim) then
                  ll = 3
                  m  = 2
               end if
               kst = kcsi(nbl,is)
               kfn = kcsf(nbl,is)
               jst = jcsi(nbl,is)
               jfn = jcsf(nbl,is)
c
               do 3080 k = kst,kfn
               do 3080 j = jst,jfn
c
c              Flag point to indicate it has already been displaced, 
c              to eliminate doubling of displacement
c
               do ii = 1,nbndry
                   if(i.eq.nsrfflg(ii,1).and.j.eq.nsrfflg(ii,2).and.
     .                k.eq.nsrfflg(ii,3)) then
                      if(nsrfflg(ii,4).eq.1) goto 3075
                      nsrfflg(ii,4) = 1
                      goto 3069
                   end if
               enddo
3069           continue
c              save any existing displacements
               delti1 = delti(j,k,1,m)
               delti2 = delti(j,k,2,m)
               delti3 = delti(j,k,3,m)
c              initialize modal displacements
               delti(j,k,1,m) = 0.
               delti(j,k,2,m) = 0.
               delti(j,k,3,m) = 0.
c              sum delta contributions from all modes
               do 3070 n = 1,nmodes
               delti(j,k,1,m) = delti(j,k,1,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdi(j,k,ll+1,n,iaes)*rgrefl
               delti(j,k,2,m) = delti(j,k,2,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdi(j,k,ll+2,n,iaes)*rgrefl
               delti(j,k,3,m) = delti(j,k,3,m) + (xs(2*n-1,iaes)
     .                - xxn(2*n-1,iaes))*xmdi(j,k,ll+3,n,iaes)*rgrefl
3070           continue
c              add modal displacements to existing displacements
               delti(j,k,1,m)=delti(j,k,1,m) + delti1
               delti(j,k,2,m)=delti(j,k,2,m) + delti2
               delti(j,k,3,m)=delti(j,k,3,m) + delti3
3075           continue
3080           continue
3090           continue
c
            end if
c
         end if
c
         end do
c
      end do
c
      return
      end
