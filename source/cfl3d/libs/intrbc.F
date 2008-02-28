      subroutine intrbc(q,jdim,kdim,idim,nbl,ldim,maxbl,iitot,lig,
     .                  iipntsg,dxintg,dyintg,dzintg,iiig,jjig,kkig,
     .                  qb,qj0,qk0,qi0,qq,bcj,bck,bci,nou,bou,nbuf,
     .                  ibufdim,icorr,iset)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Interpolate corrections for boundary values for all 
c     grids overlapped in the current mesh.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension q(jdim,kdim,idim,ldim),qb(iitot,5,3),lig(maxbl),
     .          iipntsg(maxbl),dxintg(iitot),dyintg(iitot),
     .          dzintg(iitot),iiig(iitot),jjig(iitot),kkig(iitot)
      dimension qi0(jdim,kdim,ldim,4),qj0(kdim,idim-1,ldim,4),
     .          qk0(jdim,idim-1,ldim,4)
      dimension qq(jdim+1,kdim+1,idim+1,ldim)
      dimension bcj(kdim,idim-1,2),bck(jdim,idim-1,2),bci(jdim,kdim,2)
c
      common /sklton/ isklton
      common /twod/ i2d
      common /zero/iexp
c
c     10.**(-iexp) is machine zero
c
      qlimit  = 10.**(-iexp)
c
c   Note: iset is index for loading into appropriate qb array
c        (1 for q, 2 for vk0,vj0,vi0, 3 for tj0,tk0,ti0)
c
      lsta = lig(nbl)
      lend = lsta+iipntsg(nbl)-1
      if (iipntsg(nbl) .eq. 0) return
c
c     qq is "augmented" q array of dimensions jdim+1 x kdim+1 x idim+1
c     containing cell-center data at interior points and cell-face 
c     center data at boundary points
c
      call augmntq(q,jdim,kdim,idim,nbl,ldim,qj0,qk0,qi0,qq,
     .             bcj,bck,bci,nou,bou,nbuf,ibufdim,icorr)
c
      if (i2d.eq.0) then
         do 12 n=1,ldim
         do 12 l=lsta,lend
c
c        set up interpolation coefficients 
c
         s1 = qq(jjig(l)+1, kkig(l)+1, iiig(l)+1, n)
         s2 = qq(jjig(l)+1, kkig(l)+1, iiig(l)+2, n)
         s3 = qq(jjig(l)+2, kkig(l)+1, iiig(l)+2, n)
         s4 = qq(jjig(l)+2, kkig(l)+1, iiig(l)+1, n)
         s5 = qq(jjig(l)+1, kkig(l)+2, iiig(l)+1, n)
         s6 = qq(jjig(l)+1, kkig(l)+2, iiig(l)+2, n)
         s7 = qq(jjig(l)+2, kkig(l)+2, iiig(l)+2, n)
         s8 = qq(jjig(l)+2, kkig(l)+2, iiig(l)+1, n)
c
         a1 =  s1 
         a2 = -s1+s2 
         a3 = -s1+s4 
         a4 = -s1+s5 
         a5 =  s1-s2+s3-s4
         a6 =  s1-s2-s5+s6
         a7 =  s1-s4-s5+s8
         a8 = -s1+s2-s3+s4+s5-s6+s7-s8 
c
c        interpolate and store in qb array 
c
         qb(l,n,iset) = a1 + a2*dxintg(l) 
     .                + a3*dyintg(l)
     .                + a4*dzintg(l)
     .                + a5*dxintg(l)*dyintg(l)
     .                + a6*dxintg(l)*dzintg(l)
     .                + a7*dyintg(l)*dzintg(l)
     .                + a8*dxintg(l)*dyintg(l)*dzintg(l) 
   12    continue
      else
         do 121 n=1,ldim
         do 121 l=lsta,lend
c
c        set up interpolation coefficients 
c
         s1 = qq(jjig(l)+1, kkig(l)+1, 2,n)
         s4 = qq(jjig(l)+2, kkig(l)+1, 2,n)
         s5 = qq(jjig(l)+1, kkig(l)+2, 2,n)
         s8 = qq(jjig(l)+2, kkig(l)+2, 2,n)
c 
         a3 = -s1+s4 
         a4 = -s1+s5 
         a7 =  s1-s4-s5+s8
c
c        interpolate and store in qb array 
c
         qb(l,n,iset) = s1+a3*dyintg(l)+a4*dzintg(l)+
     .                     a7*dyintg(l)*dzintg(l)
  121    continue
      end if
c
c     limit negative rho,p,vist3d,turb values to a small positive value
c
      if (iset.eq.1) then
         do n=1,5,4
            do l=lsta,lend
               if (real(qb(l,n,iset)) .lt. 0.) qb(l,n,iset) = qlimit
            end do
         end do
      else
         do n=1,ldim
            do l=lsta,lend
               if (real(qb(l,n,iset)) .lt. 0.) then
                  qb(l,n,iset) = qlimit
               end if
            end do
         end do
      end if
      if (isklton.eq.1 .and. iset.eq.1) then
         nou(1) = min(nou(1)+1,ibufdim)
         write(bou(nou(1),1),*)'  qb points updated = ',iipntsg(nbl)
      end if
      return
      end 
