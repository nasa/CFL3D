      subroutine rotateqb(nbl,dthtx,dthty,dthtz,maxbl,iitot,ibcg,
     .                    lig,lbg,ibpntsg,iipntsg,qb)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Rotate solution in qb array through angle dthtx/y/z for 
c     chimera scheme with rotating grids
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension ibcg(iitot),lig(maxbl),lbg(maxbl),ibpntsg(maxbl,4),
     .          iipntsg(maxbl),qb(iitot,5,3)
c
      iset=1
c
c     lsta = lbg(nbl)
c     lend = lsta-1
      lsta = lig(nbl)
      lend = lsta+iipntsg(nbl)-1
c
      if (abs(real(dthtx)) .gt. 0.) then
c
c        rotate qb about an axis parallel to the x-axis
c
         ca = cos(dthtx)
         sa = sin(dthtx)
c
c        if (ibpntsg(nbl,1).gt.0) then
c           lend = lsta+ibpntsg(nbl,1)-1
            do 10 l=lsta,lend
c           qb3 = qb(ibcg(l),3,iset)
c           qb(ibcg(l),3,iset) = qb3*ca - qb(ibcg(l),4,iset)*sa
c           qb(ibcg(l),4,iset) = qb3*sa + qb(ibcg(l),4,iset)*ca
            qb3 = qb(l,3,iset)
            qb(l,3,iset) = qb3*ca - qb(l,4,iset)*sa
            qb(l,4,iset) = qb3*sa + qb(l,4,iset)*ca
   10       continue
c        end if
c
      end if
c
      if (abs(real(dthty)) .gt. 0.) then
c
c        rotate qb about an axis parallel to the y-axis
c
         ca = cos(dthty)
         sa = sin(dthty)
c
c        if (ibpntsg(nbl,1).gt.0) then
c           lend = lsta+ibpntsg(nbl,1)-1
            do 20 l=lsta,lend
c           qb2 = qb(ibcg(l),2,iset)
c           qb(ibcg(l),2,iset) =  qb2*ca + qb(ibcg(l),4,iset)*sa
c           qb(ibcg(l),4,iset) = -qb2*sa + qb(ibcg(l),4,iset)*ca
            qb2 = qb(l,2,iset)
            qb(l,2,iset) =  qb2*ca + qb(l,4,iset)*sa
            qb(l,4,iset) = -qb2*sa + qb(l,4,iset)*ca
   20       continue
c        end if
c
      end if
c
      if (abs(real(dthtz)) .gt. 0.) then
c
c        rotate qb about an axis parallel to the z-axis
c
         ca = cos(dthtz)
         sa = sin(dthtz)
c
c        if (ibpntsg(nbl,1).gt.0) then
c           lend = lsta+ibpntsg(nbl,1)-1
            do 30 l=lsta,lend
c           qb2 = qb(ibcg(l),2,iset)
c           qb(ibcg(l),2,iset) = qb2*ca - qb(ibcg(l),3,iset)*sa
c           qb(ibcg(l),3,iset) = qb2*sa + qb(ibcg(l),3,iset)*ca
            qb2 = qb(l,2,iset)
            qb(l,2,iset) = qb2*ca - qb(l,3,iset)*sa
            qb(l,3,iset) = qb2*sa + qb(l,3,iset)*ca
   30       continue
c        end if
c
      end if
c
      return
      end
