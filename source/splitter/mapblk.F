      subroutine mapblk(nbl,id,imin,imax,jmin,jmax,kmin,kmax,mbloc,
     .                  msegt,mtot,imap,idbloc,ivisb,itrb,val,xdum,
     .                  iold,nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                  nseg,idno,ijk,idseg,idnext)
c
c     $Id$
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot)
      dimension imin(*),imax(*),jmin(*),jmax(*),kmin(*),kmax(*)
      dimension id(*)
c
      common /oldbl/ noldbl
c
      nbl = 0
      do 100 n=1,noldbl
         ibl = idobl(n)
         do 101 n2=1,nsubbl(n)
            nbl = nbl + 1
            id(nbl) = n
            imin(nbl) = ijk(1,ibl)
            imax(nbl) = ijk(2,ibl)
            jmin(nbl) = ijk(3,ibl)
            jmax(nbl) = ijk(4,ibl)
            kmin(nbl) = ijk(5,ibl)
            kmax(nbl) = ijk(6,ibl)
            ibl = idnext(ibl)
 101     continue
 100  continue
c
      return
      end
