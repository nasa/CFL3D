      subroutine mkseg(iseg,it,ns,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
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
c
      common /segment/ nsgtop
c
      nsgtop = nsgtop + 1
      ns = nsgtop
      idbloc(ns) = it
      nseg(it) = nseg(it) + 1
      nxtseg(ns) = idseg(it)
      idseg(it) = ns
      intrfc(ns) = intrfc(iseg)
      do 10 i=1,msegt
         imap(i,ns) = imap(i,iseg)
 10   continue
      do 11 i=1,4
         iold(i,ns) = iold(i,iseg)
 11   continue
      do 12 i=1,msegt
         xdum(i,ns) = xdum(i,iseg)
 12   continue
      val(ns) = val(iseg)
c
      return
      end
