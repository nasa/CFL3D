      subroutine movseg(iseg,is,it,mbloc,msegt,mtot,
     .                  imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                  nxtseg,intrfc,ipatch,nsubbl,idobl,
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
c
      nseg(is) = nseg(is) - 1
      nseg(it) = nseg(it) + 1
      if (iseg.eq.idseg(is)) then
         idseg(is) = nxtseg(iseg)
      else
         itmp = iprvsg(iseg,idseg(is),mbloc,msegt,mtot,
     .                        imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                        nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                        nseg,idno,ijk,idseg,idnext)
         nxtseg(itmp) = nxtseg(iseg)
      endif
      nxtseg(iseg) = idseg(it)
      idseg(it) = iseg
      idbloc(iseg) = it
      itmp = intrfc(iseg)
      if (itmp.ne.0) then
         imap(7,itmp) = it
      endif
c
      return
      end
