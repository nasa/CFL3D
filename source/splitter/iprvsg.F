      integer function iprvsg(iseg,istart,mbloc,msegt,mtot,
     .                        imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                        nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                        nseg,idno,ijk,idseg,idnext)
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot)
c
      i = istart
 100  if (nxtseg(i).eq.iseg) goto 999
      i = nxtseg(i)
      goto 100
 999  continue
      iprvsg = i
c
      return
      end
