      subroutine ipadd(ibl,iifc,mbloc,msegt,mtot,
     .                 imap,idbloc,ivisb,itrb,val,xdum,iold,
     .                 nxtseg,intrfc,ipatch,nsubbl,idobl,
     .                 nseg,idno,ijk,idseg,idnext)
c
      dimension nsubbl(mbloc),idobl(mbloc),nseg(mbloc),idno(mbloc),
     .          ijk(6,mbloc),idseg(mbloc),idnext(mbloc)
      dimension imap(msegt,mtot),idbloc(mtot),
     .          ivisb(msegt,mbloc),itrb(7,mbloc),
     .          val(mtot),xdum(msegt,mtot),
     .          iold(4,mtot),nxtseg(mtot),intrfc(mtot),
     .          ipatch(mtot)
c
      ifc = abs(iifc)
      itop  = ipatch(1)
      i = 2
 1    if (i.ge.itop) goto 2
      j = i + 1
      if (ipatch(i).eq.ibl .and. ipatch(j).eq.ifc) return
      i = j + 1
      goto 1
 2    continue
      
      ipatch(itop) = ibl
      itop = itop + 1
      ipatch(itop) = ifc
      itop = itop + 1
      ipatch(1) = itop
c      write(6,*) 'adding block ',ibl,', face ',ifc
      
      return
      end
