      subroutine shortinp(nb1,ne1,nb2,ne2,idm,jdm,kdm,iface)
c
c     $Id$
c
c     converts cfl3d "shortcut" input to full range of index values
c
      if (nb1.eq.0) nb1 = 1
      if (nb2.eq.0) nb2 = 1
      if (ne1.eq.0) then
         if (iface.eq.1 .or. iface.eq.2) ne1 = jdm
         if (iface.eq.3 .or. iface.eq.4) ne1 = kdm
         if (iface.eq.5 .or. iface.eq.6) ne1 = idm
      end if
      if (ne2.eq.0) then
         if (iface.eq.1 .or. iface.eq.2) ne2 = kdm
         if (iface.eq.3 .or. iface.eq.4) ne2 = idm
         if (iface.eq.5 .or. iface.eq.6) ne2 = jdm
      end if
c 
      return
      end 
