      subroutine gradinfo(delh,iunit)
c***********************************************************************
c     Purpose:  Set derivative step size delh based on which of the 
c               available imaginary components is non zero; will also
c               print out a message to unit "iunit" as to which
c               derivative is active if iunit > 0.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /complx/ xmach_img,alpha_img,beta_img,reue_img,tinf_img,
     .                geom_img,surf_img
c
      if (real(xmach_img) .gt. 0.) then
         delh = xmach_img
         if (iunit .gt. 0) then
            write(iunit,'(''                    derivatives with '',
     .                 ''respect to Mach number'')')
         end if
      else if (real(alpha_img) .gt. 0.) then
         delh = alpha_img
         if (iunit .gt. 0) then
            write(iunit,'(''                       derivatives with '',
     .                 ''respect to alpha'')')
         end if
      else if (real(beta_img) .gt. 0.) then
         delh = beta_img
         if (iunit .gt. 0) then
            write(iunit,'(''                       derivatives with '',
     .                 ''respect to beta'')')
         end if
      else if (real(reue_img) .gt. 0.) then
         delh = reue_img
         if (iunit .gt. 0) then
            write(iunit,'(''                       derivatives with '',
     .                 ''respect to reue'')')
         end if
      else if (real(tinf_img) .gt. 0.) then
         delh = tinf_img
         if (iunit .gt. 0) then
            write(iunit,'(''                       derivatives with '',
     .                 ''respect to Tinf'')')
         end if
      else if (real(geom_img) .gt. 0.) then
         delh = real(geom_img)
         if (iunit .gt. 0) then
            write(iunit,'(''                      derivatives with '',
     .                 ''respect to geometry'')')
         end if
      else if (real(surf_img) .gt. 0.) then
         delh = real(surf_img)
         if (iunit .gt. 0) then
            write(iunit,'(''                      derivatives with '',
     .                 ''respect to (surface) geometry'')')
         end if
      end if
c
      return
      end
