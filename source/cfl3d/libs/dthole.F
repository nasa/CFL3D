      subroutine dthole(jdim,kdim,idim,dtj,vol,blank,dtmin,
     .                  nou,bou,nbuf,ibufdim)
c***********************************************************************
c     Purpose:  Update the dt values for the hole and fringe cells;
c     the values will be replaced by dtmin.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension blank(jdim,kdim,idim)
      dimension dtj(jdim,kdim,idim-1),vol(jdim,kdim,idim-1) 
c
c     nou(1) = min(nou(1)+1,ibufdim)
c     if (isklton.gt.0) then
c        nou(1) = min(nou(1)+1,ibufdim)
c        write(bou(nou(1),1),*)' in dthole, dtmin= ',real(dtmin)
c     end if
c
      kn = jdim*kdim*(idim-1)
cdir$ ivdep
      do 12 n=1,kn
      dtj(n,1,1) = ccvmgt(vol(n,1,1)/dtmin,dtj(n,1,1),
     .                  (real(blank(n,1,1)).eq.0.e0))
   12 continue
      return
      end 
