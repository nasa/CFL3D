      subroutine setqc0(jdim,kdim,idim,q,qc0)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Store conservative variables for use in 2nd
c     order temporal differencing and subiteration
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),qc0(jdim,kdim,idim-1,5)
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
c 
      idim1 = idim-1
      nt    = jdim*kdim
      nplq  = min(idim1,999000/nt)
      npl   = nplq
      do 20 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      n = nt*npl - jdim -1
c
cdir$ ivdep 
      do 10 izz=1,n 
c
      qc0(izz,1,i,1) = q(izz,1,i,1)
c
      qc0(izz,1,i,2) = q(izz,1,i,1)*q(izz,1,i,2)
c
      qc0(izz,1,i,3) = q(izz,1,i,1)*q(izz,1,i,3)
c
      qc0(izz,1,i,4) = q(izz,1,i,1)*q(izz,1,i,4)
c
      qc0(izz,1,i,5) = q(izz,1,i,5)/gm1+0.5*q(izz,1,i,1)
     . *(q(izz,1,i,2)**2+q(izz,1,i,3)**2+q(izz,1,i,4)**2)
   10 continue
   20 continue    
      return
      end
