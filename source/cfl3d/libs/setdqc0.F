      subroutine setdqc0(jdim,kdim,idim,q,qc0,dqc0)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Store conservative variables: Q(n)-Q(n-1)
c     for use in 2nd order temporal differencing
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q(jdim,kdim,idim,5),qc0(jdim,kdim,idim-1,5),
     .          dqc0(jdim,kdim,idim-1,5)
      common /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
c 
      idim1 = idim-1
      nt    = jdim*kdim
      nplq  = min(idim1,999000/nt)
      npl   = nplq
c
      do 50 i=1,idim1,nplq
      if (i+npl-1.gt.idim1) npl = idim1-i+1
      n = nt*npl - jdim -1
c
c     if rho at 1st point is zero, then we know that we do not have qc0
c     data from the restart; no choice but to set dqc0=0 for 1st iteration
c
      if (real(qc0(1,1,1,1)) .eq. 0.) then
         do 20 l=1,5
cdir$ ivdep 
         do 10 izz=1,n 
   10    dqc0(izz,1,i,l) = 0.0
   20    continue
c
      else
c
cdir$ ivdep 
         do 30 izz=1,n
c
c        construct and store conserved variables
c
         dqc0(izz,1,i,1) = q(izz,1,i,1) - qc0(izz,1,i,1)
         dqc0(izz,1,i,2) = q(izz,1,i,1)*q(izz,1,i,2) - qc0(izz,1,i,2)
         dqc0(izz,1,i,3) = q(izz,1,i,1)*q(izz,1,i,3) - qc0(izz,1,i,3)
         dqc0(izz,1,i,4) = q(izz,1,i,1)*q(izz,1,i,4) - qc0(izz,1,i,4)
         dqc0(izz,1,i,5) = q(izz,1,i,5)/gm1+0.5*q(izz,1,i,1)
     .    *(q(izz,1,i,2)**2+q(izz,1,i,3)**2+q(izz,1,i,4)**2)
     .    - qc0(izz,1,i,5)
   30    continue
      end if
   50 continue    
      return
      end
