      subroutine rotateq0(ld,md,q0,q0rot,lsta,lend,msta,mend,
     .                    dthtx,dthty,dthtz)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Rotate solution at ghost points contained in array 
c     q0 (qi0/qj0/qk0) through a specified angle;                
c     rotated solution stored in q0rot
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension q0(ld,md,5,4),q0rot(ld,md,5,4)
c
      if (abs(real(dthtx)) .gt. 0.) then
c
c        rotate q0 about an axis parallel to the x-axis
c
         ca = cos(dthtx)
         sa = sin(dthtx)
c
         do 110 n=1,4
         do 120 l=lsta,lend
         do 130 m=msta,mend
         q0rot(l,m,1,n) =  q0(l,m,1,n)
         q0rot(l,m,2,n) =  q0(l,m,2,n)
         temp           =  q0(l,m,3,n)
         q0rot(l,m,3,n) =  q0(l,m,3,n)*ca - q0(l,m,4,n)*sa
         q0rot(l,m,4,n) =  temp       *sa + q0(l,m,4,n)*ca
         q0rot(l,m,5,n) =  q0(l,m,5,n)
 130     continue
 120     continue
 110     continue
c
      else if (abs(real(dthty)) .gt. 0.) then
c
c        rotate q0 about an axis parallel to the y-axis
c
         ca = cos(dthty)
         sa = sin(dthty)
c
         do 210 n=1,4
         do 220 l=lsta,lend
         do 230 m=msta,mend
         q0rot(l,m,1,n) =  q0(l,m,1,n)
         temp           =  q0(l,m,2,n)
         q0rot(l,m,2,n) =  q0(l,m,2,n)*ca + q0(l,m,4,n)*sa
         q0rot(l,m,3,n) =  q0(l,m,3,n)
         q0rot(l,m,4,n) = -temp       *sa + q0(l,m,4,n)*ca
         q0rot(l,m,5,n) =  q0(l,m,5,n)
 230     continue
 220     continue
 210     continue
c
      else if (abs(real(dthtz)) .gt. 0.) then
c
c        rotate q0 about an axis parallel to the z-axis
c
         ca = cos(dthtz)
         sa = sin(dthtz)
c
         do 310 n=1,4
         do 320 l=lsta,lend
         do 330 m=msta,mend
         q0rot(l,m,1,n) =  q0(l,m,1,n)
         temp           =  q0(l,m,2,n)
         q0rot(l,m,2,n) =  q0(l,m,2,n)*ca - q0(l,m,3,n)*sa
         q0rot(l,m,3,n) =  temp       *sa + q0(l,m,3,n)*ca
         q0rot(l,m,4,n) =  q0(l,m,4,n)
         q0rot(l,m,5,n) =  q0(l,m,5,n)
 330     continue
 320     continue
 310     continue
c
      else
         do 410 n=1,4
         do 420 l=lsta,lend
         do 430 m=msta,mend
         q0rot(l,m,1,n) =  q0(l,m,1,n)
         q0rot(l,m,2,n) =  q0(l,m,2,n)
         q0rot(l,m,3,n) =  q0(l,m,3,n)
         q0rot(l,m,4,n) =  q0(l,m,4,n)
         q0rot(l,m,5,n) =  q0(l,m,5,n)
 430     continue
 420     continue
 410     continue
c
      end if
c
      return
      end
