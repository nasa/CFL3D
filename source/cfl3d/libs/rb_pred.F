      subroutine rb_pred(zrg1,zrg0,alf1,alf0,thtd1,thtd0,dzdt) 
c***********************************************************************
c     Purpose: Update the rigid body displacements and velocities via
c              a predictor step of the rigid body equations of motion.
c
c              
c        Reference: Vinh, Lam-Son, Edwards, J.W., Seidel, D. A., Batina,
c                   J. T., "Transonic Stability and Control of Aircraft
c                   Using CFD Methods," AIAA Paper 88-4374.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      common /rbstmt1/ bmatrb(4,4),gforcnrb(4),gforcnmrb(4),gforcsrb(4),
     .                 stmrb(4,4),stmirb(4,4),xsrb(4),xxnrb(4),xnmrb(4),
     .                 x0rb(4)
      common /rbstmt2/ tmass,yinert,uinfrb,qinfrb,greflrb,gaccel,crefrb,
     .                 xtmref,areat 
c
c     rigid body displacement and velocity prediction via rigid body
c     equations of motion
c
      do j = 1,4
         gforcnmrb(j) = gforcnrb(j)
      end do
      do j = 1,4
         gforcnrb(j) = gforcsrb(j)
      end do
c
      zrg0 = -xnmrb(1)/greflrb 
      alf0 =  xnmrb(3)
      thtd0=  xnmrb(4)
      do n=1,4
         xxnrb(n) = xsrb(n)
      end do
c
      do n=1,4
         xsrb(n) = 0.
         do j=1,4
            xsrb(n) = xsrb(n) + stmrb(n,j)*xxnrb(j)
     .              + .5*stmirb(n,j)*(3.*gforcnrb(j) - gforcnmrb(j))
         end do
      end do
c
      zrg1 = -xsrb(1)/greflrb 
      dzdt = -xsrb(2)/greflrb
      alf1 =  xsrb(3)
      thtd1=  xsrb(4)
c
      return
      end
