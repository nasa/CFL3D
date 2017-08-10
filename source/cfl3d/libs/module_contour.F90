!  ---------------------------------------------------------------------------
!  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
!  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
!  with point-matched, patched, overset, or embedded connectivities. Both
!  multigrid and mesh sequencing are available in time-accurate or
!  steady-state modes.
!
!  Copyright 2001 United States Government as represented by the Administrator
!  of the National Aeronautics and Space Administration. All Rights Reserved.
! 
!  The CFL3D platform is licensed under the Apache License, Version 2.0 
!  (the "License"); you may not use this file except in compliance with the 
!  License. You may obtain a copy of the License at 
!  http://www.apache.org/licenses/LICENSE-2.0. 
! 
!  Unless required by applicable law or agreed to in writing, software 
!  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
!  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
!  License for the specific language governing permissions and limitations 
!  under the License.
!  ---------------------------------------------------------------------------
!
MODULE module_contour
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: cont_plot
CONTAINS
  SUBROUTINE  cont_plot(&
                        nframe,jdim,kdim,idim,nummem,nblk,q,x,y,z,vist3d,ux,turre, smin)
    INTEGER,INTENT(in)::jdim,kdim,idim,nummem, nblk,nframe
    REAL,INTENT(in) :: x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),&
         q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim),&
         ux(jdim-1,kdim-1,idim-1,9),turre(jdim,kdim,idim,nummem),&
         smin(jdim-1,kdim-1,idim-1)

    ! common block
    COMMON /twod/ i2d
    INTEGER ::i2d
    COMMON /fluid2/ pr,prt,cbar
    REAL :: pr,prt,cbar
    COMMON /reyue/ reue,tinf,ivisc(3)
    REAL::reue, tinf
    INTEGER ::ivisc
    COMMON /fluid/gamma,gm1,gp1,gm1g,gp1g,ggm1
    REAL::gamma,gm1,gp1,gm1g,gp1g,ggm1
    COMMON /info/ title(20),rkap(3),xmach
    REAL :: title,rkap,xmach

    ! local variable
    CHARACTER(len=20) :: str,strframe
    character(len=120)::plt_title
    INTEGER :: i,j,k,m
    REAL :: xc,yc,zc
    REAL :: c2b,c2bp, fnu(jdim),tt,re,utau(jdim),vel&
         ,tau_w,fnuloc,Rel(jdim),cf(jdim)
    REAL :: deltabl, coef
    
    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10,vel1

    IF(nummem/=2.and.nummem/=7) RETURN
    WRITE(str,'(I5.5)')nblk
    WRITE(strframe,'(I4.4)')nframe
    OPEN(1234,file="turre-"//TRIM(ADJUSTL(strframe))//"-blk"//TRIM(ADJUSTL(str))//".plt")

    WRITE(plt_title,'(2(1x,A,es12.5))')'Re=',reue,'Mach=',xmach
    WRITE(1234,'(A)')'title="'//TRIM(ADJUSTL(plt_title))//'"'
    IF(nummem==7) THEN
       WRITE(1234,'(A)')'variables="x","y","z","dist","q1", "q2", "q3", "q4", "q5", '//&
            '"dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz", '//&
            '"mu", "mut", "t11", "t22", "t33", "t12", "t23", "t13", "omega/eps/zeta"'
    ELSEIF(nummem==2) THEN
       WRITE(1234,'(A)')'variables="x","y","z","dist","q1", "q2", "q3", "q4", "q5", '//&
            '"dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz", '//&
            ',"mu","mut", "omg/eps/zeta", "tke"'
    ENDIF
    IF(i2d==1) THEN
       WRITE(1234,'(A,i5,A,i5)')"ZONE  I=", jdim-1," J=", kdim-1
    ELSE
       WRITE(1234,'(A,i5,A,i5,A,i5)')"ZONE  I=", idim-1, " J=",jdim-1," K=", kdim-1
    ENDIF

    DO k=1,kdim-1
       DO j=1,jdim-1
          DO i=1,idim-1
             xc = 0.125*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j,k,i+1) + &
                  x(j+1,k+1,i) + x(j+1,k,i+1) + x(j,k+1,i+1)+x(j+1,k+1,i+1))
             yc = 0.125*(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+y(j,k,i+1) + &
                  y(j+1,k+1,i) + y(j+1,k,i+1) + y(j,k+1,i+1)+y(j+1,k+1,i+1))
             zc = 0.125*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j,k,i+1) + &
                  z(j+1,k+1,i) + z(j+1,k,i+1) + z(j,k+1,i+1)+z(j+1,k+1,i+1))
             c2b=cbar/tinf
             c2bp=c2b+1.0
             tt=gamma*q(j,k,i,5)/q(j,k,i,1)
             fnuloc=c2bp*tt*SQRT(tt)/(c2b+tt)
             WRITE(1234,'(25(1x,es12.5))')xc,yc,zc,smin(j,k,i),q(j,k,i,1:5), ux(j,k,i,1:9),&
                  fnuloc,vist3d(j,k,i), turre(j,k,i,1:nummem)
          ENDDO
       ENDDO
    ENDDO
    CLOSE(1234)

  END SUBROUTINE cont_plot
END MODULE module_contour
