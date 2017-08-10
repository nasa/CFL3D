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
!
! An output module for profiles at specified locations.
!
! Author: Xudong Xiao
! funded under NASA Cooperative Agreement NNX11AI56A
! Date: June-11-2012
!
!
MODULE module_profileout
  IMPLICIT NONE
  
  TYPE output_parm
     character(len=40) :: name
     INTEGER :: iblk
     INTEGER :: istr,iend,jstr,jend,kstr,kend
     INTEGER :: wall_str,wall_end,wall_dir
  END TYPE output_parm
 
  TYPE(output_parm),ALLOCATABLE :: oparm(:)
  TYPE(output_parm),ALLOCATABLE :: cfparm(:) ! for Cf cp output
  PRIVATE :: set_default
CONTAINS
  SUBROUTINE profile_read_input(filename)
    CHARACTER(len=*),INTENT(in) :: filename
    
    LOGICAL :: lexist
    INTEGER :: iost,iost1
    CHARACTER(len=170) :: cline,str,key
    REAL  :: rval
    INTEGER :: ncp_profs,nprofiles
    
    INQUIRE(file=filename, exist=lexist)
    IF(.NOT.lexist) RETURN
    
    ncp_profs = 0
    nprofiles = 0
    
    DO WHILE(.TRUE.)
       OPEN(1235,file=filename, status='old')
       READ(1235,'(A)',iostat=iost)cline
       IF(iost/=0) EXIT
       IF(LEN_TRIM(cline)==0) CYCLE
       str = ADJUSTL(cline)
       IF(str(1:1)=='!') CYCLE
       READ(cline,*,iostat=iost1)key,rval
       CALL toupper(key)
       IF(TRIM(ADJUSTL(key))=="PROFILE") THEN
          nprofiles= nprofiles + 1
       ENDIF
       IF(TRIM(ADJUSTL(key))=="CP_PROFILE") THEN
          ncp_profs = ncp_profs+1
       ENDIF
    ENDDO
    IF(nprofiles==0.AND.ncp_profs==0) THEN
       CLOSE(1235)
       RETURN
    ENDIF
    
    IF(nprofiles>0) ALLOCATE(oparm(nprofiles))
    IF(ncp_profs>0) ALLOCATE(cfparm(ncp_profs))
    
    REWIND(1235)
    nprofiles = 0
    ncp_profs = 0
    DO WHILE(.TRUE.)
       OPEN(1235,file=filename, status='old')
       READ(1235,'(A)',iostat=iost)cline
       IF(iost/=0) EXIT
       IF(LEN_TRIM(cline)==0) CYCLE
       str = ADJUSTL(cline)
       IF(str(1:1)=='!') CYCLE
       READ(cline,*,iostat=iost1)key,rval
       CALL toupper(key)
       SELECT CASE(TRIM(ADJUSTL(key)))
       CASE("PROFILE")
          !WRITE(*,*)TRIM(cline)
          nprofiles = nprofiles + 1
          READ(cline,*)key,oparm(nprofiles)%name, &
               oparm(nprofiles)%iblk,&
               oparm(nprofiles)%istr,&
               oparm(nprofiles)%iend,&
               oparm(nprofiles)%jstr,&
               oparm(nprofiles)%jend,&
               oparm(nprofiles)%kstr,&
               oparm(nprofiles)%kend
       CASE("CP_PROFILE")
          !WRITE(*,*)TRIM(cline)
          ncp_profs = ncp_profs + 1
          READ(cline,*)key,cfparm(ncp_profs)%name, &
               cfparm(ncp_profs)%iblk,&
               cfparm(ncp_profs)%istr,&
               cfparm(ncp_profs)%iend,&
               cfparm(ncp_profs)%jstr,&
               cfparm(ncp_profs)%jend,&
               cfparm(ncp_profs)%kstr,&
               cfparm(ncp_profs)%kend,&
               cfparm(ncp_profs)%wall_dir,&
               cfparm(ncp_profs)%wall_str,&
               cfparm(ncp_profs)%wall_end
       END SELECT
    ENDDO
    CLOSE(1235)
    RETURN
  CONTAINS
    SUBROUTINE toupper(str)
      CHARACTER(len=*)::str
      INTEGER :: ibig,ismall
      INTEGER :: i
      ismall = ichar('a')
      ibig   = ichar('A')
      DO i=1,LEN(str)
         IF(str(i:i)>='a' .AND. str(i:i)<='z') THEN
            str(i:i) = CHAR(ICHAR(str(i:i))-ismall + ibig)
         ENDIF
      ENDDO
    END SUBROUTINE toupper
  END SUBROUTINE profile_read_input
    
  SUBROUTINE set_default(istr,iend,imin,imax)
    INTEGER :: istr,iend
    INTEGER,intent(in) ::imin,imax

    IF(istr==0) THEN 
       istr = imin
    ENDIF
    IF(iend==0) THEN
       iend = imax
    ENDIF
    
  END SUBROUTINE set_default
    

  SUBROUTINE profile_plot(jdim,kdim,idim,nummem,nblk,q,x,y,z,vist3d,ux,turre,smin)
    INTEGER,INTENT(in)::jdim,kdim,idim,nummem, nblk
    REAL,INTENT(in) :: x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),&
         q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim),&
         ux(jdim,kdim,idim,9),turre(jdim,kdim,idim,nummem),&
         smin(jdim-1,kdim-1,idim-1)
    
    INTEGER :: i,j,k,m,istr,jstr,kstr,iend,jend,kend,n
    REAL ::xc,yc,zc
    IF(.NOT.ALLOCATED(oparm)) RETURN
    DO n=1,SIZE(oparm,1)
       IF(oparm(n)%iblk/=nblk) CYCLE
       OPEN(1234,file=TRIM(ADJUSTL(oparm(n)%name))//".plt")
       IF(nummem==7) THEN
          WRITE(1234,'(A)')'variables="x","y","z","dist","q1", "q2", "q3", "q4", "q5", '//&
               '"dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz", '//&
               '"mut", "t11", "t22", "t33", "t12", "t23", "t13", "zeta/omega/eps"'
       ELSEIF(nummem==2) THEN
          WRITE(1234,'(A)')'variables="x","y","z","dist","q1", "q2", "q3", "q4", "q5", '//&
            '"dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz", '//&
            ',"mut", "omega/eps/zeta", "tke"'
       ELSE
          CLOSE(1234,status='delete')
          WRITE(*,*)"profile_plot: unknown nummem",oparm(n)%name,nummem
          CYCLE
       ENDIF
       
       WRITE(1234,'(A)')'TITLE="'//TRIM(ADJUSTL(oparm(n)%name))//'"'

       istr =oparm(n)%istr
       jstr =oparm(n)%jstr
       kstr =oparm(n)%kstr

       iend =oparm(n)%iend
       jend =oparm(n)%jend
       kend =oparm(n)%kend

       CALL set_default(istr,iend,1,idim-1)
       CALL set_default(jstr,jend,1,jdim-1)
       CALL set_default(kstr,kend,1,kdim-1)
       
       WRITE(1234,'(A)')'ZONE T="'//TRIM(ADJUSTL(oparm(n)%name))
 !            '" I=', iend-istr+1, " J=",jend-jstr+1," K=", kend-kstr+1
       DO k=kstr,kend
          DO j=jstr,jend
             DO i=istr,iend
                
                xc = 0.125*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j,k,i+1) + &
                     x(j+1,k+1,i) + x(j+1,k,i+1) + x(j,k+1,i+1)+x(j+1,k+1,i+1))
                yc = 0.125*(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+y(j,k,i+1) + &
                     y(j+1,k+1,i) + y(j+1,k,i+1) + y(j,k+1,i+1)+y(j+1,k+1,i+1))
                zc = 0.125*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j,k,i+1) + &
                     z(j+1,k+1,i) + z(j+1,k,i+1) + z(j,k+1,i+1)+z(j+1,k+1,i+1))
                WRITE(1234,'(35(1x,es12.5))')xc,yc,zc,smin(j,k,i),q(j,k,i,1:5), ux(j,k,i,1:9),&
                     vist3d(j,k,i), turre(j,k,i,1:nummem)
             ENDDO
          ENDDO
       ENDDO
       CLOSE(1234)
    END DO
    
  END SUBROUTINE profile_plot

  
  SUBROUTINE cfcp_plot(jdim,kdim,idim,nummem,nblk,q,x,y,z,vist3d,ux,turre,smin)
    INTEGER,INTENT(in)::jdim,kdim,idim,nummem, nblk
    REAL,INTENT(in) :: x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),&
         q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim),&
         ux(jdim,kdim,idim,9),turre(jdim,kdim,idim,nummem),&
         smin(jdim-1,kdim-1,idim-1)
    
    INTEGER :: i,j,k,m,istr,jstr,kstr,iend,jend,kend,n
    REAL ::xc,yc,zc
    
    COMMON /fluid2/ pr,prt,cbar
    REAL :: pr,prt,cbar
    COMMON /reyue/ reue,tinf,ivisc(3)
    REAL::reue, tinf
    INTEGER ::ivisc
    COMMON /fluid/gamma,gm1,gp1,gm1g,gp1g,ggm1
    REAL::gamma,gm1,gp1,gm1g,gp1g,ggm1
    COMMON /info/ title(20),rkap(3),xmach
    REAL :: title,rkap,xmach
    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10,vel1


    REAL :: c2b, c2bp,tt,fnu,vel,tau_w,re
    INTEGER :: idx(3),k1
    REAL :: dprev,uprev,ruprev,rho,rhou_ave,u_ave,theta,dstar,uinf
    REAL, parameter :: rinf=1.0

    


    Re = Reue/xmach
    
    IF(.NOT.ALLOCATED(cfparm)) RETURN

    DO n=1,SIZE(cfparm,1)
       IF(cfparm(n)%iblk/=nblk) CYCLE
       istr =cfparm(n)%istr
       jstr =cfparm(n)%jstr
       kstr =cfparm(n)%kstr
       iend =cfparm(n)%iend
       jend =cfparm(n)%jend
       kend =cfparm(n)%kend
       CALL set_default(istr,iend,1,idim-1)
       CALL set_default(jstr,jend,1,jdim-1)
       CALL set_default(kstr,kend,1,kdim-1)
       select case(cfparm(n)%wall_dir)
       case(1)
         call set_default(cfparm(n)%wall_str,cfparm(n)%wall_end,1,idim-1)
       case(2)
         call set_default(cfparm(n)%wall_str,cfparm(n)%wall_end,1,jdim-1)
       case(3)
         call set_default(cfparm(n)%wall_str,cfparm(n)%wall_end,1,kdim-1)
       end select

       OPEN(1234,file=TRIM(ADJUSTL(cfparm(n)%name))//".plt")
       WRITE(1234,'(A)')'variables=x,y,z,cf,cp,theta,Re_theta,dstar,Re_dstar'
       WRITE(1234,'(A)')'ZONE T="'//TRIM(ADJUSTL(cfparm(n)%name))//'"' !//&
!            '" I=', iend-istr+1, " J=",jend-jstr+1," K=", kend-kstr+1
       
       DO i=istr,iend
          DO k=kstr,kend
             DO j=jstr,jend
                c2b=cbar/tinf
                c2bp=c2b+1.0
                tt=gamma*q(j,k,i,5)/q(j,k,i,1)
                fnu=c2bp*tt*SQRT(tt)/(c2b+tt)
                
                vel = SIGN(SQRT(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2),q(j,k,i,2))
                tau_w = fnu*ABS(vel)/ABS(smin(j,k,i))/Re
                !utau = SQRT(tau_w/q(j,1,i,1))
                xc = 0.125*(x(j,k,i)+x(j+1,k,i)+x(j,k+1,i)+x(j,k,i+1) + &
                     x(j+1,k+1,i) + x(j+1,k,i+1) + x(j,k+1,i+1)+x(j+1,k+1,i+1))
                yc = 0.125*(y(j,k,i)+y(j+1,k,i)+y(j,k+1,i)+y(j,k,i+1) + &
                     y(j+1,k+1,i) + y(j+1,k,i+1) + y(j,k+1,i+1)+y(j+1,k+1,i+1))
                zc = 0.125*(z(j,k,i)+z(j+1,k,i)+z(j,k+1,i)+z(j,k,i+1) + &
                     z(j+1,k+1,i) + z(j+1,k,i+1) + z(j,k+1,i+1)+z(j+1,k+1,i+1))
                idx=(/i,j,k/)
                dprev = 0
                uprev = 0
                ruprev = 0
                theta = 0
                dstar = 0
                idx(cfparm(n)%wall_dir) = cfparm(n)%wall_end
                uinf = 0.
                ! pick the maximum velocity as the edge velocity
                DO k1=cfparm(n)%wall_str,cfparm(n)%wall_end,&
                     SIGN(1,cfparm(n)%wall_end - cfparm(n)%wall_str)
                   idx(cfparm(n)%wall_dir) = k1
                   uinf = q(idx(2),idx(3),idx(1),2)
                   vel = SQRT(q(idx(2),idx(3),idx(1),2)**2+&
                                   q(idx(2),idx(3),idx(1),3)**2+&
                                   q(idx(2),idx(3),idx(1),4)**2)
                   uinf = max(uinf,vel) 
                enddo
                DO k1=cfparm(n)%wall_str,cfparm(n)%wall_end,&
                     SIGN(1,cfparm(n)%wall_end - cfparm(n)%wall_str)
                   idx(cfparm(n)%wall_dir) = k1
                  
                   vel = SQRT(q(idx(2),idx(3),idx(1),2)**2+&
                                   q(idx(2),idx(3),idx(1),3)**2+&
                                   q(idx(2),idx(3),idx(1),4)**2)
                   rho = q(idx(2),idx(3),idx(1),1)
                   rhou_ave = 0.5*(rho*vel + ruprev)
                   u_ave = 0.5*(vel+uprev)

                   theta = theta + rhou_ave/(rinf*uinf)*MAX((1.-u_ave/uinf),0.)*(smin(idx(2),idx(3),idx(1))-dprev)
                   dstar = dstar + MAX(1.-rhou_ave/(rinf*uinf),0.)*(smin(idx(2),idx(3),idx(1))-dprev)
                   dprev = smin(idx(2),idx(3),idx(1))
                   uprev = vel
                   ruprev = rho*vel
                   WRITE(*,*)k1,theta,dstar,vel,xmach
                ENDDO
                
                !Rel  = (0.5*(x(j,1,i)+x(j+1,1,i)) - x(jstart,1,i))*Reue
                WRITE(1234,'(14(1x,es12.5))') xc,yc,zc,tau_w/(0.5*rho0*xmach**2),&
                     (q(j,k,i,5)-p0)/(0.5*rho0*xmach**2),theta, reue*theta,dstar,dstar*reue
             ENDDO
          ENDDO
       ENDDO
       CLOSE(1234)
    ENDDO

  END SUBROUTINE cfcp_plot
END MODULE module_profileout
  
