!     
!     The k-omega Stress Transport Model Module
!
!     Author:  Xudong Xiao, Corvid Technologies, Inc
!
!  Reference:
!    1) Wilcox, "Turbulence Modeling for CFD", 3rd edition, pp.332-333
MODULE module_kwstm
  IMPLICIT NONE
  PRIVATE

  ! total of 8 model constants, c.f. 
  REAL,SAVE :: &  !     
       c1=9./5., &       !1
       c2=10./19., &     !2
       alpha=0.52, &  !3
       beta_0 = 0.0708,& !4
       beta_star = 0.09, &!5
       sigma=0.5, &      !6
       sigma_star = 0.6, & ! 7 
       sigma_d0 = 0.125    ! 8
  
  REAL :: alpha_hat, beta_hat, gamma_hat

  INTEGER,save :: nsubit = 1
  REAL,SAVE :: factor=1.0,factor1,factor2
  INTEGER, save :: irealizability=1
  REAL,save :: pdratio=10.0
  integer :: nfreq = 100000
  LOGICAL,SAVE :: need_contplot=.false., need_kwsplot=.false.

  PUBLIC:: kws_init,kws_main,kws_get_nummem, &
       kws_initvist,kws_read_parm,kws_plot,&
       kws_dump_movie,&
       kws_init_from_2eq,kws_save_2eq,need_contplot
       
  
  INTEGER,SAVE :: ux_signal = 0

  INTEGER ::myrank, myroot
  INTEGER :: mylevel, myicyc
  INTEGER,SAVE :: level_o=-1,icyc_o =1,debug_j=-11
  INTEGER,SAVE,public :: jcut = 64 ,jstart = 17

  integer,save :: kdiff= 0
  integer,save :: iopt=0  !
   
  ! source term for k-equation
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: src_k

  ! For Local Dual Time Stepping Method
  REAL,PUBLIC, SAVE :: cfl_psd=1.0, cfl_loc =10.0!
  !       cfl number for the pseudo time step, and local time step

#define CHECK_SOURCE_TERM
#ifdef CHECK_SOURCE_TERM
  REAL, ALLOCATABLE :: source_items(:,:,:,:,:)
#endif
  REAL, ALLOCATABLE:: source(:,:,:,:)
  REAL, ALLOCATABLE:: rhs(:,:,:,:)
  REAL, PUBLIC, SAVE :: xmutratio=0.0   ! used to add a small portion of eddy viscosity to the laminar viscosity, it is used in the [fgh]fluxv.F
  REAL, save,public:: ildts=0

  LOGICAL, save :: lkzstm_from_2eq = .false.
CONTAINS
!   The following subroutine - NO LONGER USED!!!!
! SUBROUTINE kws_set_tu(tu_inf_percent)
!   REAL, INTENT(in) :: tu_inf_percent
!   IF(tu_inf_percent>0) THEN
!      tu = tu_inf_percent/100.
!   ENDIF
! END SUBROUTINE kws_set_tu
  
!   The following subroutine - NO LONGER USED!!!!
! SUBROUTINE kws_set_edinf(eddy_visc_inf)
!   REAL, intent(in) :: eddy_visc_inf
!   
!   IF(eddy_visc_inf>0) THEN
!      xnu_ratio = eddy_visc_inf
!   ENDIF
! END SUBROUTINE kws_set_edinf
  
  !
  ! return the number of history variables needed for
  ! the stress model
  !
  FUNCTION kws_get_nummem()RESULT(iret)
    INTEGER :: iret
    iret = 7
  END FUNCTION kws_get_nummem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize the stress components and the zeta, eddy viscosity
  ! called by: libs/initvist.F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   The following subroutine - NO LONGER USED!!!!
  SUBROUTINE kws_initvist(jdim,kdim,idim,nummem,cij, vist3d)
    INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
    REAL,INTENT(out) :: cij(jdim,kdim,idim,nummem),vist3d(jdim,kdim,idim)

    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10
    INTEGER :: i,j,k

    
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             cij(j,k,i,1:nummem) = tur10(1:nummem)
             vist3d(j,k,i) = 0.01
          ENDDO
       ENDDO
    ENDDO
    
    RETURN

  END SUBROUTINE kws_initvist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read some parameters for RSM and write out
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE kws_init(nou,bou,nbuf,ibufdim)
    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
    REAL ::  title,rkap,xmach,alpha__,beta__,dt,fmax
    INTEGER :: nit,ntt, &
         idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
         iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv
    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10
    real :: tke, tu, xnu_ratio
    INTEGER :: i,j,k,m
    COMMON /mydist2/ nnodes,myhost,myid,mycomm
    INTEGER::nnodes,myhost,myid,mycomm
    INTEGER :: nbuf, ibufdim, nou(nbuf)
    character(len=120)::bou(ibufdim,nbuf)


    INTEGER,SAVE :: ivisited=0

    if(ivisited==0) then 
      call kws_read_parm(nou,bou,nbuf,ibufdim)
      ivisited = 1
    end if

    return

  END SUBROUTINE kws_init
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     read parameters from the input file: kz_consts.dat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE kws_read_parm(nou,bou,nbuf,ibufdim)
    USE module_profileout,ONLY:profile_read_input
    !
    ! common blocks needed for this subroutine
    ! 
    COMMON /mydist2/ nnodes,myhost,myid,mycomm
    INTEGER::nnodes,myhost,myid,mycomm
    
    COMMON /maxiv/ ivmx  ! to get the turbulence model flag
    INTEGER :: ivmx

    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
    REAL:: title,rkap, xmach,alpha__, beta__, dt, fmax
    INTEGER ::nit,ntt, &
         idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
         iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv

    COMMON /reyue/ reue,tinf,ivisc(3)
    REAL::reue,tinf
    INTEGER :: ivisc

    COMMON /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze, &
         iwarneddy,itime2read,itaturb,tur1cut,tur2cut,& 
         iturbord,tur1cutlev,tur2cutlev
    REAL :: cflturb,edvislim,tur1cut,tur2cut, tur1cutlev,tur2cutlev
    INTEGER :: iturbprod,nsubturb,nfreeze, iwarneddy,itime2read,itaturb,iturbord

    INTEGER :: nbuf, ibufdim, nou(nbuf)
    character(len=120)::bou(ibufdim,nbuf)

    CHARACTER(len=70) :: cline,key,str
    REAL :: rval
    !
    ! local variables
    !
    REAL:: reOMa, maORe
    INTEGER :: iost,iost1
    LOGICAL :: lexist  ! used to test if kwstm.dat exists
    INTEGER :: ival

    INQUIRE(file="kwstm.dat",exist=lexist)
    nsubit= nsubturb
    IF(lexist) THEN       
       OPEN(1234,file="kwstm.dat")
       DO WHILE(.TRUE.) 
          READ(1234,'(A)',iostat=iost)cline
          IF(iost/=0) EXIT
          IF(LEN_TRIM(cline)==0) cycle
          str = adjustl(cline)
          if(str(1:1)=='!') cycle
          READ(cline,*,iostat=iost1)key,rval
          CALL toupper(key)
          IF(iost1/=0) THEN
             IF(myid==myhost .or. myid==1) THEN
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' Unknown Input line in kwstm.dat:'')')
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(a70)') cline
             ENDIF
          ENDIF

          SELECT CASE(TRIM(ADJUSTL(key)))
          CASE("C1")
             c1=rval
          CASE("C2")
             c2 = rval
          CASE("ALPHA")
             alpha = rval
          CASE("BETA_0")
             beta_0 = rval
          CASE("BETA_STAR")
             beta_star= rval
          CASE("SIGMA")
             SIGMA= rval
          CASE("SIGMA_STAR")
             sigma_star= rval
          CASE("SIGMA_D0")
             sigma_d0= rval
!         No longer allow tu and xnu_ratio because they are input via keyword
!         (as turbintensity_inf_percent and eddy_visc_inf (see readkey.F))
!         CASE("TU","tu")
!            tu = rval
!         CASE("XNU_RATIO")
!            xnu_ratio = rval
          case("JCUT")
             jcut = nint(rval)
          case("JSTART")
             jstart = nint(rval)
          case("KDIFF")
             KDIFF = nint(rval)
          case("IREAL")
             irealizability= nint(rval)
          case("IOPT")
             iopt= nint(rval)
          case("CFLLOC")
             cfl_loc= rval
          case("PDRATIO")
             pdratio= rval
          case("MUTRATIO")
             xmutratio= rval ! used in the [fgh]fluxv subroutine
          case("ILDTS")
             ildts= rval ! used in the [fgh]fluxv subroutine
          CASE("FACTOR")
             factor= rval 
          CASE("NSUBIT")
             nsubit= nint(rval) 
          CASE("NFREQ")
             nfreq= NINT(rval) ! used in the [fgh]fluxv subroutine
          CASE("PROFILE")
          CASE("CP_PROFILE")
          CASE("CF_YPLUS_PLOT")
             need_kwsplot = nint(rval)==1
          CASE("CONTOUR_PLOT")
             need_contplot= nint(rval)==1
          CASE default
             IF(myid==myhost .or. myid==1) THEN
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'('' ERROR: Unknown Input keyword in kwstm.dat:'')')
                nou(1) = min(nou(1)+1,ibufdim)
                write(bou(nou(1),1),'(a70)') TRIM(ADJUSTL(cline))
                call termn8(myid,-1,ibufdim,nbuf,bou,nou)
             ENDIF
          END SELECT
       ENDDO
100    CLOSE(1234)
       CALL profile_read_input("kwstm.dat")
    ENDIF
    alpha_hat = (8.+c2)/11.
    beta_hat = (8*c2 - 2.)/11.
    gamma_hat = (60*c2 - 4.)/55.
    IF(myid == myhost .or. myid==1) THEN

       open(120,file='kwstm.out',form='formatted',status='unknown')

       write(120,'('' -----------------------------'')')
       write(120,'('' List of Wilcox Stress-Omega RSM model Constants'')')
       write(120,'('' -----------------------------'')')
       call wr_parm("c1",c1)
       CALL wr_parm("c2",c2)
       call wr_parm("alpha",alpha)
       call wr_parm("beta_0", beta_0)
       call wr_parm("beta_star", beta_star)
       CALL wr_parm("sigma", sigma)
       CALL wr_parm("sigma_star",sigma_star)
       CALL wr_parm("sigma_d0", sigma_d0)

       write(120,'('' -----------------------------'')')
       write(120,'('' B.L. Profile Location'')')
       write(120,'('' -----------------------------'')')

       call wr_parm("jcut", real(jcut))
       call wr_parm("jstart", real(jstart))

       write(120,'('' -----------------------------'')')
       write(120,'('' K diffusion option 0-T_ij, 1-Robs K version(not supported), 2-(6.50)  3-(6.49)'')')
       write(120,'('' -----------------------------'')')

       call wr_parm("kdiff", real(kdiff))
       
       write(120,'('' -----------------------------'')')
       write(120,'('' Realizability  0-off  1-on'')')
       write(120,'('' -----------------------------'')')
       
       call wr_parm("ireal", real(irealizability))
       
       write(120,'('' -----------------------------'')')
       write(120,'('' Miscellaneous'')')
       write(120,'('' -----------------------------'')')
       
       CALL wr_parm("cflloc", REAL(cfl_loc))
       CALL wr_parm("pdratio", REAL(pdratio))
       CALL wr_parm("mutratio", REAL(xmutratio))
       CALL wr_parm("factor", REAL(factor))
       CALL wr_parm("nsubit", REAL(nsubit))
       CALL wr_parm("nfreq", REAL(nfreq))
       CALL wr_parm("ildts", REAL(ildts))
       ival = 0
       if(need_contplot) ival = 1
       call wr_parm("contour_plot",real(ival))

       ival = 0
       if(need_kwsplot) ival = 1         
       call wr_parm("CF_YPLUS_PLOT",real(ival))
       
       write(120,'('' -----------------------------'')')

       close (120)

    ENDIF
  CONTAINS
    SUBROUTINE wr_parm(key,rval)
      CHARACTER(len=*),INTENT(in)::key
      REAL,intent(in) :: rval
      write(120,'(a20,1x,es12.5)') key,rval
    END SUBROUTINE wr_parm
    
    subroutine toupper(str)
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
    end subroutine
  END SUBROUTINE kws_read_parm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solve seven equations for Wilcox stress-omega Reynolds Stress model
  ! (Wilcox, D. C., "Turbulence Modeling For CFD", 3rd edition, 2006)
  !
  ! 1=-u'u', 2=-v'v', 3=-w'w', 4=-u'v', 5=-v'w', 6=-u'w', 7=omega
  !
  ! Input:
  !    Turb:  the turbulence variables saved in W array, pointed by
  !           lw(19,nbl), c.f. pointers.F for details.
  ! output:
  !    turb 
  !    sumn: RMS of the residual of each turbulence quantity.
  !    negn: number of occurances of negative quantities (only applies to 
  !           a) t11=-u'u', t22=-v'v', t33=-w'w', and omega
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE kws_main(jdim,kdim,idim,myid,myhost,nummem,&
       tj0,  tk0, ti0,turb,&   ! stress components + enstropy
       ux, &    ! velocity gradients
       sj,sk,si,vol,vj0,vk0,vi0, &         ! metrics
       q,qj0,qk0,qi0,dtj,vist3d,level,icyc,sumn,negn, &
       zksav2,smin,x,y,z,nbl)

    ! --- input ---
    INTEGER,INTENT(in) :: jdim,kdim,idim,myid,myhost,nummem
    !   icyc :: sub-iteration id in the timeaccurate calculation
    REAL,intent(inout) :: tj0(kdim,idim-1,nummem,4),tk0(jdim,idim-1,nummem,4), &
         ti0(jdim,kdim,nummem,4)
    ! velocity gradients
    REAL,INTENT(in) :: ux(jdim-1,kdim-1,idim-1,9)
    ! metrics
    REAL,INTENT(in) :: sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1), q(jdim,kdim,idim,5),&
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4),&
         dtj(jdim,kdim,idim-1)

    REAL,intent(in) :: vj0(kdim,idim-1,4),vk0(jdim,idim-1,4),&
         vi0(jdim,kdim,4)
    REAL,INTENT(in) :: x(jdim,kdim,idim),y(jdim,kdim,idim), &
         z(jdim,kdim,idim)
    INTEGER,intent(in) :: nbl

    REAL, INTENT(inout)::turb(jdim,kdim,idim,nummem),vist3d(jdim,kdim,idim),&
         zksav2(jdim,kdim,idim,2*nummem), smin(jdim-1,kdim-1,idim-1)
    INTEGER,intent(in) :: level,icyc
    REAL, INTENT(out) :: sumn(nummem)
    INTEGER,intent(out) :: negn(nummem)

    ! common block 
    common /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
             idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
             iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv 
    REAL :: title,rkap, xmach, alpha__,beta__,dt,fmax
    integer ::nit,ntt, &
             idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
             iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv 
    ! local variables, scalar ----
    INTEGER :: j,k,i,iv

    ! ---- local static allocatable arrays ----
    ! By defining the following arrays as "SAVE" type, the code
    ! may terminate by itself if the array has not been deallocated.
    ! This approach eliminates otherwise possible memory leak.

    ! local working array for stresses and zeta
    ! it includes the ghost cells
    REAL, SAVE, ALLOCATABLE:: turre(:,:,:,:) ! 6 (uiuj) + 1 (zeta)
    REAL, SAVE, ALLOCATABLE:: tke(:,:,:)   !turbulent kinetic energy
    REAL, SAVE, ALLOCATABLE:: omega(:,:,:,:) !vorticity vector
    REAL, SAVE, ALLOCATABLE:: bij(:,:,:,:)
    REAL, SAVE, ALLOCATABLE:: fmu(:,:,:)! laminar dynamic viscosity(with density units)
    REAL, SAVE, ALLOCATABLE:: timestep(:,:,:)
    REAL, SAVE, allocatable:: d(:,:,:,:,:) ! source term Jacobian
    REAL, SAVE, ALLOCATABLE,DIMENSION(:,:,:,:):: al,ar,bl,br
    REAL, SAVE, ALLOCATABLE,DIMENSION(:,:,:,:,:) :: dbijdx,dbijdxx
    INTEGER:: isub
    INTEGER,SAVE  :: icycle = 0
    INTEGER,save :: ntt_start=0
    INTEGER :: m
   
    if(allocated(src_k)) then   
       deallocate(src_k)
    endif
    allocate(src_k(jdim,kdim,idim))
    IF(icycle==0.AND.lkzstm_from_2eq) THEN
       ntt_start=ntt
    ENDIF
    icycle = icycle + 1
    IF(ntt==ntt_start.and.lkzstm_from_2eq) THEN
       CALL kws_init_stress_from_2eq(nbl,jdim,kdim,idim, nummem,ux,vist3d,q(1,1,1,1),turb)
       CALL kws_dump_movie(1,jdim,kdim,idim,nummem,nbl,q,x,y,z,vist3d,ux,turb,turb,smin)
       ! let the cfl3d to set up correct ghost cell values
       RETURN
    ENDIF
    
    
    !if(nres<2000) return
    !write(300,*)level,icyc,jdim
    mylevel =level
    myicyc = icyc

    myrank = myid
    myroot = myhost

#ifdef CHECK_SOURCE_TERM
    IF(ALLOCATED(source_items)) THEN
       deallocate(source_items)
    ENDIF
    ALLOCATE(source_items(9,nummem,jdim,kdim,idim))
#endif
    ! 6 stress components + 1 zeta, including ghost cells
    ALLOCATE(turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem))
    ALLOCATE(d(nummem,nummem,jdim,kdim,idim))
    ALLOCATE(al(2,jdim,kdim,idim),ar(2,jdim,kdim,idim),&
         bl(2,jdim,kdim,idim),br(2,jdim,kdim,idim))

    ! turbulent kinetic energy
    ALLOCATE(tke(0:jdim,0:kdim,0:idim))

    ! vorticity
    ALLOCATE(omega(0:jdim,0:kdim,0:idim,3))

    ! source terms
    IF(ALLOCATED(source)) THEN
       DEALLOCATE(source)
    ENDIF
    ALLOCATE(source(1:jdim-1,1:kdim-1,1:idim-1,nummem));
    if(allocated(rhs)) then
       deallocate(rhs)
    endif
    ALLOCATE(rhs(1:jdim-1,1:kdim-1,1:idim-1,nummem))


    ! bij terms
    ALLOCATE(bij(0:jdim, 0:kdim, 0:idim, 6))

    ! fmu array
    ALLOCATE(fmu(0:jdim,0:kdim,0:idim))

    ALLOCATE(timestep(jdim-1,kdim-1,idim-1))
    
    ! derivatives for c5 terms
    ALLOCATE(dbijdx (jdim-1,kdim-1,idim-1,6,3))
    ALLOCATE(dbijdxx(jdim-1,kdim-1,idim-1,6,3))
    


    ! copy inner cell data from turb array to the turre

    CALL fill_omega(jdim,kdim,idim, ux,omega)
    CALL fill_fmu(jdim,kdim,idim,q,qj0,qk0,qi0,fmu)
    CALL get_timestep(jdim,kdim,idim,dtj,vol,timestep,icyc)
       
    IF(icyc==1) THEN
       CALL save_lasttimestep(jdim,kdim,idim,nummem, turb,&
            zksav2(1,1,1,1),zksav2(1,1,1,nummem+1))
    ENDIF
    ! calculate the source terms for the 
    DO isub=1,nsubit
       source=0.0;
!      CALL set_wallbc(jdim,kdim,idim,nummem,q, turb, smin, vist3d, tj0,tk0,ti0)
       CALL fill_turre(jdim,kdim,idim, nummem, turb,tj0,tk0,ti0,turre)
       CALL fill_tke(jdim,kdim,idim, nummem, turre, tke, bij)
       call kws_dbij_dx(jdim,kdim,idim,nummem, turre,tke,sj,sk,si,vol,dbijdx)
       call kws_dbij_dxx(jdim,kdim,idim,nummem,dbijdx,sj,sk,si,vol,dbijdxx)
       CALL get_source(jdim,kdim,idim,nummem,q,qj0,qk0,qi0,turre,tke,bij,omega,&
            sj,sk,si,vol,vj0,vk0,vi0,ux,fmu,source,rhs,d,zksav2,timestep,dbijdxx,smin)

       al = 0;ar= 0;bl = 0; br=0
       select case(kdiff)
       case(0)  
         CALL get_diffusion(jdim,kdim,idim,nummem,q,qj0,qk0,qi0,turre,tke,&
              sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs,d,al,ar,bl,br)
       case(2)
         CALL kws_cijk(jdim,kdim,idim,nummem,q,qj0,qk0,qi0,turre,tke,&
              sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs,d,al,ar,bl,br)
       case(3)
         CALL kws_cijk0(jdim,kdim,idim,nummem,q,qj0,qk0,qi0,turre,tke,&
              sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs,d,al,ar,bl,br)
       case default
         stop "not supported diffusion options"
       end select
       
       CALL get_advection (jdim,kdim,idim,nummem,q,turre,sj,sk,si,vol,qj0,qk0,qi0,rhs,d,al,ar,bl,br)
       !CALL get_advection1 (jdim,kdim,idim,nummem,q,turre,sj,sk,si,vol,qj0,qk0,qi0,rhs)
       ! solve the (I - J dt )*\Delta Q  = Source*Dt
       IF(iopt<2.OR..TRUE.) THEN
       DO m=1,nummem
          DO i=1,idim-1
             DO k=1,kdim-1
                DO j=1,jdim-1
                   rhs(j,k,i,m) = rhs(j,k,i,m)*timestep(j,k,i)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
          CALL afsolver_k(jdim,kdim,idim,nummem,turre,q,qk0,fmu,tke,sk,vol,vk0,timestep,source,rhs)
          CALL afsolver_j(jdim,kdim,idim,nummem,turre,q,qj0,fmu,tke,sj,vol,vj0,timestep,source,rhs)
          CALL afsolver_i(jdim,kdim,idim,nummem,turre,q,qi0,fmu,tke,si,vol,vi0,timestep,source,rhs)
          CALL update(jdim,kdim,idim,nummem,rhs,turre,q,fmu, turb,vist3d,sumn,negn,ux)
       ELSE
          CALL sgs_solver_2d(jdim,kdim,idim,nummem, timestep, d,al,ar,bl,br,rhs)
          CALL update(jdim,kdim,idim,nummem,rhs,turre,q,fmu, turb,vist3d,sumn,negn,ux)
       ENDIF

    ENDDO
   
    IF(nfreq>0) THEN
       IF(MOD(ntt,nfreq)==0) THEN
          CALL kws_dump_movie(ntt/nfreq,jdim,kdim,idim,nummem,nbl,q,x,y,z,vist3d,ux,turb,turb,smin)
       ENDIF
!   ***dump movie file at end of cycle
!      IF(icyc==ncyc1(1) .or. icyc==ncyc1(2) .or. icyc==ncyc1(3) .or.&
!         icyc==ncyc1(4) .or. icyc==ncyc1(5)) then
!         CALL kws_dump_movie(0,jdim,kdim,idim,nummem,nbl,q,x,y,z,vist3d,ux,turb,turb,smin)
!      ENDIF
!   ***
    ENDIF
    DEALLOCATE(omega, tke, turre,bij,fmu,timestep,d,&
         al,ar,bl,br)
    DEALLOCATE(dbijdx,dbijdxx)
    !WRITE(200+myid,*) "kws_main called"
  END SUBROUTINE kws_main

  SUBROUTINE get_source(jdim,kdim,idim,nummem,q,&
       qj0,qk0,qi0, &
       turre,tke,bb,omega,&
       sj,sk,si,vol,vj0,vk0,vi0,ux,fmu,source,rhs,d,zksav,timestep,dbijdxx,smin)
    INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
    real, intent(in) :: smin(jdim-1,kdim-1,idim-1)
    REAL,INTENT(in) :: q(jdim,kdim,idim,5),&
         omega(0:jdim,0:kdim,0:idim,3), &
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4),&
         tke(0:jdim,0:kdim,0:idim), &
         sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),&
         vj0(kdim,idim-1,4),vk0(jdim,idim-1,4),&
         vi0(jdim,kdim,4), &
         ux(jdim-1,kdim-1,idim-1,9),fmu(0:jdim,0:kdim,0:idim),&
         zksav(jdim,kdim,idim,nummem),&
         timestep(jdim-1,kdim-1,idim-1),&
         dbijdxx(jdim-1,kdim-1,idim-1,6,3)

    real,intent(inout) :: &
         bb(0:jdim,0:kdim,0:idim,6), &
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),d(nummem,nummem,jdim,kdim,idim)

    REAL,INTENT(out) :: source(jdim-1,kdim-1,idim-1,nummem),rhs(jdim-1,kdim-1,idim-1,nummem)

    !common block needed
    common /lam/ ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
    integer :: ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
    common /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /twod/ i2d
    integer :: i2d
    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax
    REAL :: title,rkap,xmach,alpha__,beta__,dt,fmax
    COMMON /unst/ time,cfltau,ntstep,ita,iunst
    REAL :: time, cfltau
    INTEGER ::ntstep,ita,iunst
    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10


    ! tensor index to vector index translator for
    ! symmetry(and skew-symmetry) tensors.
    INTEGER, PARAMETER :: idxto6(3,3) = RESHAPE((/1,4,6,4,2,5,6,5,3/),(/3,3/))

    ! the sign associated with the Wij tensor
    REAL, PARAMETER:: WSGN(3,3) = RESHAPE((/0.,-1.,-1.,1.,0.,-1.,1.,1.,0./),(/3,3/))

    ! uij matrix index to ux vector index translator
    INTEGER, PARAMETER::idxux(3,3) = RESHAPE( (/1,4,7,2,5,8,3,6,9/),(/3,3/))

    ! Kroneck Delta
    REAL, PARAMETER :: krndelta(3,3) = RESHAPE((/1.,0.,0.,0.,1.,0.,0.,0.,1./),(/3,3/))
    ! vector index to tensor index translator
    INTEGER, PARAMETER:: idxto33_i(6) = (/1,2,3,1,2,1/)
    INTEGER, PARAMETER:: idxto33_j(6) = (/1,2,3,2,3,3/)

    REAL ::  wij(jdim,6) ! vorticity tensor (skew symmetric)
    REAL ::  sij(jdim,6) ! traceless deformation tensor
    REAL ::  sij_hat(jdim,6) ! sij - 0.5*duidxi * delta_ij (used in the 
    REAL ::  sij0(jdim,6)! deformation tensor
    REAL ::  eps(jdim)   ! dissipation of tke,  \mu * zeta  (with density units in it)
    REAL :: mu_t(jdim)    ! turbulent eddyviscosity
    REAL:: dpdi,dpdk,dpdj
    REAL,DIMENSION(jdim,3):: dpdx,dtkedx,dwdx
    REAL :: drdi,drdk,drdj,dtkedi,dtkedj,dtkedk,&
         dwdi,dwdj,dwdk
    REAL,DIMENSION(jdim,3):: drdx
    REAL :: dodi(3),dodj(3),dodk(3)
    REAL,DIMENSION(3,3,jdim) :: dodx
    INTEGER :: im,ip,jm,jp,km,kp


    INTEGER:: i,j,k,m,n,icur,jcur  ! icur and jcur are current stress tensor indices
    REAL :: prod,prdij(6),phi_ij,dissip
    REAL :: xma_re   ! Mach/Re
    REAL :: re_xma   ! Re/Mach
    REAL :: siave(3),sjave(3),skave(3)
    REAL :: coef_i ! related to i2d
    REAL :: bbij(6)   ! tau_ij + 2./3. tke * delta (ij)
    REAL :: dij,pij,pp,c1term, alpha_term, beta_term, gamma_term,beta
    REAL :: divg
    REAL :: pip,pim,pjp,pjm,pkp,pkm,rim,rip,rjm,rjp,rkm,rkp
    REAL :: phi(6)
    REAL :: prd_real(jdim)
    ! -- for omega equation --
    REAL :: tsum, kcur, x_w, f_beta, kdotw, sd_term
    integer,save:: ivisited=0
    ! -- for time step --
    real :: dtmp,cutoff

    coef_i=1.
    if(i2d==1) coef_i = 0

    xma_re = xmach/reue
    re_xma = reue/xmach

    wij(1:jdim,1:3) = 0
    if(ivisited==0) then
       ivisited=1
#ifdef DEBUG_KZSTM_IDXTO6
       write(*,*)'idxto6:'
       write(*,'(3(1x,I8))') idxto6(1,:) 
       write(*,'(3(1x,I8))') idxto6(2,:) 
       write(*,'(3(1x,I8))') idxto6(3,:) 
       write(*,*)'wsgn'
       write(*,'(3(1x,es12.5))') wsgn(1,:) 
       write(*,'(3(1x,es12.5))') wsgn(2,:) 
       write(*,'(3(1x,es12.5))') wsgn(3,:) 
       write(*,*)"idxux:"
       write(*,'(3(1x,I8))') idxux(1,:) 
       write(*,'(3(1x,I8))') idxux(2,:) 
       write(*,'(3(1x,I8))') idxux(3,:) 
       write(*,*)"krndelta:"
       write(*,'(3(1x,es12.5))') krndelta(1,:) 
       write(*,'(3(1x,es12.5))') krndelta(2,:) 
       write(*,'(3(1x,es12.5))') krndelta(3,:) 
#endif
    endif


    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             d(:,:, j,k,i) = 0  ! reset the jacobian matrix
             pkp = qk0(j,i,5,3)
             pkm = qk0(j,i,5,1)

             IF(k+1<=kdim-1) pkp = q(j,k+1,i,5)
             IF(k-1>=1) pkm=q(j,k-1,i,5)

             pjp = qj0(k,i,5,3)
             pjm = qj0(k,i,5,1)

             IF(j+1<=jdim-1) pjp = q(j+1,k,i,5)
             IF(j-1>=1) pjm=q(j-1,k,i,5)

             pip = qi0(j,k,5,3)
             pim = qi0(j,k,5,1)
             IF(i+1<=idim-1) pip = q(j,k,i+1,5)
             IF(i-1>=1) pim=q(j,k,i-1,5)

             dpdi=0.5*(pip-pim)*coef_i
             dpdj=0.5*(pjp-pjm)
             dpdk=0.5*(pkp-pkm)


             rkp = qk0(j,i,1,3)
             rkm = qk0(j,i,1,1)

             IF(k+1<=kdim-1) rkp = q(j,k+1,i,1)
             IF(k-1>=1) rkm=q(j,k-1,i,1)

             rjp = qj0(k,i,1,3)
             rjm = qj0(k,i,1,1)

             IF(j+1<=jdim-1) rjp = q(j+1,k,i,1)
             IF(j-1>=1) rjm=q(j-1,k,i,1)

             rip = qi0(j,k,1,3)
             rim = qi0(j,k,1,1)
             IF(i+1<=idim-1) rip = q(j,k,i+1,1)
             IF(i-1>=1) rim=q(j,k,i-1,1)


             drdi=0.5*(rip-rim)*coef_i
             drdj=0.5*(rjp-rjm)
             drdk=0.5*(rkp-rkm)

             dtkedi = 0.5*(tke(j,k,i+1)-tke(j,k,i-1))*coef_i
             dtkedj = 0.5*(tke(j+1,k,i)-tke(j-1,k,i))
             dtkedk = 0.5*(tke(j,k+1,i)-tke(j,k-1,i))

             dwdi = 0.5*(turre(j,k,i+1,7)-turre(j,k,i-1,7))*coef_i
             dwdj = 0.5*(turre(j+1,k,i,7)-turre(j-1,k,i,7))
             dwdk = 0.5*(turre(j,k+1,i,7)-turre(j,k-1,i,7))

             siave(1:3)=0.5*(si(j,k,i,1:3)*si(j,k,i,4)+si(j,k,i+1,1:3)*si(j,k,i+1,4))
             sjave(1:3)=0.5*(sj(j,k,i,1:3)*sj(j,k,i,4)+sj(j+1,k,i,1:3)*sj(j+1,k,i,4))
             skave(1:3)=0.5*(sk(j,k,i,1:3)*sk(j,k,i,4)+sk(j,k+1,i,1:3)*sk(j,k+1,i,4))

             siave = siave/vol(j,k,i)
             sjave = sjave/vol(j,k,i)
             skave = skave/vol(j,k,i)

             dpdx(j,:) = siave*dpdi + sjave*dpdj + skave*dpdk
             drdx(j,:) = siave*drdi + sjave*drdj + skave*drdk

             dtkedx(j,:) = siave*dtkedi+sjave*dtkedj+skave*dtkedk
             dwdx(j,:) =   siave*  dwdi+sjave*  dwdj+skave*  dwdk

             mu_t(j) = q(j,k,i,1)*tke(j,k,i)/turre(j,k,i,7)
             im=1;ip=1;jm=1;jp=1;km=1;kp=1
             IF(i==1) im = 0 !.AND.turre(j,k,i-2,1) == 1e20) im = 0
             IF(i==idim-1) ip = 0 !.AND.turre(j,k,i+2,1) == 1e20) ip = 0
             IF(j==1) jm = 0 !.AND.turre(j-2,k,i,1) == 1e20) jm = 0
             IF(j==jdim-1) jp = 0 !.AND.turre(j+2,k,i,1) == 1e20) jp = 0

             IF(k==1) km = 0!.AND.turre(j,k-2,i,1) == 1e20) km = 0
             IF(k==kdim-1) kp = 0!.AND.turre(j,k+2,i,1) == 1e20) kp = 0
          ENDDO

          DO j=1,jdim-1
             ! fill wij 
             wij(j,4) = 0.5*(ux(j,k,i,2) - ux(j,k,i,4)) ! w12=0.5*(uy - vx)
             wij(j,5) = 0.5*(ux(j,k,i,6) - ux(j,k,i,8)) ! w23=0.5*(vz - wy)
             wij(j,6) = 0.5*(ux(j,k,i,3) - ux(j,k,i,7)) ! w13=0.5*(uz - wx)

             sij(j,1) = ux(j,k,i,1)
             sij(j,2) = ux(j,k,i,5)
             sij(j,3) = ux(j,k,i,9)
             sij(j,4) = 0.5*(ux(j,k,i,2) + ux(j,k,i,4)) ! s12=0.5*(uy + vx)
             sij(j,5) = 0.5*(ux(j,k,i,6) + ux(j,k,i,8)) ! s23=0.5*(vz + wy)
             sij(j,6) = 0.5*(ux(j,k,i,3) + ux(j,k,i,7)) ! s13=0.5*(uz + wx)
             sij0(j,:) = sij(j,:)
             eps(j)   = beta_star*turre(j,k,i,7)*tke(j,k,i)
             divg     = ux(j,k,i,1)+ux(j,k,i,5)+ux(j,k,i,9)
             sij(j,1) = sij(j,1) - 1./3.*divg 
             sij(j,2) = sij(j,2) - 1./3.*divg 
             sij(j,3) = sij(j,3) - 1./3.*divg 
             
             sij_hat(j,:) = sij0(j,:)
             sij_hat(j,1) = sij_hat(j,1) - 0.5*divg 
             sij_hat(j,2) = sij_hat(j,2) - 0.5*divg 
             sij_hat(j,3) = sij_hat(j,3) - 0.5*divg 
          ENDDO
          DO j=1,jdim-1
             if ((i.ge.ilamlo .and. i.lt.ilamhi .and. &
                  j.ge.jlamlo .and. j.lt.jlamhi .and. &
                  k.ge.klamlo .and. k.lt.klamhi) .or. &
                  real(smin(j,k,i)) .lt. 0.) then
                cutoff=0.
             else
                cutoff=1.
             endif
             DO m=1,6
                icur= idxto33_i(m)
                jcur= idxto33_j(m)
                ! calculate production term
                prdij(m) = 0
                DO n=1,3
                   prdij(m) =&
                        turre(j,k,i,idxto6(icur,n))*ux(j,k,i,idxux(jcur,n))+&
                        turre(j,k,i,idxto6(jcur,n))*ux(j,k,i,idxux(icur,n))+prdij(m)
                ENDDO
                prdij(m) = prdij(m)*cutoff 
             enddo
             prd_real(j) = (prdij(1)+prdij(2)+prdij(3))*0.5*cutoff

             BBij(1:6)=turre(j,k,i,1:6)
             BBij(1:3)=BBij(1:3) + 2./3.*tke(j,k,i)

             DO m=1,6
                !IF(m>3) 
                !prd_coef = 1.0
                icur= idxto33_i(m)
                jcur= idxto33_j(m) 
                dij = 0
                pij = prdij(m)
                pp  = prd_real(j)
                DO n = 1,3
                   dij = dij + turre(j,k,i,idxto6(icur,n))*ux(j,k,i,idxux(n,jcur)) + &
                        turre(j,k,i,idxto6(jcur,n))*ux(j,k,i,idxux(n,icur))
                ENDDO
                if(m<=3) then
                  pij = pij - 2./3. *pp
                  dij = dij - 2./3. *pp 
                endif

                ! phi_ij term

                c1term = re_xma*c1*turre(j,k,i,7)*bbij(m)*beta_star
                alpha_term = alpha_hat*pij
                beta_term  = beta_hat*dij
                gamma_term = gamma_hat * tke(j,k,i)*sij(j,m)

                dissip = 2./3 *re_xma*eps(j)*krndelta(icur,jcur)
                phi_ij = (c1term - alpha_term - beta_term - gamma_term)
                phi(m) = phi_ij 
                source(j,k,i,m) = (-prdij(m)-phi_ij+dissip)!*vol(j,k,i)
#ifdef CHECK_SOURCE_TERM
                source_items(1,m,j,k,i) = -prdij(m)
                source_items(2,m,j,k,i) = dissip
                source_items(3,m,j,k,i) = c1term
                source_items(4,m,j,k,i) = phi_ij
                source_items(5,m,j,k,i) = alpha_term
                source_items(6,m,j,k,i) = beta_term
                source_items(7,m,j,k,i) = gamma_term
#endif

             ENDDO
          ENDDO
          DO j=1,jdim-1
             src_k(j,k,i) = -0.5*(source(j,k,i,1) +source(j,k,i,2) + source(j,k,i,3))*vol(j,k,i)
          enddo
          ! omega-equation source term
          DO j=1,jdim-1

             prod   = alpha *turre(j,k,i,7)/tke(j,k,i)*prd_real(j)
             tsum = 0
             DO icur = 1,3
                DO jcur = 1,3
                   DO kcur = 1,3
                      tsum = tsum + wij(j,idxto6(icur,jcur))*wsgn(icur,jcur)*&
                                    sij_hat(j,idxto6(kcur,icur))*&
                                    wij(j,idxto6(jcur,kcur))*wsgn(jcur,kcur)
                   ENDDO
                ENDDO
             ENDDO

             x_w = xma_re **3 *ABS( tsum/(beta_star*turre(j,k,i,7)**3))

             f_beta  = (1. + 85.*x_w)/(1.+100.*x_w)

             beta = beta_0 * f_beta

             dissip = Re_xma * beta * turre(j,k,i,7)**2

             kdotw = dtkedx(j,1)*dwdx(j,1) + dtkedx(j,2)*dwdx(j,2) + dtkedx(j,3)*dwdx(j,3)

             sd_term  = sigma_d0 *xma_re /turre(j,k,i,7) *max(0., kdotw)

             source(j,k,i,7)=(prod - dissip + sd_term )!*q(j,k,i,1)
          ENDDO
       ENDDO
    ENDDO
    rhs = source
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             !source(j,k,i,1:nummem) = &
             !source(j,k,i,1:nummem)*timestep(j,k,i)*(cfl_loc+cfl_psd)/cfl_psd
             !d(1:6,:,j,k,i) = 0
             d(7,:,j,k,i) = 0
             d(7,7,j,k,i) = 0.5*(rhs(j,k,i,7) - ABS(rhs(j,k,i,7)))/(turre(j,k,i,7))
          ENDDO
       ENDDO
    ENDDO
    ! check if it is time accurate calculation
    IF(dt>0) THEN
       IF(ita<0) THEN
          DO m=1,nummem
             DO i=1,idim-1
                DO k=1,kdim-1
                   DO j=1,jdim-1
                      if(ildts==1) then
                         dtmp = timestep(j,k,i)*(cfl_loc+cfl_psd)/cfl_psd
                      else
                         dtmp = dt
                      endif
                      rhs(j,k,i,m) = rhs(j,k,i,m) + (zksav(j,k,i,m)-turre(j,k,i,m))/dtmp
                      !d(m,m,j,k,i) =  d(m,m,j,k,i) + 1./dt
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO m=1,nummem
             DO i=1,idim-1
                DO k=1,kdim-1
                   DO j=1,jdim-1
                      rhs(j,k,i,m) = rhs(j,k,i,m) + (zksav(j,k,i,m)-turre(j,k,i,m))/dt
                      d(m,m,j,k,i) =  d(m,m,j,k,i)-1./dt
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
       ENDIF
    ENDIF
    !#endif

#ifdef CRAP
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0) THEN
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,1), "c","c11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,2), "c","c22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,3), "c","c33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,4), "c","c12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,5), "c","c23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,6), "c","c13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, source(1,1,1,7), "zet","zet.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,1), "ux1","ux1.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,2), "ux2","ux2.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,3), "ux3","ux3.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,4), "ux4","ux4.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,5), "ux5","ux5.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,6), "ux6","ux6.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,7), "ux7","ux7.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,8), "ux8","ux8.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim,1,kdim,1, ux(1,1,1,9), "ux9","ux9.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim-1,1,kdim-1,1, q(1,1,1,1), "q1","q1.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim-1,1,kdim-1,1, q(1,1,1,2), "q2","q2.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim-1,1,kdim-1,1, q(1,1,1,3), "q3","q3.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim-1,1,kdim-1,1, q(1,1,1,4), "q4","q4.plt")
       CALL dump_tecplot(jdim, kdim, idim, 1,jdim-1,1,kdim-1,1, q(1,1,1,5), "q5","q5.plt")
       CALL dump_tecplot(jdim+1, kdim+1, idim+1, 1,jdim+1,1,kdim+1,2, omega(0,0,0,1), "o1","omega1.plt")
       CALL dump_tecplot(jdim+1, kdim+1, idim+1, 1,jdim+1,1,kdim+1,2, omega(0,0,0,2), "o2","omega2.plt")
       CALL dump_tecplot(jdim+1, kdim+1, idim+1, 1,jdim+1,1,kdim+1,2, omega(0,0,0,3), "o3","omega3.plt")
    ENDIF
#endif
  END SUBROUTINE get_source
  
  SUBROUTINE get_diffusion(jdim,kdim,idim,nummem,q, &
       qj0,qk0,qi0, &
       turre,tke,&
       sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs,d,&
       al,ar,bl,br)
    INTEGER,INTENT(in) :: jdim,idim,kdim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5),&
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),&
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4),&
         tke(0:jdim,0:kdim,0:idim), &
         sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),&
         vj0(kdim,idim-1,4),vk0(jdim,idim-1,4),&
         vi0(jdim,kdim,4), &
         fmu(0:jdim,0:kdim,0:idim)
    REAL, INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)
    REAL, INTENT(inout) :: d(nummem,nummem,jdim,kdim,idim)
    REAL, INTENT(out), DIMENSION(2,jdim,kdim,idim) :: al,ar,bl,br

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
    COMMON /twod/ i2d
    INTEGER :: i2d

    !local variables
    INTEGER :: i,j,k,m
  
    REAL:: diff(MAX(jdim,kdim,idim),nummem)
    REAL:: xmut(MAX(jdim,kdim,idim))
    REAL:: xmu_ave,xmut_ave,vleft,vright,tke_ave,ome_ave,rho_ave
    REAL:: areax2,areay2,areaz2,area2,rvol,xma_re
    
    REAL, parameter :: one_third = 0.3333333333333333
    REAL, PARAMETER :: two_third = 2./3.
    REAL :: rv_l, rv_r,rv_c,xl,yl,zl, xr,yr,zr,xc,yc,zc
    REAL :: xmu,flux(max(jdim,kdim,idim),nummem)
    REAL :: xcoef(max(jdim,kdim,idim),nummem)

    xma_re= xmach/reue
    ! diffusion terms in the j-direction
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             IF(j==1) THEN
                rho_ave = 0.5*(qj0(k,i,1,1)+q(j,k,i,1))
             ELSEif(j==jdim) then
                rho_ave  = 0.5*(qj0(k,i,1,3)+q(j-1,k,i,1))
             ELSE
                rho_ave  = 0.5*(q(j-1,k,i,1)+q(j,k,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j-1,k,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j-1,k,i,7))
             xmut(j) = MAX(0.0,rho_ave*tke_ave/ome_ave)
             diff(j,:) =turre(j,k,i,:) -  turre(j-1,k,i,:)
             IF(j==1) THEN
                vleft = vj0(k,i,1)
             ELSE
                vleft = vol(j-1,k,i)
             ENDIF
             IF(j==jdim) THEN
                vright = vj0(k,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             xmu = xmu_ave+xmut(j)*sigma_star
             flux(j,1:6) = diff(j,1:6)*sj(j,k,i,4)**2*rvol*xmu
             xcoef(j,1) = sj(j,k,i,4)**2*rvol*xmu 

             xmu = xmu_ave+xmut(j)*sigma
             flux(j,7) = diff(j,7)*sj(j,k,i,4)**2*rvol*xmu
             xcoef(j,2) = sj(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !CALL cap_xmut(xmut,0,jdim) 
          DO j=1,jdim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(j+1,:)-flux(j,:))/(q(j,k,i,1)*vol(j,k,i))

             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(j+1,2) + xcoef(j,2))/vol(j,k,i) 
             al(2,j,k,i) = -xma_re*xcoef(j,2)/vol(j,k,i)
             ar(2,j,k,i) = -xma_re*xcoef(j+1,2)/vol(j,k,i)
             
             al(1,j,k,i) = -xma_re*xcoef(j,1)/vol(j,k,i)
             ar(1,j,k,i) = -xma_re*xcoef(j+1,1)/vol(j,k,i)
             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(j,1)+xcoef(j+1,1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! diffusion terms in the k-direction
    DO i=1,idim-1
       DO j=1,jdim-1
          DO k=1,kdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             IF(k==1) THEN
                rho_ave = 0.5*(qk0(j,i,1,1)+q(j,k,i,1))
             ELSEIF(k==kdim) THEN
                rho_ave  = 0.5*(qk0(j,i,1,3)+q(j,k-1,i,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k-1,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k-1,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k-1,i,7))
             xmut(k) = max(0.0,rho_ave*tke_ave/ome_ave)
             diff(k,:) =turre(j,k,i,:) -  turre(j,k-1,i,:)
             IF(k==1) THEN
                vleft = vk0(j,i,1)
             ELSE
                vleft = vol(j,k-1,i)
             ENDIF
             IF(k==kdim) THEN
                vright = vk0(j,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             xmu = xmu_ave+xmut(k)*sigma_star
             flux(k,1:6) = diff(k,1:6)*sk(j,k,i,4)**2*rvol*xmu
             xcoef(k,1) = sk(j,k,i,4)**2*rvol*xmu 
             ! c5 term
             !xmu = xmu_ave*one_third*0.0 + c5*xmut(k)
             !flux(k,1:3) = flux(k,1:3)+0.5*(turre(j,k,i,1:3)-turre(j,k-1,i,1:3) +two_third* &
             !     (tke(j,k,i)-tke(j,k-1,i)))*sk(j,k,i,4)**2*rvol*xmu
             !flux(k,4:6) = flux(k,4:6)+ 0.5*(turre(j,k,i,4:6)-turre(j,k-1,i,4:6))*sk(j,k,i,4)**2*rvol*xmu
             !xcoef(k,1) = xcoef(k,1) + sk(j,k,i,4)**2*rvol*xmu

             xmu = xmu_ave+xmut(k)*sigma
             flux(k,7) = diff(k,7)*sk(j,k,i,4)**2*rvol*xmu
             xcoef(k,2) = sk(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO k=1,kdim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(k+1,:)-flux(k,:))/(q(j,k,i,1)*vol(j,k,i))
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(k+1,2) + xcoef(k,2))/vol(j,k,i) 
             bl(2,j,k,i) = -xma_re*xcoef(k,2)/vol(j,k,i)
             br(2,j,k,i) = -xma_re*xcoef(k+1,2)/vol(j,k,i)
             
             bl(1,j,k,i) = -xma_re*xcoef(k,1)/vol(j,k,i)
             br(1,j,k,i) = -xma_re*xcoef(k+1,1)/vol(j,k,i)
             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(k,1)+xcoef(k+1,1))/vol(j,k,i)
             ENDDO
                          
          ENDDO
       ENDDO
    ENDDO
    IF(i2d/=1) THEN
    ! diffusion terms in the i-direction
    DO k=1,kdim-1
       DO j=1,jdim-1
          DO i=1,idim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             IF(i==1) THEN
                rho_ave = 0.5*(qi0(j,k,1,1)+q(j,k,i,1))
             ELSEIF(i==idim) THEN
                rho_ave  = 0.5*(qi0(j,k,1,3)+q(j,k,i-1,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k,i-1,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k,i-1))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k,i-1,7))
             xmut(i) = MAX(0.,rho_ave*tke_ave/ome_ave)

             diff(i,:) =turre(j,k,i,:) -  turre(j,k,i-1,:)
             IF(i==1) THEN
                vleft = vi0(j,k,1)
             ELSE
                vleft = vol(j,k,i-1)
             ENDIF
             IF(i==idim) THEN
                vright = vi0(j,k,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             xmu = xmu_ave+xmut(i)*sigma_star
             flux(i,1:6) = diff(i,1:6)*si(j,k,i,4)**2*rvol*xmu
             xcoef(i,1) = si(j,k,i,4)**2*rvol*xmu 
             xmu = xmu_ave+xmut(i)*sigma
             flux(i,7) = diff(i,7)*si(j,k,i,4)**2*rvol*xmu
             xcoef(i,2) = si(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO i=1,idim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(i+1,:)-flux(i,:))/(q(j,k,i,1)*vol(j,k,i))
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(i+1,2) + xcoef(i,2))/vol(j,k,i) 

             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(i,1)+xcoef(i+1,1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO    
    ENDIF
#ifdef DEBUG_KWSTM_GETDIFFU
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.or..false.) THEN
       write(*,*)"zeta=", turre(20,1:3,1,7)
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rdif","r-dif11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rdif","r-dif22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rdif","r-dif33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rdif","r-dif12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rdif","r-dif23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rdif","r-dif13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rdif","r-dif44.plt")
    ENDIF
#endif
    RETURN
  END SUBROUTINE get_diffusion
  SUBROUTINE get_diffusion_old(jdim,kdim,idim,nummem,q, &
       qj0,qk0,qi0, &
       turre,tke,&
       sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs)
    INTEGER,INTENT(in) :: jdim,idim,kdim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5),&
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),&
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4),&
         tke(0:jdim,0:kdim,0:idim), &
         sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),&
         vj0(kdim,idim-1,4),vk0(jdim,idim-1,4),&
         vi0(jdim,kdim,4), &
         fmu(0:jdim,0:kdim,0:idim)
    REAL, INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
    COMMON /twod/ i2d
    INTEGER :: i2d

    !local variables
    INTEGER :: i,j,k,m
  
    REAL:: diff(MAX(jdim,kdim,idim),nummem)
    REAL:: xmut(MAX(jdim,kdim,idim))
    REAL:: xmu_ave,xmut_ave,vleft,vright,tke_ave,ome_ave,rho_ave
    REAL:: areax2,areay2,areaz2,area2,rvol,xma_re
    
    REAL, parameter :: one_third = 0.3333333333333333
    REAL, PARAMETER :: two_third = 2./3.

    xma_re= xmach/reue
    ! diffusion terms in the j-direction
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             IF(j==1) THEN
                rho_ave = 0.5*(qj0(k,i,1,1)+q(j,k,i,1))
             ELSEif(j==jdim) then
                rho_ave  = 0.5*(qj0(k,i,1,3)+q(j-1,k,i,1))
             ELSE
                rho_ave  = 0.5*(q(j-1,k,i,1)+q(j,k,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j-1,k,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j-1,k,i,7))
             xmut(j) =MAX(0.0, rho_ave*tke_ave/ome_ave)

          ENDDO

          !CALL cap_xmut(xmut,0,jdim) 
          DO j=1,jdim
             IF(j==1) THEN
                vleft = vj0(k,i,1)
             ELSE
                vleft = vol(j-1,k,i)
             ENDIF
             IF(j==jdim) THEN
                vright = vj0(k,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             
             areax2= (sj(j,k,i,1)*sj(j,k,i,4))**2
             areay2= (sj(j,k,i,2)*sj(j,k,i,4))**2
             areaz2= (sj(j,k,i,3)*sj(j,k,i,4))**2
             area2 = (areax2 + areay2 + areaz2)*rvol
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             xmut_ave= xmut(j)
             ! diffusion term in enstrophy equation
             diff(j,7) = xma_re*(xmu_ave+sigma*xmut_ave)*(turre(j,k,i,7) - turre(j-1,k,i,7))*area2
             
             diff(j,1:6) = xma_re*(xmu_ave+sigma_star*xmut_ave)*(turre(j,k,i,1:6)-turre(j-1,k,i,1:6))*area2
             
             ! the diffusion term hiding in the Phi_ij
             !diff(j,1:3) = diff(j,1:3)+0.5*xma_re*(xmu_ave+c5*xmut_ave)*(turre(j,k,i,1:3)-turre(j-1,k,i,1:3) +two_third* &
             !     (tke(j,k,i)-tke(j-1,k,i)))*area2  
             !diff(j,4:6) = diff(j,4:6)+ 0.5*xma_re*(xmu_ave+c5*xmut_ave)*(turre(j,k,i,4:6)-turre(j-1,k,i,4:6))*area2
          ENDDO
          DO m=1,7
             DO j=1,jdim-1
                rhs(j,k,i,m) = rhs(j,k,i,m)+(diff(j+1,m)-diff(j,m))/(q(j,k,i,1)*vol(j,k,i))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! diffusion terms in the k-direction
    DO i=1,idim-1
       DO j=1,jdim-1
          DO k=1,kdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             IF(k==1) THEN
                rho_ave = 0.5*(qk0(j,i,1,1)+q(j,k,i,1))
             ELSEIF(k==kdim) THEN
                rho_ave  = 0.5*(qk0(j,i,1,3)+q(j,k-1,i,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k-1,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k-1,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k-1,i,7))
             xmut(k) = rho_ave*tke_ave/(ome_ave)
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO k=1,kdim
             IF(k==1) THEN
                vleft = vk0(j,i,1)
             ELSE
                vleft = vol(j,k-1,i)
             ENDIF
             IF(k==kdim) THEN
                vright = vk0(j,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             
             areax2= (sk(j,k,i,1)*sk(j,k,i,4))**2
             areay2= (sk(j,k,i,2)*sk(j,k,i,4))**2
             areaz2= (sk(j,k,i,3)*sk(j,k,i,4))**2
             area2 = (areax2 + areay2 + areaz2)*rvol
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             xmut_ave= xmut(k)

             ! diffusion term in enstrophy equation
             diff(k,7) = xma_re*(xmu_ave+sigma*xmut_ave)*(turre(j,k,i,7) - turre(j,k-1,i,7))*area2
             diff(k,1:6) = xma_re*(xmu_ave+sigma_star*xmut_ave)*(turre(j,k,i,1:6)-turre(j,k-1,i,1:6))*area2
                          
             ! the diffusion term hiding in the Phi_ij
            ! diff(k,1:3) = diff(k,1:3)+0.5*xma_re*(xmu_ave+c5*xmut_ave)*(turre(j,k,i,1:3)-turre(j,k-1,i,1:3) + two_third*&
            !      (tke(j,k,i)-tke(j,k-1,i)))*area2  
            ! diff(k,4:6) = diff(k,4:6)+ 0.5*xma_re*(xmu_ave+c5*xmut_ave)*(turre(j,k,i,4:6)-turre(j,k-1,i,4:6))*area2
          ENDDO
          DO m=1,7
             DO k=1,kdim-1
                rhs(j,k,i,m) = rhs(j,k,i,m)+(diff(k+1,m)-diff(k,m))/(q(j,k,i,1)*vol(j,k,i))
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    IF(i2d/=1) THEN
    ! diffusion terms in the i-direction
    DO k=1,kdim-1
       DO j=1,jdim-1
          DO i=1,idim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             IF(i==1) THEN
                rho_ave = 0.5*(qi0(j,k,1,1)+q(j,k,i,1))
             ELSEIF(i==idim) THEN
                rho_ave  = 0.5*(qi0(j,k,1,3)+q(j,k,i-1,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k,i-1,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k,i-1))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k,i-1,7))
             xmut(i) = rho_ave*tke_ave/(ome_ave)
          ENDDO
          !CALL cap_xmut(xmut,0,idim) 
          DO i=1,idim
             IF(i==1) THEN
                vleft = vi0(j,k,1)
             ELSE
                vleft = vol(j,k,i-1)
             ENDIF
             IF(i==idim) THEN
                vright = vi0(j,k,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             
             areax2= (si(j,k,i,1)*si(j,k,i,4))**2
             areay2= (si(j,k,i,2)*si(j,k,i,4))**2
             areaz2= (si(j,k,i,3)*si(j,k,i,4))**2
             area2 = (areax2 + areay2 + areaz2)*rvol
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             xmut_ave= xmut(i)
             ! diffusion term in enstrophy equation
             diff(i,7) = xma_re*(xmu_ave+sigma*xmut_ave)*(turre(j,k,i,7) - turre(j,k,i-1,7))*area2
             diff(i,1:6) = xma_re*(xmu_ave+sigma_star*xmut_ave)*(turre(j,k,i,1:6)-turre(j,k,i-1,1:6))*area2
             
             ! the diffusion term hiding in the Phi_ij
             !diff(i,1:3) = diff(i,1:3)+0.5*xma_re*(xmu_ave+c5*xmut_ave)*(turre(j,k,i,1:3)-turre(j,k,i-1,1:3) +two_third* &
             !     (tke(j,k,i)-tke(j,k,i-1)))*area2
             !diff(i,4:6) = diff(i,4:6)+0.5*xma_re*(xmu_ave+c5*xmut_ave)*(turre(j,k,i,4:6)-turre(j,k,i-1,4:6))*area2
          ENDDO
          DO m=1,7
             DO i=1,idim-1
                rhs(j,k,i,m) = rhs(j,k,i,m)+(diff(i+1,m)-diff(i,m))/(q(j,k,i,1)*vol(j,k,i))
             ENDDO
          ENDDO
       ENDDO
    ENDDO    
    ENDIF
#ifdef CRAP
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0) THEN
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rdif","r-dif11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rdif","r-dif22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rdif","r-dif33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rdif","r-dif12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rdif","r-dif23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rdif","r-dif13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rdif","r-difzt.plt")
    ENDIF
#endif
    RETURN
  END SUBROUTINE get_diffusion_old


  subroutine cap_xmut(xmut,istr,iend)
     integer,intent(in) :: istr,iend
     real,intent(inout) :: xmut(istr:iend)
     
     integer :: i
     return
     do i=istr,iend
        xmut(i) = SIGN(MAX(ABS(xmut(i)),0.01),xmut(i))!min(xmut(i),10.) 
     enddo 
  end subroutine
  SUBROUTINE afsolver_i(jdim,kdim,idim,nummem,turre,q,qi0,fmu,tke,si,vol,vi0,timestep,source,rhs)
    INTEGER, INTENT(in)::jdim,kdim,idim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5), &
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem), &
         si(jdim,kdim,idim,5),&
         vol(jdim,kdim,idim-1),&
         vi0(jdim,kdim,4),&
         qi0(jdim,kdim,5,4),&
         fmu(0:jdim,0:kdim,0:idim),source(jdim-1,kdim-1,idim-1,nummem),&
         timestep(jdim-1,kdim-1,idim-1),tke(0:jdim,0:kdim,0:idim)
    
    REAL,INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
    COMMON /twod/ i2d
    integer :: i2d
    
    ! three diagonal terms
    REAL :: rlow(idim-1,nummem),diag(idim-1,nummem),rup(idim-1,nummem),rcol(idim-1)
    REAL :: up(idim-1),um(idim-1),uu
    REAL :: xmet(idim)   ! metrics term for diffusion
    REAL :: xmu_ave,xmut(idim)
    INTEGER :: i,j,k,m
    REAL:: xc,yc,zc,tc
    REAL:: rhov, diff7(2),diff13(2),diff46(2),xma_re,sol(idim-1)
    REAL:: rho_ave,tke_ave,ome_ave
    REAL :: rv_l, rv_r,rv_c,xl,yl,zl, xr,yr,zr,vleft,vright
    REAL :: xlc,xrc,xmu_l, xmu_r, tll, trr, tcc,amul,amur

    if(i2d==1) return 
    xma_re = xmach/reue
    DO j=1,jdim-1
       DO k=1,kdim-1
          DO i=1,idim-1
             xc=0.5*(si(j,k,i+1,1)*si(j,k,i+1,4)+ &
                     si(j,k,i  ,1)*si(j,k,i  ,4))/vol(j,k,i)
             yc=0.5*(si(j,k,i+1,2)*si(j,k,i+1,4)+ &
                     si(j,k,i  ,2)*si(j,k,i  ,4))/vol(j,k,i)
             zc=0.5*(si(j,k,i+1,3)*si(j,k,i+1,4)+ &
                     si(j,k,i  ,3)*si(j,k,i  ,4))/vol(j,k,i)
             tc=0.5*(si(j,k,i+1,5)*si(j,k,i+1,4)+ &
                     si(j,k,i  ,5)*si(j,k,i  ,4))/vol(j,k,i)
             uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc

             up(i) = 0.5*(uu+abs(uu))
             um(i) = 0.5*(uu-abs(uu))
          ENDDO

          ! convection terms - 1st order(time and space)
          DO m=1,nummem
             DO i=1,idim-1
                rlow(i,m)= -up(i)
                rup(i,m) = um(i)
                diag(i,m) =up(i)-um(i)
             ENDDO
          ENDDO

          DO i=1,idim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             IF(i==1) THEN
                rho_ave = 0.5*(qi0(j,k,1,1)+q(j,k,i,1))
             ELSEIF(i==idim) THEN
                rho_ave  = 0.5*(qi0(j,k,1,3)+q(j,k,i-1,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k,i-1,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k,i-1))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k,i-1,7))
             xmut(i) = MAX(0.,rho_ave*tke_ave/ome_ave)
          ENDDO

          DO i=1,idim-1
             IF(i==1) THEN
                vleft = vi0(j,k,1)
             ELSE
                vleft = vol(j,k,i-1)
             ENDIF
             IF(i==idim-1) THEN
                vright = vi0(j,k,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rv_l = 1./(0.5*(vleft + vol(j,k,i)))
             rv_r = 1./(0.5*(vright + vol(j,k,i)))
             rv_c = 1./vol(j,k,i)
             
             xl= (si(j,k,i,1)*si(j,k,i,4))*rv_l
             yl= (si(j,k,i,2)*si(j,k,i,4))*rv_l
             zl= (si(j,k,i,3)*si(j,k,i,4))*rv_l

             xr= (si(j,k,i+1,1)*si(j,k,i+1,4))*rv_r
             yr= (si(j,k,i+1,2)*si(j,k,i+1,4))*rv_r
             zr= (si(j,k,i+1,3)*si(j,k,i+1,4))*rv_r

             xc = 0.5*(si(j,k,i,1)*si(j,k,i,4)+si(j,k,i+1,1)*si(j,k,i+1,4))*rv_c
             yc = 0.5*(si(j,k,i,2)*si(j,k,i,4)+si(j,k,i+1,2)*si(j,k,i+1,4))*rv_c
             zc = 0.5*(si(j,k,i,3)*si(j,k,i,4)+si(j,k,i+1,3)*si(j,k,i+1,4))*rv_c
      
             xlc = xl*xc + yl*yc + zl*zc
             xrc = xr*xc + yr*yc + zr*zc
 
             amul=0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             amur=0.5*(fmu(j,k,i)+fmu(j,k,i+1))
             xmu_l = amul+xmut(i)*sigma
             xmu_r = amur+xmut(i+1)*sigma
             trr   = xmu_r*xrc
             tll   = xmu_l*xlc
             tcc   = trr+tll

             rlow(i,7) = rlow(i,7) - tll*xma_re/q(j,k,i,1)
             diag(i,7) = diag(i,7) + (trr+tll)*xma_re/q(j,k,i,1)
             rup(i,7) = rup(i,7) - trr*xma_re/q(j,k,i,1)
            
             xmu_l = amul+xmut(i)*sigma_star
             xmu_r = amur+xmut(i+1)*sigma_star
             trr   = xmu_r*xrc
             tll   = xmu_l*xlc
             tcc   = trr+tll

             rlow(i,1:6) = rlow(i,1:6) - tll*xma_re/q(j,k,i,1)
             diag(i,1:6) = diag(i,1:6) + (trr+tll)*xma_re/q(j,k,i,1)
             rup(i,1:6) = rup(i,1:6) - trr*xma_re/q(j,k,i,1)
          enddo 

          DO m=1,nummem
             DO i=1,idim-1
                rlow(i,m) = rlow(i,m)*timestep(j,k,i)
                rup(i,m) = rup(i,m)*timestep(j,k,i)
                diag(i,m) = diag(i,m)*timestep(j,k,i) + 1.
             ENDDO
          ENDDO
          
          DO m=1,nummem
             DO i=1,idim-1
                rcol(i) = rhs(j,k,i,m)
             ENDDO
             CALL tri_solver(idim-1,rlow(1,m),diag(1,m),rup(1,m),rcol,sol)
             DO i=1,idim-1
                rhs(j,k,i,m) = sol(i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE afsolver_i

  ! tri-diagonal solver in the j direction
  SUBROUTINE afsolver_j(jdim,kdim,idim,nummem,turre,q,qj0,fmu,tke,sj,vol,vj0,timestep,source,rhs)
    INTEGER, INTENT(in)::jdim,kdim,idim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5), &
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem), &
         sj(jdim,kdim,idim-1,5),&
         !sk(jdim,kdim,idim-1,5),&
         !si(jdim,kdim,idim,5),&
         vol(jdim,kdim,idim-1),&
         vj0(kdim,idim-1,4),&
         qj0(kdim,idim-1,5,4),&
         !vk0(jdim,idim-1,4),&
         !vi0(jdim,kdim,4), &
         fmu(0:jdim,0:kdim,0:idim),source(jdim-1,kdim-1,idim-1,nummem),&
         timestep(jdim-1,kdim-1,idim-1),tke(0:jdim,0:kdim,0:idim)
    
    REAL,INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
   
    ! three diagonal terms
    REAL :: rlow(jdim-1,nummem),diag(jdim-1,nummem),rup(jdim-1,nummem),rcol(jdim-1)
    REAL :: up(jdim-1),um(jdim-1),uu
    REAL :: xmu_ave,xmut(jdim)
    INTEGER :: i,j,k,m
    REAL:: xc,yc,zc,tc
    REAL:: rhov, xma_re,sol(jdim-1)
    REAL:: rho_ave,tke_ave,ome_ave
    REAL :: rv_l, rv_r,rv_c,xl,yl,zl, xr,yr,zr,vleft,vright
    REAL :: xlc,xrc,xmu_l, xmu_r, tll, trr, tcc,amul,amur
    REAL :: re_xma 
    xma_re = xmach/reue; re_xma = 1./xma_re
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             xc=0.5*(sj(j+1,k,i,1)*sj(j+1,k,i,4)+ &
                  sj(j,k,i  ,1)*sj(j,k,i  ,4))/vol(j,k,i)
             yc=0.5*(sj(j+1,k,i,2)*sj(j+1,k,i,4)+ &
                  sj(j,k,i  ,2)*sj(j,k,i  ,4))/vol(j,k,i)
             zc=0.5*(sj(j+1,k,i,3)*sj(j+1,k,i,4)+ &
                  sj(j,k,i  ,3)*sj(j,k,i  ,4))/vol(j,k,i)
             tc=0.5*(sj(j+1,k,i,5)*sj(j+1,k,i,4)+ &
                  sj(j,k,i  ,5)*sj(j,k,i  ,4))/vol(j,k,i)
             uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc

             up(j) = 0.5*(uu+abs(uu))
             um(j) = 0.5*(uu-abs(uu))
          ENDDO
          DO j=1,jdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             IF(j==1) THEN
                rho_ave = 0.5*(qj0(k,i,1,1)+q(j,k,i,1))
             ELSEif(j==jdim) then
                rho_ave  = 0.5*(qj0(k,i,1,3)+q(j-1,k,i,1))
             ELSE
                rho_ave  = 0.5*(q(j-1,k,i,1)+q(j,k,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j-1,k,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j-1,k,i,7))
             xmut(j) = MAX(0.0,rho_ave*tke_ave/ome_ave)
          ENDDO

          ! convection terms - 1st order(time and space)
          DO m=1,nummem
             DO j=1,jdim-1
                rlow(j,m)= -up(j)
                rup(j,m) = um(j)
                diag(j,m) =up(j)-um(j)
             ENDDO
          ENDDO
          DO j=1,jdim-1
             IF(j==1) THEN
                vleft = vj0(k,i,1)
             ELSE
                vleft = vol(j-1,k,i)
             ENDIF
             IF(j==jdim-1) THEN
                vright = vj0(k,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF

             rv_l = 1./(0.5*(vleft + vol(j,k,i)))
             rv_r = 1./(0.5*(vright + vol(j,k,i)))
             rv_c = 1./vol(j,k,i)
             
             xl= (sj(j,k,i,1)*sj(j,k,i,4))*rv_l
             yl= (sj(j,k,i,2)*sj(j,k,i,4))*rv_l
             zl= (sj(j,k,i,3)*sj(j,k,i,4))*rv_l

             xr= (sj(j+1,k,i,1)*sj(j+1,k,i,4))*rv_r
             yr= (sj(j+1,k,i,2)*sj(j+1,k,i,4))*rv_r
             zr= (sj(j+1,k,i,3)*sj(j+1,k,i,4))*rv_r

             xc = 0.5*(sj(j,k,i,1)*sj(j,k,i,4)+sj(j+1,k,i,1)*sj(j+1,k,i,4))*rv_c
             yc = 0.5*(sj(j,k,i,2)*sj(j,k,i,4)+sj(j+1,k,i,2)*sj(j+1,k,i,4))*rv_c
             zc = 0.5*(sj(j,k,i,3)*sj(j,k,i,4)+sj(j+1,k,i,3)*sj(j+1,k,i,4))*rv_c
      
             xlc = xl*xc + yl*yc + zl*zc
             xrc = xr*xc + yr*yc + zr*zc
 
             amul=0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             amur=0.5*(fmu(j,k,i)+fmu(j+1,k,i))
             xmu_l = amul+xmut(j)*sigma
             xmu_r = amur+xmut(j+1)*sigma
             trr   = xmu_r*xrc
             tll   = xmu_l*xlc
             tcc   = trr+tll
            
             rlow(j,7) = rlow(j,7)-tll/q(j,k,i,1)*xma_re
             diag(j,7) = diag(j,7)+tcc/q(j,k,i,1)*xma_re
             rup(j,7) =   rup(j,7)-trr/q(j,k,i,1)*xma_re
       
             xmu_l = amul+xmut(j)*sigma_star
             xmu_r = amur+xmut(j+1)*sigma_star
             trr   = xmu_r*xrc
             tll   = xmu_l*xlc
             tcc   = trr+tll
             ! for stresses
             rlow(j,1:6) = rlow(j,1:6)- tll*xma_re/q(j,k,i,1)
             diag(j,1:6) = diag(j,1:6)+ tcc*xma_re/q(j,k,i,1)
             rup(j,1:6) =  rup(j,1:6) - trr*xma_re/q(j,k,i,1)
          ENDDO

          
          DO m=1,nummem
             DO j=1,jdim-1
                rlow(j,m) = rlow(j,m)*timestep(j,k,i)
                rup(j,m) = rup(j,m)*timestep(j,k,i)
                diag(j,m) = diag(j,m)*timestep(j,k,i) + 1.
             ENDDO
          ENDDO
          DO m=1,nummem
             DO j=1,jdim-1
                rcol(j) = rhs(j,k,i,m)
             ENDDO
#ifdef DEBUG_KZSTM
             WRITE(2000+myrank,*)"m=",m,i,j,k
             WRITE(2000+myrank,'(A,1x,5(1x,es12.5))')"l D u xmut=", &
                     (rlow(j,m),diag(j,m),rup(j,m),xmut_ave(j),rcol(j),j=1,jdim-1)
             WRITE(2000+myrank,'(A,1x,2(1x,es12.5))')"xmut,xmet=",(xmut_ave(j),xmet(j),j=1,jdim)
             WRITE(2000+myrank,'(A,1x,7(1x,es12.5))')"turre=",&
                  (turre(j,k,i,1:7),j=-1,jdim+1)
             
             CLOSE(2000+myrank)
#endif
             CALL tri_solver(jdim-1,rlow(1,m),diag(1,m),rup(1,m),rcol,sol)
             DO j=1,jdim-1
                rhs(j,k,i,m) = sol(j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

#ifdef DEBUG_KWSTM_SOLVER_J
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.or..false.) THEN
       write(*,*)"solver_j=", rhs(1:4,1,1,6)
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rdif","j-dif11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rdif","j-dif22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rdif","j-dif33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rdif","j-dif12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rdif","j-dif23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rdif","j-dif13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rdif","j-dif44.plt")
    ENDIF
#endif
  END SUBROUTINE afsolver_j

  ! tri-diagonal solver in the k direction
  SUBROUTINE afsolver_k(jdim,kdim,idim,nummem,turre,q,qk0,fmu,tke,sk,vol,vk0,timestep,source,rhs)
    INTEGER, INTENT(in)::jdim,kdim,idim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5), &
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem), &
         sk(jdim,kdim,idim-1,5),&
         !si(jdim,kdim,idim,5),&
         vol(jdim,kdim,idim-1),&
         vk0(jdim,idim-1,4),&
         qk0(jdim,idim-1,5,4),&
         !vk0(jdim,idim-1,4),&
         !vi0(jdim,kdim,4), &
         fmu(0:jdim,0:kdim,0:idim),source(jdim-1,kdim-1,idim-1,nummem),&
         timestep(jdim-1,kdim-1,idim-1),tke(0:jdim,0:kdim,0:idim)
    
    REAL,INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
    
    ! three diagonal terms
    REAL :: rlow(kdim-1,nummem),diag(kdim-1,nummem),rup(kdim-1,nummem),rcol(kdim-1)
    REAL :: up(kdim-1),um(kdim-1),uu
    REAL :: xmet(kdim)   ! metrics term for diffusion
    REAL :: xmu_ave,xmut(kdim)
    INTEGER :: i,j,k,m,iflag
    REAL:: xc,yc,zc,tc
    REAL:: rhov, diff7(2),diff13(2),diff46(2),xma_re,sol(kdim-1)
    REAL:: rho_ave,tke_ave,ome_ave
    REAL :: rv_l, rv_r,rv_c,xl,yl,zl, xr,yr,zr,vleft,vright
    REAL :: xlc,xrc,xmu_l, xmu_r, tll, trr, tcc,amul,amur
    REAL :: re_xma
    
    xma_re = xmach/reue; re_xma =1./xma_re
    iflag = 0
    DO i=1,idim-1
       DO j=1,jdim-1
          DO k=1,kdim-1
             xc=0.5*(sk(j,k+1,i,1)*sk(j,k+1,i,4)+ &
                     sk(j,k,i  ,1)*sk(j,k,i  ,4))/vol(j,k,i)
             yc=0.5*(sk(j,k+1,i,2)*sk(j,k+1,i,4)+ &
                     sk(j,k,i  ,2)*sk(j,k,i  ,4))/vol(j,k,i)
             zc=0.5*(sk(j,k+1,i,3)*sk(j,k+1,i,4)+ &
                     sk(j,k,i  ,3)*sk(j,k,i  ,4))/vol(j,k,i)
             tc=0.5*(sk(j,k+1,i,5)*sk(j,k+1,i,4)+ &
                     sk(j,k,i  ,5)*sk(j,k,i  ,4))/vol(j,k,i)
             uu=xc*q(j,k,i,2)+yc*q(j,k,i,3)+zc*q(j,k,i,4)+tc

             up(k) = 0.5*(uu+abs(uu))
             um(k) = 0.5*(uu-abs(uu))
          ENDDO
          ! convection terms - 1st order(time and space)
          DO m=1,nummem
             DO k=1,kdim-1
                rlow(k,m)= -up(k)
                rup(k,m) = um(k)
                diag(k,m) =up(k)-um(k)
             ENDDO
          ENDDO

          DO k=1,kdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             IF(k==1) THEN
                rho_ave = 0.5*(qk0(j,i,1,1)+q(j,k,i,1))
             ELSEIF(k==kdim) THEN
                rho_ave  = 0.5*(qk0(j,i,1,3)+q(j,k-1,i,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k-1,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k-1,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k-1,i,7))
             xmut(k) = MAX(0.0,rho_ave*tke_ave/ome_ave)
          ENDDO
          DO k=1,kdim-1
             IF(k==1) THEN
                vleft = vk0(j,i,1)
             ELSE
                vleft = vol(j,k-1,i)
             ENDIF
             IF(k==kdim-1) THEN
                vright = vk0(j,i,3)
             ELSE
                vright = vol(j,k+1,i)
             ENDIF
             rv_l = 1./(0.5*(vleft + vol(j,k,i)))
             rv_r = 1./(0.5*(vright + vol(j,k,i)))
             rv_c = 1./vol(j,k,i)
             
             xl= (sk(j,k,i,1)*sk(j,k,i,4))*rv_l
             yl= (sk(j,k,i,2)*sk(j,k,i,4))*rv_l
             zl= (sk(j,k,i,3)*sk(j,k,i,4))*rv_l

             xr= (sk(j,k+1,i,1)*sk(j,k+1,i,4))*rv_r
             yr= (sk(j,k+1,i,2)*sk(j,k+1,i,4))*rv_r
             zr= (sk(j,k+1,i,3)*sk(j,k+1,i,4))*rv_r

             xc = 0.5*(sk(j,k,i,1)*sk(j,k,i,4)+sk(j,k+1,i,1)*sk(j,k+1,i,4))*rv_c
             yc = 0.5*(sk(j,k,i,2)*sk(j,k,i,4)+sk(j,k+1,i,2)*sk(j,k+1,i,4))*rv_c
             zc = 0.5*(sk(j,k,i,3)*sk(j,k,i,4)+sk(j,k+1,i,3)*sk(j,k+1,i,4))*rv_c
      
             xlc = xl*xc + yl*yc + zl*zc
             xrc = xr*xc + yr*yc + zr*zc
 
             amul=0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             amur=0.5*(fmu(j,k,i)+fmu(j,k+1,i))
             xmu_l = amul+xmut(k)*sigma
             xmu_r = amur+xmut(k+1)*sigma
             trr   = xmu_r*xrc
             tll   = xmu_l*xlc
             tcc   = trr+tll
            
             !for zeta equations
             !rhs(j,k,i,7) = rhs(j,k,i,7) + xma_re*(trr*turre(j,k+1,i,7) - tcc*turre(j,k,i,7) +tll*turre(j,k-1,i,7))/q(j,k,i,1)
             rlow(k,7) = rlow(k,7)-tll*xma_re/q(j,k,i,1)
             diag(k,7) = diag(k,7)+tcc*xma_re/q(j,k,i,1) 
             rup(k,7) = rup(k,7)  -trr*xma_re/q(j,k,i,1)

             xmu_l = amul+xmut(k)*sigma_star
             xmu_r = amur+xmut(k+1)*sigma_star
             trr   = xmu_r*xrc
             tll   = xmu_l*xlc
             tcc   = trr+tll
             !rhs(j,k,i,1:6) = rhs(j,k,i,1:6) + xma_re*(trr*turre(j,k+1,i,1:6) - tcc*turre(j,k,i,1:6) +tll*turre(j,k-1,i,1:6))/q(j,k,i,1)
             rlow(k,1:6) = rlow(k,1:6)-tll*xma_re/q(j,k,i,1)
             diag(k,1:6) = diag(k,1:6)+tcc*xma_re/q(j,k,i,1) 
             rup(k,1:6) = rup(k,1:6)  -trr*xma_re/q(j,k,i,1)
          ENDDO

          do k=1,kdim-1
             do m=1,3  ! no special treatment for the shear stresses(12,23,13)
                diag(k,m) = diag(k,m) +2./3.*beta_star*re_xma*turre(j,k,i,7) ! - 0.5*(source(j,k,i,m)+abs(source(j,k,i,m)))/turre(j,k,i,m)
             enddo
             diag(k,7) = diag(k,7) +2.*beta_0*re_xma*turre(j,k,i,7)!- 0.5*(source(j,k,i,7) - abs(source(j,k,i,7)))/turre(j,k,i,7)
          enddo 

          DO m=1,nummem
             DO k=1,kdim-1
                rlow(k,m) = rlow(k,m)*timestep(j,k,i)!source(j,k,i,m)
                rup(k,m)  = rup(k,m) *timestep(j,k,i)!source(j,k,i,m)
                diag(k,m) = diag(k,m)*timestep(j,k,i)+1.!source(j,k,i,m) + 1.
             ENDDO
          ENDDO
          DO m=1,nummem
             DO k=1,kdim-1
                rcol(k) = rhs(j,k,i,m)
             ENDDO
             CALL tri_solver(kdim-1,rlow(1,m),diag(1,m),rup(1,m),rcol,sol)  
             DO k=1,kdim-1
                rhs(j,k,i,m) = sol(k)
             ENDDO

          ENDDO
       ENDDO
    ENDDO
#ifdef DEBUG_KWSTM_AFSOLVER_K
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.or..false.) THEN
       write(*,*)"solver_k=", rhs(1:4,1,1,6)
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rdif","k-dif11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rdif","k-dif22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rdif","k-dif33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rdif","k-dif12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rdif","k-dif23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rdif","k-dif13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rdif","k-dif44.plt")
    ENDIF
#endif
  END SUBROUTINE afsolver_k

!----------- internal (private) subroutines
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE fill_turre(jdim,kdim,idim, nummem, turb,tj0,tk0,ti0,turre)
    INTEGER,INTENT(in) :: idim,jdim,kdim,nummem
    REAL,INTENT(out) :: turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem)
    REAL,INTENT(inout) :: turb(jdim,kdim,idim,nummem),tj0(kdim,idim-1,nummem,4),&
         tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
    COMMON /twod/ i2d
    integer :: i2d
    

    INTEGER :: i,j,k,iv

    turre(:,:,:,:)=0.
    turre(1:jdim-1, 1:kdim-1,1:idim-1, 1:nummem)= &
         turb(1:jdim-1, 1:kdim-1,1:idim-1, 1:nummem)
!   IF(i2d==1) THEN
!      DO iv = 1,nummem
!         DO j=1,jdim-1
!            DO k=1,kdim-1
!               ti0(j,k,iv,1:4) = turb(j,k,1,iv)
!            ENDDO
!         ENDDO
!      ENDDO
!   else
    
    DO j=1,jdim-1
       DO k=1,kdim-1
          DO iv =1,nummem
             turre(j,k,0, iv) =      ti0(j,k,iv,1)
             turre(j,k,-1,iv) =      ti0(j,k,iv,2)
             turre(j,k,idim,  iv) =  ti0(j,k,iv,3)
             turre(j,k,idim+1,iv) =  ti0(j,k,iv,4)
          ENDDO
       ENDDO
    ENDDO
!   ENDIF
    DO i=1,idim-1
       DO k=1,kdim-1
          DO iv =1,nummem
             turre( 0,k,i,iv) =      tj0(k,i,iv,1)
             turre(-1,k,i,iv) =      tj0(k,i,iv,2)
             turre(jdim,  k,i,iv) =  tj0(k,i,iv,3)
             turre(jdim+1,k,i,iv) =  tj0(k,i,iv,4)
          ENDDO
       ENDDO
    ENDDO

    DO i=1,idim-1
       DO j=1,jdim-1
          DO iv =1,nummem
             turre(j, 0,i,iv) =      tk0(j,i,iv,1)
             turre(j,-1,i,iv) =      tk0(j,i,iv,2)
             turre(j,kdim,  i,iv) =  tk0(j,i,iv,3)
             turre(j,kdim+1,i,iv) =  tk0(j,i,iv,4)
          ENDDO
       ENDDO
    ENDDO
!   if(i2d==1) then
!       turre(:,:,-1,:)=turre(:,:,1,:)
!       turre(:,:,0,:)=turre(:,:,1,:)
!       turre(:,:,2,:)=turre(:,:,1,:)
!       turre(:,:,3,:)=turre(:,:,1,:)
!   endif
  END SUBROUTINE fill_turre
  !
  ! calculate the turbulent kinetic energy
  !
  SUBROUTINE fill_tke(jdim,kdim,idim,nummem, turre,tke,bij)
    INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
    REAL,INTENT(in) :: turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem)
    REAL,INTENT(out) :: tke(0:jdim,0:kdim,0:idim),bij(0:jdim,0:kdim,0:idim,6)
    
    INTEGER :: i,j,k,m
    REAL :: coef
    
    DO i=0,idim
       DO k=0,kdim
          DO j=0,jdim
             tke(j,k,i) =  -0.5*(turre(j,k,i,1)+turre(j,k,i,2)+turre(j,k,i,3))
          ENDDO
       ENDDO
    ENDDO

    ! calculate bij = (cij + 2/3 k delta(ij))/k
    DO m=1,6
       coef=2./3.
       IF(m>3) coef = 0.
       DO i=1,idim-1
          DO k=1,kdim-1
             DO j=1,jdim-1
                bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    ! ghost cells
    DO m=1,6
       coef=2./3.
       if(m>3) coef = 0.
       DO i=1,idim-1
          DO k=1,kdim-1
             j=0
             bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
             j=jdim
             bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
          ENDDO
       ENDDO
       
       DO j=1,jdim-1
          DO k=1,kdim-1
             i=0
             bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
             i=idim
             bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
          ENDDO
       ENDDO

       DO j=1,jdim-1
          DO i=1,idim-1
             k=0
             bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
             k=kdim
             bij(j,k,i,m) = -0.5*(turre(j,k,i,m) +coef*tke(j,k,i))/tke(j,k,i)
          ENDDO
       ENDDO
       
    ENDDO
    RETURN
  END SUBROUTINE fill_tke

  ! calculate the vorticity
  SUBROUTINE fill_omega(jdim,kdim,idim, ux,omega)
    INTEGER,intent(in) :: jdim,kdim,idim
    REAL,INTENT(in) :: ux(:,:,:,:)
    REAL,INTENT(out) :: omega(0:jdim,0:kdim,0:idim,3)

    INTEGER :: i,j,k

    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             omega(j,k,i,1) = ux(j,k,i,8)-ux(j,k,i,6)  !(wy-vz)
             omega(j,k,i,2) = ux(j,k,i,3)-ux(j,k,i,7)  !(uz-wx)
             omega(j,k,i,3) = ux(j,k,i,4)-ux(j,k,i,2)  !(vx-uy)
          ENDDO
       ENDDO
    ENDDO
    RETURN
  END SUBROUTINE fill_omega

  ! calculate the laminar dynamic viscosity
  SUBROUTINE fill_fmu(jdim,kdim,idim,q,qj0,qk0,qi0,fmu)
    INTEGER, INTENT(in) :: jdim,kdim,idim
    REAL, INTENT(in) :: q(jdim,kdim,idim,5),&
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
    REAL, INTENT(out) :: fmu(0:jdim,0:kdim,0:idim)

    COMMON /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
    REAL ::  gamma,gm1,gp1,gm1g,gp1g,ggm1
    
    COMMON /fluid2/ pr,prt,cbar
    REAL :: pr,prt,cbar

    COMMON /reyue/ reue,tinf,ivisc(3)
    REAL :: reue,tinf
    integer :: ivisc

    INTEGER :: i,j,k
    REAL :: tt,c2b,c2bp

    c2b=cbar/tinf
    c2bp=c2b+1.0

    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             tt=gamma*q(j,k,i,5)/q(j,k,i,1)
             fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
          ENDDO
       ENDDO
    ENDDO

    DO k=1,kdim-1
       DO j=1,jdim-1
          i=0
          tt=gamma*qi0(j,k,5,1)/qi0(j,k,1,1)
          fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
          
          i=idim
          tt=gamma*qi0(j,k,5,3)/qi0(j,k,1,3)
          fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
       ENDDO
    ENDDO
    
    DO j=1,jdim-1
       DO i=1,idim-1
          k=0
          tt=gamma*qk0(j,i,5,1)/qk0(j,i,1,1)
          fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
          k=kdim
          tt=gamma*qk0(j,i,5,3)/qk0(j,i,1,3)
          fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
       ENDDO
    ENDDO

    DO i=1,idim-1
       DO k=1,kdim-1
          j=0
          tt=gamma*qj0(k,i,5,1)/qj0(k,i,1,1)
          fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
          j=jdim
          tt=gamma*qj0(k,i,5,3)/qj0(k,i,1,3)
          fmu(j,k,i)=c2bp*tt*SQRT(tt)/(c2b+tt)
       ENDDO
    ENDDO
    return
  END SUBROUTINE fill_fmu


  SUBROUTINE get_timestep(jdim,kdim,idim,dtj,vol,timestep,icyc)
    INTEGER,INTENT(in) :: jdim,kdim,idim,icyc
    REAL,INTENT(in) :: dtj(jdim,kdim,idim-1),vol(jdim,kdim,idim-1)
    REAL,INTENT(out) :: timestep(jdim-1,kdim-1,idim-1)

    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
    REAL ::  title,rkap,xmach,alpha__,beta__,dt,fmax
    INTEGER :: nit,ntt, &
         idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
         iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv


    COMMON /turbconv/ cflturb(7)
    real :: cflturb

    COMMON /unst/ time,cfltau,ntstep,ita,iunst
    REAL :: time, cfltau
    INTEGER ::ntstep,ita,iunst


    INTEGER :: i,j,k
    integer,save :: iflag=0
    REAL :: rmax,rmin

    factor1=factor
    factor2=factor
    !
    ! Overwrite factors with keyword value "cflturb()" if nonzero
    !
    IF (REAL(cflturb(1)).NE.0.) THEN
       factor1 = cflturb(1)
    END IF
    IF (REAL(cflturb(2)).NE.0.) THEN
       factor2 = cflturb(2)
    END IF
    !factor2 is set relative to factor1
    factor2=factor2/factor1
    !
    ! Timestep for turb model
    !
     
    IF (REAL(dt).LT.0) THEN
       DO i=1,idim-1
          DO k=1,kdim-1
             DO j=1,jdim-1
                timestep(j,k,i)=factor1*vol(j,k,i)/dtj(j,k,i)
                timestep(j,k,i)=min(timestep(j,k,i),100.)
             ENDDO
          ENDDO
       ENDDO
       iflag =1  
    ELSE
       !turbulence model advanced WITH physical time ONLY
       !(pseudo-time term NOT included, even for tau-TS in mean-
       !flow equations, since multigrid is not used for turb. eq.)
       IF(ita>0) THEN
          DO i=1,idim-1
             DO k=1,kdim-1
                DO j=1,jdim-1
                   timestep(j,k,i)=dt
                   factor2=1.
                ENDDO
             ENDDO
          ENDDO
       ELSE
          DO i=1,idim-1
             DO k=1,kdim-1
                DO j=1,jdim-1
                   timestep(j,k,i)=factor1*vol(j,k,i)/dtj(j,k,i)
                   timestep(j,k,i)=MIN(timestep(j,k,i),100.)
                ENDDO
             ENDDO
          ENDDO
       END IF
       IF(MOD(ntt,50)==0.AND.icyc==1) THEN
          rmax = 0
          rmin = huge(rmin)
          DO i=1,idim-1
             DO k=1,kdim-1
                DO j=1,jdim-1
                   rmax=max(dt/timestep(j,k,i),rmax)
                   rmin=min(dt/timestep(j,k,i),rmin)
                ENDDO
             ENDDO
          ENDDO
          
       ENDIF
    ENDIF
  END SUBROUTINE get_timestep

  SUBROUTINE get_advection(jdim,kdim,idim,nummem,q,turre,sj,sk,si,vol,qj0,qk0,qi0,rhs,d,al,ar,bl,br)
    INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5),&
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),&
         sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1), &
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4)
    REAL,INTENT(inout):: rhs(jdim-1,kdim-1,idim-1,nummem),& 
         d(nummem,nummem,jdim,kdim,idim)
    REAL, INTENT(inout),DIMENSION(2,jdim,kdim,idim) :: al,ar,bl,br
         

    COMMON /twod/ i2d
    INTEGER :: i2d
    COMMON /turbconv/ cflturb(7),edvislim,iturbprod,nsubturb,nfreeze, &
         iwarneddy,itime2read,itaturb,tur1cut,tur2cut,& 
         iturbord,tur1cutlev,tur2cutlev
    REAL :: cflturb,edvislim,tur1cut,tur2cut, tur1cutlev,tur2cutlev
    INTEGER :: iturbprod,nsubturb,nfreeze, iwarneddy,itime2read,itaturb,iturbord


    REAL,SAVE :: order2=0.0  ! =1.0, second order scheme
    REAL :: a, b, c
    integer :: i,j,k,m
    REAL :: xc,yc,zc,tc,sgnu,app,apm,uu
    ! vel: cell interface velocity
    REAL :: order,flux(nummem, MAX(idim,jdim,kdim)),vel(3,MAX(idim,kdim,jdim)), &
         qc(5,0:MAX(jdim,kdim,idim)),ul(MAX(idim,jdim,kdim)),ur(MAX(idim,jdim,kdim)),&
         lamb(MAX(jdim,kdim,idim))

    
    order2 = 0
    if(iturbord==2) order2=1
    a = 1.0 + order2*0.5
    b = 1.0 + order2
    c = order2*0.5
    ! j-direction
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim-1
             qc(1:5,j) = q(j,k,i,1:5)
          ENDDO
          qc(1:5,0) = qj0(k,i,1:5,1)
          qc(1:5,jdim) = qj0(k,i,1:5,3)
          ! -----------------------------
          DO j=1,jdim
             UL(j) = (qc(2,j-1)*sj(j,k,i,1)+qc(3,j-1)*sj(j,k,i,2)+qc(4,j-1)*sj(j,k,i,3))*sj(j,k,i,4)
             UR(j) = (qc(2,j)*sj(j,k,i,1)+qc(3,j)*sj(j,k,i,2)+qc(4,j)*sj(j,k,i,3))*sj(j,k,i,4)
             !Rusanov Scheme
             lamb(j) = MAX(ABS(ur(j)),ABS(ul(j)))
             flux(:,j) = 0.5*(turre(j-1,k,i,:)*qc(1,j-1)*UL(j) +turre(j,k,i,:)*qc(1,j)*UR(j)-&
                  lamb(j)*(turre(j,k,i,:)*qc(1,j) - turre(j-1,k,i,:)*qc(1,j-1)))
          ENDDO

          IF(turre(-1,k,i,1)==HUGE(turre(-1,k,i,1))) THEN
             flux(:,1) = 0. 
          ENDIF
          
          IF(turre(jdim+1,k,i,1)==HUGE(turre(jdim+1,k,i,1))) THEN
             flux(:,jdim) = 0.
          ENDIF
             
          DO j=1,jdim-1
             al(:,j,k,i) = al(:,j,k,i) - 0.5*(ul(j)+lamb(j))/vol(j,k,i)
             ar(:,j,k,i) = ar(:,j,k,i) + 0.5*(ur(j+1)-lamb(j+1))/vol(j,k,i)
             DO m=1,nummem
                rhs(j,k,i,m)=rhs(j,k,i,m)- (flux(m,j+1)-flux(m,j))/(q(j,k,i,1)*vol(j,k,i))
                d(m,m,j,k,i) = d(m,m,j,k,i) - 0.5*(ul(j+1) - ur(j)+lamb(j)+lamb(j+1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    !k-direction
    DO i=1,idim-1
       DO j=1,jdim-1
          DO k=1,kdim-1
             qc(1:5,k) = q(j,k,i,1:5)
          ENDDO
          qc(1:5,0) = qk0(j,i,1:5,1)
          qc(1:5,kdim) = qk0(j,i,1:5,3)
          ! -----------------------------
          DO k=1,kdim
             UL(k) = (qc(2,k-1)*sk(j,k,i,1)+qc(3,k-1)*sk(j,k,i,2)+qc(4,k-1)*sk(j,k,i,3))*sk(j,k,i,4)
             UR(k) = (qc(2,k)*sk(j,k,i,1)+qc(3,k)*sk(j,k,i,2)+qc(4,k)*sk(j,k,i,3))*sk(j,k,i,4)
             !Rusanov Scheme
             lamb(k) = MAX(ABS(ur(k)),ABS(ul(k)))
             flux(:,k) = 0.5*(turre(j,k-1,i,:)*qc(1,k-1)*UL(k) +turre(j,k,i,:)*qc(1,k)*UR(k)-&
                  lamb(k)*(turre(j,k,i,:)*qc(1,k) - turre(j,k-1,i,:)*qc(1,k-1)))
          ENDDO
          IF(turre(j,-1,i,1)==HUGE(turre(j,-1,i,1))) THEN
             flux(:,1) = 0.
          ENDIF
          
          IF(turre(j,kdim+1,i,1)==HUGE(turre(j,-1,i,1))) THEN
             flux(:,kdim) = 0.
          ENDIF
          
          DO k=1,kdim-1
             bl(:,j,k,i) = bl(:,j,k,i) - 0.5*(ul(k)+lamb(k))/vol(j,k,i)
             br(:,j,k,i) = br(:,j,k,i) + 0.5*(ur(k+1)-lamb(k+1))/vol(j,k,i)
             DO m=1,nummem
                rhs(j,k,i,m)=rhs(j,k,i,m)-(flux(m,k+1)-flux(m,k))/(q(j,k,i,1)*vol(j,k,i))
                d(m,m,j,k,i) = d(m,m,j,k,i) - 0.5*(ul(k+1)-ur(k)+lamb(k)+lamb(k+1))/vol(j,k,i)
             ENDDO
          END DO
       ENDDO
    ENDDO

    IF(i2d==0) THEN
       !i-direction
       DO j=1,jdim-1
          DO k=1,kdim-1
             DO i=1,kdim-1
                qc(1:5,i) = q(j,k,i,1:5)
             ENDDO
             qc(1:5,0) = qi0(j,k,1:5,1)
             qc(1:5,idim) = qi0(j,k,1:5,3)
             ! -----------------------------
             DO i=1,idim
                UL(i) = (qc(2,i-1)*si(j,k,i,1)+qc(3,i-1)*si(j,k,i,2)+qc(4,k-1)*si(j,k,i,3))*si(j,k,i,4)
                UR(i) = (qc(2,i)*si(j,k,i,1)+qc(3,i)*si(j,k,i,2)+qc(4,k)*si(j,k,i,3))*si(j,k,i,4)
                !Rusanov Scheme
                flux(:,i) = 0.5*(turre(j,k,i-1,:)*qc(1,i-1)*UL(i) +turre(j,k,i,:)*qc(1,i)*UR(i)-&
                     MAX(ABS(ur(i)),ABS(ul(i)))*(turre(j,k,i,:)*qc(1,i) - turre(j,k,i-1,:)*qc(1,i-1)))
             ENDDO

             DO i=1,idim-1
                do m=1,nummem
                   rhs(j,k,i,m)=rhs(j,k,i,m)-(flux(m,i+1)-flux(m,i))/(q(j,k,i,1)*vol(j,k,i))
                   d(m,m,j,k,i) = d(m,m,j,k,i) - 0.25*(ABS(ul(i))+ABS(ul(i+1))+ABS(ur(i))+ABS(ur(i+1)))/(q(j,k,i,1)*vol(j,k,i))
                enddo
             END DO
          ENDDO
       ENDDO
    ENDIF
#ifdef DEBUG_KWSTM_ADVECTION
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.or..false.) THEN
       write(*,*)"adv=",rhs(1,kdim-3:kdim-1,1,1)
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "radv","r-adv11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "radv","r-adv22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "radv","r-adv33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "radv","r-adv12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "radv","r-adv23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "radv","r-adv13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "radv","r-adv44.plt")
    ENDIF
#endif
  END SUBROUTINE get_advection

  ! use fractional approach to update the source term
  SUBROUTINE fractional(jdim,kdim,idim,nummem,rhs,turre,timestep,d,zksav) 
    INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
    REAL,INTENT(in) :: timestep(jdim-1,kdim-1,idim-1),zksav(jdim,kdim,idim,nummem)
    REAL, intent(inout) :: rhs(jdim-1,kdim-1,idim-1,nummem),d(nummem,nummem, jdim,kdim,idim),&
        turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem)
    integer :: i,j,k,m
    real :: v1(nummem),v2(nummem)
    
    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax
    REAL :: title,rkap,xmach,alpha__,beta__,dt,fmax
    COMMON /unst/ time,cfltau,ntstep,ita,iunst
    REAL :: time, cfltau
    INTEGER ::ntstep,ita,iunst
   
   
    DO m=1,nummem
       DO i=1,idim-1
          DO k=1,kdim-1
             DO j=1,jdim-1
                rhs(j,k,i,m) = rhs(j,k,i,m)*timestep(j,k,i)
                d(:,m,j,k,i) = -d(:,m,j,k,i)*timestep(j,k,i)
                d(m,m,j,k,i) = 1+ d(m,m,j,k,i)
             ENDDO
          enddo
       enddo
    enddo
   

    !diagonal solver
    
    DO i=1,idim-1
       do k=1,kdim-1
          do j=1,jdim-1
             CALL lu_ijk(d(1,1,j,k,i),7)
             v1(1:7) = rhs(j,k,i,1:7)
             call lu_solv(d(1,1,j,k,i),v1,7,v2)
             turre(j,k,i,:) = turre(j,k,i,:)+v2(1:7)
             rhs(j,k,i,1:7) = v2(1:7)
          enddo
       enddo
    enddo

#ifdef DEBUG_KWSTM_FRACTIONAL
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.OR..false.) THEN
       write(*,*)"aftsrc=",rhs(1,kdim-3:kdim-1,1,1)
       write(*,*)"zeta=", rhs(1,1:3,1,7)
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rhs","b-rhs11.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rhs","b-rhs22.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rhs","b-rhs33.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rhs","b-rhs12.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rhs","b-rhs23.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rhs","b-rhs13.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rhs","b-rhs44.plt")
    ENDIF
#endif
    rhs = 0.
    return
  END SUBROUTINE fractional


  SUBROUTINE preprocess_rhs(jdim,kdim,idim,nummem,q,turre,timestep,source, rhs,d,iopt,zksav)
    INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
    REAL,INTENT(in) :: q(jdim,kdim,idim,5),turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),&
         timestep(jdim-1,kdim-1,idim-1),zksav(jdim,kdim,idim,nummem)
    REAL, intent(inout) :: source(jdim-1,kdim-1,idim-1,nummem),&
         rhs(jdim-1,kdim-1,idim-1,nummem),d(nummem,nummem, jdim,kdim,idim)
    integer, intent(in) :: iopt

    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10

    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax
    real :: title,rkap,xmach,alpha__,beta__,dt,fmax
    COMMON /unst/ time,cfltau,ntstep,ita,iunst
    REAL :: time, cfltau
    INTEGER ::ntstep,ita,iunst

    INTEGER :: i,j,k,m
    REAL :: dsdcij,dsdcijmax(nummem),cmm,v1(7),v2(7)
    ! rhs term times time step
    ! RHS = RHS / ( 1/dt  - dSource/dCij)
    dsdcijmax =0
    IF(iopt==0) THEN
       DO m=1,nummem
          DO i=1,idim-1
             DO k=1,kdim-1
                DO j=1,jdim-1
                   d(:,m,j,k,i) = -d(:,m,j,k,i)*timestep(j,k,i)
                   d(m,m,j,k,i) = 1+ d(m,m,j,k,i)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
    endif

    DO m=1,nummem
       DO i=1,idim-1
          DO k=1,kdim-1
             DO j=1,jdim-1
                rhs(j,k,i,m) = rhs(j,k,i,m)*timestep(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    
    IF(iopt==0) THEN
#ifdef DEBUG_KWSTM_PREPROCESS
       IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.OR..false.) THEN
          WRITE(*,'(A,7(1x,es12.5))')"bf_source=",((d(1,m,1,k,1),m=1,7),k=1,kdim-1)
          WRITE(*,*)"zeta=", rhs(1,1:3,1,7)
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rhs","a-rhs11.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rhs","a-rhs22.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rhs","a-rhs33.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rhs","a-rhs12.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rhs","a-rhs23.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rhs","a-rhs13.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rhs","a-rhs44.plt")

       ENDIF
#endif
       !diagonal solver
       DO i=1,idim-1
          DO k=1,kdim-1
             DO j=1,jdim-1
                CALL lu_ijk(d(1,1,j,k,i),7)
                v1(1:7) = rhs(j,k,i,1:7)
                CALL lu_solv(d(1,1,j,k,i),v1,7,v2)
                rhs(j,k,i,1:7) = v2(1:7)
             ENDDO
          ENDDO
       ENDDO

#ifdef DEBUG_KWSTM_PREPROCESS
       IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.OR..false.) THEN
          WRITE(*,*)"aftsrc=",rhs(1,kdim-3:kdim-1,1,1)
          WRITE(*,*)"zeta=", rhs(1,1:3,1,7)
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rhs","b-rhs11.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rhs","b-rhs22.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rhs","b-rhs33.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rhs","b-rhs12.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rhs","b-rhs23.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rhs","b-rhs13.plt")
          CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rhs","b-rhs44.plt")
       ENDIF
#endif
    ENDIF
    !write(*,'(A,7(1x,es12.5))')"dsdcijmax=",dsdcijmax
    RETURN
  END SUBROUTINE preprocess_rhs

  SUBROUTINE update(jdim,kdim,idim,nummem,rhs,turre,rho,fmu, turb,vist3d,sumn,negn,ux)
    INTEGER,INTENT(in):: jdim,kdim,idim,nummem
    REAL, INTENT(in) :: rhs(jdim-1,kdim-1,idim-1,nummem),&
         rho(jdim,kdim,idim),fmu(0:jdim,0:kdim,0:idim),ux(jdim-1,kdim-1,idim-1,9)
    REAL,INTENT(out) :: turb(jdim,kdim,idim,nummem),&
         vist3d(jdim,kdim,idim)
    REAL,INTENT(out) :: sumn(nummem)
    REAL,INTENT(inout):: turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem)
    INTEGER,intent(out) :: negn(nummem)
    
    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
    REAL ::  title,rkap,xmach,alpha__,beta__,dt,fmax
    INTEGER :: nit,ntt, &
         idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
         iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv
    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10
    INTEGER :: i,j,k,m
    integer, save ::imax,jmax,kmax
    REAL :: tke,zeta_min,t1122,t2233,t1133
    REAL,save :: vtmax =0.0,tij_max(7)
    REAL, PARAMETER :: mut_top = 1e4
    INTEGER,SAVE :: counter=0,nlast =-1
    real :: s12,s23,s13

    COMMON /turbconv/ cflturb(7),edvislim
    REAL :: cflturb,edvislim
    COMMON /mydist2/ nnodes,myhost,myid,mycomm
    INTEGER::nnodes,myhost,myid,mycomm
    !return
    sumn = 0
    negn = 0
    counter = counter + 1
    IF(mylevel==level_o.AND.MOD(myicyc,icyc_o)==0.or..false.) THEN
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,1), "rhs","r1.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,2), "rhs","r2.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,3), "rhs","r3.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,4), "rhs","r4.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,5), "rhs","r5.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,6), "rhs","r6.plt")
       CALL dump_tecplot(jdim-1, kdim-1, idim-1, 1,jdim-1,1,kdim-1,1, rhs(1,1,1,7), "rhs","r7.plt")
    ENDIF
    DO m=1,nummem;
       DO i=1,idim-1; 
          DO k=1,kdim-1;
             DO j=1,jdim-1
                turb(j,k,i,m) = turre(j,k,i,m)+rhs(j,k,i,m)
                ! force turb(j,k,i,1:3) to be negative
                IF(m<=3) THEN
                   IF(turb(j,k,i,m)>0.) THEN
                      negn(m) = negn(m)+1
                   ENDIF
                   turb(j,k,i,m)=MIN(turb(j,k,i,m),tur10(m)*1e-7)
                ENDIF
                if(m==4.and.irealizability==1) then
                   t1122= sqrt(turb(j,k,i,1)*turb(j,k,i,2))
                   turb(j,k,i,m) = min(t1122, max(-t1122,turb(j,k,i,m)))
                   !s12=ux(j,k,i,2)+ux(j,k,i,4)
                   !if(turb(j,k,i,m)*s12<0) turb(j,k,i,m) = 0
                endif
                if(m==5.and.irealizability==1) then
                   t2233= sqrt(turb(j,k,i,2)*turb(j,k,i,3))
                   turb(j,k,i,m) = min(t2233, max(-t2233,turb(j,k,i,m)))
                   !s23=ux(j,k,i,6)+ux(j,k,i,8)
                   !if(turb(j,k,i,m)*s23<0) turb(j,k,i,m) = 0
                endif
                if(m==6.and.irealizability==1) then
                   t1133= sqrt(turb(j,k,i,1)*turb(j,k,i,3))
                   turb(j,k,i,m) = min(t1133, max(-t1133,turb(j,k,i,m)))
                   !s13 = ux(j,k,i,3)+ux(j,k,i,7) 
                   !if(turb(j,k,i,m)*s13<0) turb(j,k,i,m) = 0
                endif
                IF(m==7) THEN
                   IF(turb(j,k,i,m)<0) THEN
                      negn(m) = negn(m)+1
                   ENDIF
                   turb(j,k,i,m)=MAX(turb(j,k,i,m),tur10(7)*1e-7)
                ENDIF
                sumn(m) = sumn(m) + (turb(j,k,i,m)-turre(j,k,i,m))**2
                turre(j,k,i,m) = turb(j,k,i,m)
                tke =-0.5*(turb(j,k,i,1)+turb(j,k,i,2)+turb(j,k,i,3))
                vist3d(j,k,i) = MIN(edvislim,rho(j,k,i)*tke/turre(j,k,i,7))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    DO m=1,nummem
       sumn(m) = SQRT(sumn(m))/float((jdim-1)*(kdim-1)*(idim-1))
    ENDDO

#ifdef DEBUG_KWSTM
    ! get viscosity
    IF(MOD(ntt,50)==0) THEN
       WRITE(200+myid,'(A,10(1x,es12.5))')"summ=",sumn(1:nummem)
       if(vtmax/=0.0) then
          WRITE(200+myid,'(A,1x,I8,18(1x,es12.5))')"vtmax=",ntt,vtmax,tij_max(1:7),dt,cfl_psd,cfl_loc
       endif
       vtmax = 0
    endif
    IF(MOD(ntt,50)==49) THEN
       DO i=1,idim-1; 
          DO k=1,kdim-1;
             DO j=1,jdim-1
                IF(vist3d(j,k,i)>vtmax) THEN
                   vtmax = MAX(vtmax,vist3d(j,k,i))
                   tij_max = turre(j,k,i,1:7)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
    ENDIF
#endif
  END SUBROUTINE update

  ! tridiagonal solver
  ! input:
  !   tridiagonal system
  !    | d1 u1  ......|
  !    ! l2 d2 u2.....|
  !    ! 0  l3 d3 u3..|
  !    |   ...        |
  !    |         ln dn|
  SUBROUTINE tri_solver(nsize,rlow,d,rup,rhs,sol)
    INTEGER,intent(in) :: nsize
    REAL, INTENT(in) :: rlow(nsize),rup(nsize),rhs(nsize)
    REAL, intent(inout) :: d(nsize)
    REAL, intent(out) :: sol(nsize)

    INTEGER :: i
    real :: sum

    ! factorization
    DO i=2,nsize
       d(i) = d(i) - rlow(i)*rup(i-1)/d(i-1)
    ENDDO
    
    ! forward sweep
    sol(1) = rhs(1)/d(1)
    DO i=2,nsize
       sol(i) = (rhs(i) - rlow(i)*sol(i-1))/d(i)
    ENDDO
     
    DO i=1,nsize-1
       sol(i) = sol(i)*d(i)
    ENDDO
    
    DO i=nsize-1,1,-1
       sol(i) = (sol(i) - rup(i)*sol(i+1))/d(i)
    ENDDO
    
    RETURN
  END SUBROUTINE tri_solver

  ! tecplot file generation subroutine
  !    - dump j-k slice of a 3D array
  
  SUBROUTINE dump_tecplot(jdim, kdim, idim, jstr,jend,kstr,kend,iplane, u, vname,fname)
    INTEGER,INTENT(in) :: jdim, kdim, idim,jstr,jend,kstr,kend,iplane
    REAL, INTENT(in) :: u(jdim,kdim,idim)
    CHARACTER(len=*), INTENT(in) :: vname, fname

    INTEGER :: j,k
    
    OPEN(1234,file=TRIM(ADJUSTL(fname)))
    WRITE(1234,*)'Title="', TRIM(ADJUSTL(vname)),'"'
    WRITE(1234,*)'Variables="I", "J","', TRIM(ADJUSTL(vname)),'"'
    WRITE(1234,*)'Zone I=',jend-jstr+1, " J=", kend-kstr+1
    DO k=kstr,kend
       DO j=jstr,jend
          WRITE(1234,"(2(1x,I8),1x,es18.5)")j,k, real(u(j,k,iplane),4)
       ENDDO
    ENDDO
    close(1234)
  END SUBROUTINE dump_tecplot
          
  SUBROUTINE dump_tecplot_turb(jdim, kdim, idim,nummem, jstr,jend,kstr,kend,iplane, u,vist,fname)
    INTEGER,INTENT(in) :: jdim, kdim, idim,jstr,jend,kstr,kend,iplane,nummem
    REAL, INTENT(in) :: u(jdim,kdim,idim,nummem),vist(jdim,kdim,idim)
    CHARACTER(len=*), INTENT(in) :: fname
    
    INTEGER :: j,k
    
    OPEN(1234,file=TRIM(ADJUSTL(fname)))
    WRITE(1234,*)'Title="', TRIM(ADJUSTL(fname)),'"'
    WRITE(1234,*)'Variables="I", "J", "t11", "t22", "t33", "t12", "t23", "t13","zeta","mut"'
    WRITE(1234,*)'Zone I=',jend-jstr+1, " J=", kend-kstr+1
    DO k=kstr,kend
       DO j=jstr,jend
          WRITE(1234,"(2(1x,I8),8(1x,es12.5))")j,k, u(j,k,iplane,1:7),vist(j,k,iplane)
       ENDDO
    ENDDO
    CLOSE(1234)
  END SUBROUTINE dump_tecplot_turb
  !
  !only dump at the finest level of the grid
  !   level.ge.lglobal
  !
  SUBROUTINE  kws_plot( &
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
    REAL :: c2b,c2bp, fnu(jdim),tt,re,utau(jdim),vel,tau_w,fnuloc,Rel(jdim),cf(jdim)
    REAL :: deltabl, coef
    COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
         tur10(7)
    REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
         tur10,vel1
    

    if(.not.need_kwsplot) return

    WRITE(str,'(I4.4)')nblk
#ifndef DIST_MPI
    IF(i2d==1) THEN
       OPEN(1235,file="yplus-"//TRIM(ADJUSTL(str))//".plt")
       OPEN(1236, file="Cf-ReL-"//trim(adjustl(str))//".plt")
       WRITE(1236,'(A)')'variables=rel,cf,cp,x'
       WRITE(1236,'(A,1x,I4)')'ZONE I=',jdim-1
       i=1
       re=reue/xmach
!      WRITE(11,*)"Generate Y+ vs U+ profile for block:",nblk
       DO j=1,jdim-1
          c2b=cbar/tinf
          c2bp=c2b+1.0
          tt=gamma*q(j,1,i,5)/q(j,1,i,1)
          fnu(j)=c2bp*tt*SQRT(tt)/(c2b+tt)
          vel = sign(SQRT(q(j,1,i,2)**2+q(j,1,i,3)**2+q(j,1,i,4)**2),1.)!q(j,1,i,2))
          tau_w = fnu(j)*abs(vel)/ABS(smin(j,1,i))/Re
          utau(j) = SQRT(tau_w/q(j,1,i,1))
          Rel(j) = (0.5*(x(j,1,i)+x(j+1,1,i)) - x(jstart,1,i))*Reue
          write(1236,'(4(1x,es12.5))') rel(j),tau_w/(0.5*q(j,kdim-1,i,1)*xmach**2),&
               (q(j,1,i,5)-p0)/(0.5*rho0*xmach**2), 0.5*(x(j,1,i)+x(j+1,1,i))
       ENDDO
       close(1236) 
       j = jcut
       i = 1
       DO k=1,kdim-2
          vel = SQRT(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2)
          vel1 = SQRT(q(j,k+1,i,2)**2+q(j,k+1,i,3)**2+q(j,k+1,i,4)**2)
          IF(vel<0.99*u0.AND.vel1>0.99*u0) THEN
             coef = (0.99*u0 - q(j,k,i,2))/(q(j,k+1,i,2)-q(j,k,i,2))
             deltabl = smin(j,k,i)+coef*(smin(j,k+1,i)-smin(j,k,i))
          ENDIF
       ENDDO
       WRITE(plt_title,'(5(1x,A,es12.5))')'utau=',utau(j),'fnuw=',fnu(j),'Re=',reue,'Mach=',xmach,&
            "delta=",deltabl
       WRITE(1235,'(A)')'title="'//TRIM(ADJUSTL(plt_title))//'"'
       IF(nummem==2) THEN
          WRITE(1235,'(A)')'variables="y+", "u+", "y","rho", "u", "v", "w", "p", "ome", "tke", "mut","dist"'
       ELSEIF(nummem==7) THEN
          WRITE(1235,'(A)')'variables="y+", "u+", "y","rho", "u", "v", "w", "p", "t11", "t22",'//&
               '"t33", "t12", "t23", "t13", "zeta", "mut","dist"'
       elseif(nummem==0) then
          WRITE(1235,'(A)')'variables="y+", "u+", "y","rho", "u", "v", "w", "p", "t11", "t22",'//&
               '"t33", "t12", "t23", "t13", "zeta", "mut","dist"'
       ENDIF
       WRITE(1235,'(A,1x,I4)')'ZONE I=',kdim-1
       ! search for the 0.99 v_inf
       DO k=1,kdim-1
          vel = SQRT(q(j,k,i,2)**2+q(j,k,i,3)**2+q(j,k,i,4)**2)
          if(nummem/=0) then
             WRITE(1235,'(20(1x,es12.5))') &
                  q(j,1,i,1)*ABS(smin(j,k,i))*utau(j)*re/fnu(j), &
                  vel/utau(j),smin(j,k,i)/deltabl,q(j,k,i,1:5),turre(j,k,i,1:nummem),vist3d(j,k,i),smin(j,k,i)
          else
             WRITE(1235,'(20(1x,es12.5))') &
                  q(j,1,i,1)*ABS(smin(j,k,i))*utau(j)*re/fnu(j), &
                  vel/utau(j),smin(j,k,i)/deltabl,q(j,k,i,1:5),(0.0, m=1,8),smin(j,k,i)
          endif
       ENDDO
       close(1235)
    endif  ! finish 2d if-block
#endif
    RETURN

  END SUBROUTINE kws_plot

  ! fix the energy equation by adding the src_k to the residual vector
  subroutine kws_fix_eng(jdim,kdim,idim, res,vol)
    integer, intent(in) :: jdim,kdim,idim
    real, intent(inout) :: res(jdim,kdim,idim-1,5)
    real, intent(in) :: vol(jdim,kdim,idim-1)
    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
    REAL ::  title,rkap,xmach,alpha__,beta__,dt,fmax
    INTEGER :: nit,ntt, &
         idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
         iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv

    integer :: i,j,k
    do i=1,idim-1
       do k=1,kdim-1
          do j=1,jdim-1
             res(j,k,i,5)  = res(j,k,i,5)+src_k(j,k,i)*vol(j,k,i)
          enddo
       enddo
    enddo
  end subroutine

   ! lu decomposition
   SUBROUTINE lu_ijk(a,n)
     INTEGER,INTENT(in) ::n
     REAL, INTENT(inout) :: a(n,n)
 
     INTEGER:: i,j,k
 
     DO i=1, n
        DO j=2,i
           a(i,j-1) = a(i,j-1)/a(j-1,j-1)
           DO k=1,j-1
              a(i,j) = a(i,j) - a(i,k)*a(k,j)
           ENDDO
        ENDDO
 
        DO j=i+1,n
           DO k=1,i-1
              a(i,j) = a(i,j)-a(i,k)*a(k,j)
           ENDDO
        ENDDO
     ENDDO
 
   END SUBROUTINE lu_ijk
 
   ! L(i,i) = 1.0 
   SUBROUTINE lu_solv(a,rhs,n,sol)
     INTEGER,INTENT(in) :: n
     REAL, INTENT(in) :: a(n,n)
     REAL, intent(in) :: rhs(n)
     REAL, intent(out) :: sol(n)
 
     REAL :: sum
     INTEGER :: i,j
 
     ! Lx = rhs
     DO i=1,n
        sum = 0
        DO j=1,i-1
           sum = sum + a(i,j)*sol(j)
        ENDDO
        sol(i) = rhs(i) - sum
     ENDDO
     ! u y = x
     DO i=n,1,-1
        sum = 0
        DO j=n,i+1,-1
           sum = sum+a(i,j)*sol(j)
        ENDDO
        sol(i) = (sol(i) - sum)/a(i,i)
     ENDDO
     return
   END SUBROUTINE lu_solv

   !  get terms needed for time-accurate calculation 
   !  This subroutine performed:
   !      1) dq = qcur -qold
   !      2) qold <= qcur
   SUBROUTINE save_lasttimestep(jdim,kdim,idim, nummem, qcur, qold,dq)
     INTEGER,INTENT(in) :: jdim,kdim,idim,nummem
     REAL, INTENT(in) :: qcur(jdim,kdim,idim,nummem)
     REAL, INTENT(inout)  :: qold(jdim,kdim,idim,nummem)
     REAL, INTENT(out)    :: dq(jdim,kdim,idim,nummem)

     COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax
     REAL :: title,rkap,xmach,alpha__,beta__,dt,fmax

     INTEGER :: i

     IF(dt<=0) RETURN
     
     IF(qold(1,1,1,1)==0.0) THEN
        dq = 0
     ELSE
        dq = qcur - qold
     ENDIF
     
     qold = qcur

   END SUBROUTINE save_lasttimestep

   ! a : for j-direction (advection and diffusion)
   ! b : for k-direction (advection and diffusion)
   
   SUBROUTINE sgs_solver_2d(jdim,kdim,idim,nummem,timestep, d,al,ar,bl,br,rhs)
     INTEGER, INTENT(in):: jdim,kdim,idim,nummem
     REAL, INTENT(inout) :: &
          al(2,jdim,kdim,idim),ar(2,jdim,kdim,idim),&
          bl(2,jdim,kdim,idim),br(2,jdim,kdim,idim)
     REAL,INTENT(in) :: timestep(jdim-1,kdim-1,idim-1)
     REAL, INTENT(inout) :: rhs(jdim-1,kdim-1,idim-1,nummem),d(nummem,nummem, jdim,kdim,idim)

     INTEGER :: j,k,m
     INTEGER, parameter :: i=1
     REAL  :: y(nummem),v(nummem)


     !DO i=1,idim-1
        DO k=1,kdim-1
           DO j=1,jdim-1 
              !rhs(j,k,i,:) = rhs(j,k,i,:)!*timestep(j,k,i)
              d(:,:,j,k,i) = -d(:,:,j,k,i)!*timestep(j,k,i)
              do m=1,nummem
                 d(m,m,j,k,i) = 1./timestep(j,k,i) + d(m,m,j,k,i)
              enddo
#ifdef CRAP
              al(:,j,k,i) = al(:,j,k,i)*timestep(j,k,i)
              ar(:,j,k,i) = ar(:,j,k,i)*timestep(j,k,i)
              bl(:,j,k,i) = bl(:,j,k,i)*timestep(j,k,i)
              br(:,j,k,i) = br(:,j,k,i)*timestep(j,k,i)
#endif
           ENDDO
        ENDDO
     !ENDDO
     !DO i=1,idim-1
        DO k=1,kdim-1
           DO j=1,jdim-1
              CALL lu_ijk(d(1,1,j,k,i),7)
           ENDDO
        ENDDO
     !ENDDO
     
    
     ! forward sweep
     DO k=1,kdim-1
        DO j=1,jdim-1
           y = 0
           IF(j>1) THEN
              y(1:6) = -al(1,j,k,i)*rhs(j-1,k,i,1:6)
              y(7) = -al(2,j,k,i)*rhs(j-1,k,i,7)
           ENDIF
           IF(k>1) THEN
              y(1:6) = y(1:6) - bl(1,j,k,i)*rhs(j,k-1,i,1:6)
              y(7) = y(7)  - bl(2,j,k,i)*rhs(j,k-1,i,7)
           ENDIF
           y(1:7) = rhs(j,k,i,1:7) + y(1:7)
           CALL lu_solv(d(1,1,j,k,i),y,7,v)
           rhs(j,k,i,:) = v(:)
        ENDDO
     ENDDO

     ! backward
     DO k=kdim-1,1,-1
        DO j=jdim-1,1,-1
           y  = 0 
           IF(j<jdim-1) THEN
              y(1:6)  = y(1:6) + ar(1,j,k,i)*rhs(j+1,k,i,1:6)
              y(7) = y(7) + ar(2,j,k,i)*rhs(j+1,k,i,7)
           ENDIF
           IF(k<kdim-1) THEN
              y(1:6)= y(1:6) + br(1,j,k,i)*rhs(j,k+1,i,1:6)
              y(7)  = y(7)   + br(2,j,k,i)*rhs(j,k+1,i,7)
           ENDIF
           CALL lu_solv(d(1,1,j,k,i),y,7,v)
           rhs(j,k,i,1:7) = rhs(j,k,i,1:7) - v(1:7)
        ENDDO
     ENDDO
     RETURN
   END SUBROUTINE sgs_solver_2d

   ! dump movie_block
   SUBROUTINE  kws_mov_block(jdim,kdim,idim,nummem,nblk,q,x,y,z,vist3d,ux,turre,tursav,&
        smin,iunit)
     INTEGER,INTENT(in) :: iunit
     INTEGER,INTENT(in)::jdim,kdim,idim,nummem, nblk
     REAL,INTENT(in) :: x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),&
          q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim),&
          ux(jdim-1,kdim-1,idim-1,9),turre(jdim,kdim,idim,nummem),&
          smin(jdim-1,kdim-1,idim-1),tursav(jdim,kdim,idim,nummem)
     
     
     ! common block
     COMMON /twod/ i2d
     INTEGER ::i2d
     COMMON /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
     INTEGER ::movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie
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
     REAL :: c2b,c2bp,tt
     !REAL ::  fnu(jdim),re,utau(jdim),vel,tau_w
     REAL :: fnuloc,Rel(jdim),cf(jdim)
     !REAL :: deltabl, coef
     COMMON /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),&
          tur10(7)
     REAL :: p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv,&
          tur10

#ifdef CHECK_SOURCE_TERM
     !#ifdef DIST_MPI
     !IF(mblk2nd(nblk)==myid) THEN
        !#endi 
        WRITE(plt_title,'(2(1x,A,es12.5))')'Re=',reue,'Mach=',xmach
        WRITE(iunit,'(A)')'title="'//TRIM(ADJUSTL(plt_title))//'"'
        IF(nummem==7.OR.nummem==0) THEN
           WRITE(iunit,'(A)')'variables="x","y","z","dist","q1", "q2", "q3", "q4", "q5", '//&
                '"dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz", '//&
                '"mu", "mut", "t11", "t22", "t33", "t12", "t23", "t13", "zeta",'//&
                '"11_prd","11_dis", "11_ck", "11_phi", "11_c7", "11_c4", "11_c3", "11_c1",, "11_c2",'//&
                '"22_prd","22_dis", "22_ck", "22_phi", "22_c7", "22_c4", "22_c3", "22_c1",, "22_c2",'//&
                '"33_prd","33_dis", "33_ck", "33_phi", "33_c7", "33_c4", "33_c3", "33_c1",, "33_c2",'//&
                '"12_prd","12_dis", "12_ck", "12_phi", "12_c7", "12_c4", "12_c3", "12_c1",, "12_c2",'//&
                '"23_prd","23_dis", "23_ck", "23_phi", "23_c7", "23_c4", "23_c3", "23_c1",, "23_c2",'//&
                '"13_prd","13_dis", "13_ck", "13_phi", "13_c7", "13_c4", "13_c3", "13_c1" , "13_c2",'//&
                '"src1", "src2", "src3", "src4", "src5", "src6","src7",' // &
                '"rhs1", "rhs2", "rhs3", "rhs4", "rhs5", "rhs6","rhs7"'
                
        ELSEIF(nummem==2) THEN
           WRITE(iunit,'(A)')'variables="x","y","z","dist","q1", "q2", "q3", "q4", "q5", '//&
                '"dudx","dudy","dudz","dvdx","dvdy","dvdz","dwdx","dwdy","dwdz", '//&
                ',"mu","mut", "omg", "tke"'
        ELSE
           STOP"only support nummem=2 or 7,0"
        ENDIF
        IF(i2d==1) THEN
           WRITE(iunit,'(A,i5,A,i5)')"ZONE  I=", jdim-1," J=", kdim-1
        ELSE
           WRITE(iunit,'(A,i5,A,i5,A,i5)')"ZONE  I=", idim-1, " J=",jdim-1," K=", kdim-1
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
                 IF(nummem==0) THEN
                    WRITE(iunit,'(25(1x,es12.5))')xc,yc,zc,smin(j,k,i),q(j,k,i,1:5), (0.0, m=1,9),&
                         fnuloc,0.0, (0.0, m=1,7)
                 ELSE
                    if(allocated(source_items).and.allocated(source)) then
                     
                    WRITE(iunit,'(25(1x,es12.5))')xc,yc,zc,smin(j,k,i),q(j,k,i,1:5), ux(j,k,i,1:9),&
                         fnuloc,vist3d(j,k,i), turre(j,k,i,1:nummem),source_items(1:9,1:6,j,k,i),&
                         source(j,k,i,1:7),turre(j,k,i,1:7)-tursav(j,k,i,1:7)
                    else
                    WRITE(iunit,'(25(1x,es12.5))')xc,yc,zc,smin(j,k,i),q(j,k,i,1:5), ux(j,k,i,1:9),&
                         fnuloc,vist3d(j,k,i), turre(j,k,i,1:nummem),(0.0,m=1,54),&
                         (0.0,m=1,7),turre(j,k,i,1:7)-tursav(j,k,i,1:7)
                    endif
 
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !#ifdef DIST_MPI
     !ENDIF
     !#endif
#endif
   END SUBROUTINE kws_mov_block

   SUBROUTINE kws_dump_movie(nframe,jdim,kdim,idim,nummem,nblk,q,x,y,z,vist3d,ux,turre,tursav,&
     smin)
     INTEGER,INTENT(in)::jdim,kdim,idim,nummem, nblk,nframe
     REAL,INTENT(in) :: x(jdim,kdim,idim),y(jdim,kdim,idim),z(jdim,kdim,idim),&
          q(jdim,kdim,idim,5),vist3d(jdim,kdim,idim),&
          ux(jdim-1,kdim-1,idim-1,9),turre(jdim,kdim,idim,nummem),tursav(jdim,kdim,idim,nummem),&
          smin(jdim-1,kdim-1,idim-1)
     
     COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
     REAL ::  title,rkap,xmach,alpha__,beta__,dt,fmax
     INTEGER :: nit,ntt, &
          idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
          iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv
     
     INTEGER,SAVE :: iunit = 2009
     INTEGER,SAVE :: opened=-1
     CHARACTER(len=30) :: fname,cblk
     
!   ***dump movie file at end of cycle (with -END- in the name)
     IF (nframe==0) then
        fname='END'
     ELSE
        WRITE(fname,'(I4.4)')nframe
     ENDIF
!   ***
     WRITE(cblk, '(I4.4)')nblk
     
     OPEN(iunit,file="mov-"//TRIM(ADJUSTL(fname))//"-"//TRIM(ADJUSTL(cblk))//".plt")
     CALL kws_mov_block(jdim,kdim,idim,nummem,nblk,q,x,y,z,vist3d,ux,turre,tursav,&
          smin,iunit)
     close(iunit)
     
   END SUBROUTINE kws_dump_movie

   ! the new approach to implement the triple correlation terms
   ! using Launder-Reese-Rodi formulation. Equation(6.50) Wilcox Ed 3rd.
   SUBROUTINE kws_cijk(jdim,kdim,idim,nummem,q,qj0,qk0,qi0,turre,tke,&
        sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs,d,al,ar,bl,br)
     
     INTEGER,INTENT(in) :: jdim,idim,kdim,nummem
     REAL,INTENT(in) :: q(jdim,kdim,idim,5),&
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),&
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4),&
         tke(0:jdim,0:kdim,0:idim), &
         sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),&
         vj0(kdim,idim-1,4),vk0(jdim,idim-1,4),&
         vi0(jdim,kdim,4), &
         fmu(0:jdim,0:kdim,0:idim)
    REAL, INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)
    REAL, INTENT(inout) :: d(nummem,nummem,jdim,kdim,idim)
    REAL, INTENT(out), DIMENSION(2,jdim,kdim,idim) :: al,ar,bl,br

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
    COMMON /twod/ i2d
    INTEGER :: i2d

    !local variables
    INTEGER :: i,j,k,m,n,ii,jj,kk,mm,ik,jk,im,jm,km,ij
  
    REAL:: diff(MAX(jdim,kdim,idim),nummem)
    REAL:: tijm(3,6),cijk(3,6),thalf(6)
    REAL:: xmut(MAX(jdim,kdim,idim))
    REAL:: xmu_ave,xmut_ave,vleft,vright,tke_ave,ome_ave,rho_ave
    REAL:: areax2,areay2,areaz2,area2,rvol,xma_re
    
    REAL, parameter :: one_third = 0.3333333333333333
    REAL, PARAMETER :: two_third = 2./3.
    REAL :: rv_l, rv_r,rv_c,xl,yl,zl, xr,yr,zr,xc,yc,zc
    REAL :: xmu,flux(max(jdim,kdim,idim),nummem)
    REAL :: xcoef(max(jdim,kdim,idim),nummem)
    INTEGER,parameter :: m2i(6)=(/1,2,3,1,2,1/),m2j(6)=(/1,2,3,2,3,3/)
    INTEGER,parameter :: ij2m(3,3)=RESHAPE((/1,4,6,4,2,5,6,5,3/),(/3,3/))

    xma_re= xmach/reue
    ! diffusion terms in the j-direction
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             IF(j==1) THEN
                rho_ave = 0.5*(qj0(k,i,1,1)+q(j,k,i,1))
             ELSEif(j==jdim) then
                rho_ave  = 0.5*(qj0(k,i,1,3)+q(j-1,k,i,1))
             ELSE
                rho_ave  = 0.5*(q(j-1,k,i,1)+q(j,k,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j-1,k,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j-1,k,i,7))
             xmut(j) = MAX(0.0,rho_ave*tke_ave/ome_ave)
             diff(j,:) =turre(j,k,i,:) -  turre(j-1,k,i,:)
             IF(j==1) THEN
                vleft = vj0(k,i,1)
             ELSE
                vleft = vol(j-1,k,i)
             ENDIF
             IF(j==jdim) THEN
                vright = vj0(k,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             
             xmu = xmu_ave             
             flux(j,1:6) = diff(j,1:6)*sj(j,k,i,4)**2*rvol*xmu  !laminar part diffusion
             xcoef(j,1) =  sj(j,k,i,4)**2*rvol*(xmu+xmut(j)*sigma_star)
             DO m=1,6
                thalf(m) = 0.5*(turre(j-1,k,i,m)+turre(j,k,i,m))
                DO n=1,3  ! derivatives of tij
                   tijm(n,m) = diff(j,m)*sj(j,k,i,n)*sj(j,k,i,4)*rvol
                ENDDO
             ENDDO
             DO m=1,6
                ii=m2i(m)
                jj=m2j(m)
                ij = ij2m(ii,jj)
                IF(ij/=m) STOP "ij/=m"
                DO kk=1,3
                   cijk(kk,m) = 0.
                   jk = ij2m(jj,kk)
                   ik = ij2m(ii,kk)
                   DO mm=1,3
                      im = ij2m(ii,mm)
                      jm = ij2m(jj,mm)
                      km = ij2m(kk,mm)
                      cijk(kk,m) =  cijk(kk,m)+&
                           thalf(im)*tijm(mm,jk)+&
                           thalf(jm)*tijm(mm,ik)+&
                           thalf(km)*tijm(mm,ij)
                   ENDDO
                ENDDO
             ENDDO
             DO m=1,6
                flux(j,m) = flux(j,m)-sigma_star*xmut(j)/(tke_ave+1e-20)*&
                     (cijk(1,m)*sj(j,k,i,1)+cijk(2,m)*sj(j,k,i,2)+cijk(3,m)*sj(j,k,i,3))*sj(j,k,i,4)
             ENDDO
             
             xmu = xmu_ave+xmut(j)*sigma
             flux(j,7) = diff(j,7)*sj(j,k,i,4)**2*rvol*xmu
             xcoef(j,2) = sj(j,k,i,4)**2*rvol*xmu 
          ENDDO
          
          !CALL cap_xmut(xmut,0,jdim) 
          DO j=1,jdim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(j+1,:)-flux(j,:))/(q(j,k,i,1)*vol(j,k,i))
             
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(j+1,2) + xcoef(j,2))/vol(j,k,i) 
             al(2,j,k,i) = -xma_re*xcoef(j,2)/vol(j,k,i)
             ar(2,j,k,i) = -xma_re*xcoef(j+1,2)/vol(j,k,i)
             
             al(1,j,k,i) = -xma_re*xcoef(j,1)/vol(j,k,i)
             ar(1,j,k,i) = -xma_re*xcoef(j+1,1)/vol(j,k,i)
             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(j,1)+xcoef(j+1,1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! diffusion terms in the k-direction
    DO i=1,idim-1
       DO j=1,jdim-1
          DO k=1,kdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             IF(k==1) THEN
                rho_ave = 0.5*(qk0(j,i,1,1)+q(j,k,i,1))
             ELSEIF(k==kdim) THEN
                rho_ave  = 0.5*(qk0(j,i,1,3)+q(j,k-1,i,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k-1,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k-1,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k-1,i,7))
             xmut(k) = max(0.0,rho_ave*tke_ave/ome_ave)
             diff(k,:) =turre(j,k,i,:) -  turre(j,k-1,i,:)
             IF(k==1) THEN
                vleft = vk0(j,i,1)
             ELSE
                vleft = vol(j,k-1,i)
             ENDIF
             IF(k==kdim) THEN
                vright = vk0(j,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             !-----------------------------------------
             xmu = xmu_ave             
             flux(k,1:6) = diff(k,1:6)*sk(j,k,i,4)**2*rvol*xmu  !laminar part diffusion
             xcoef(k,1) =  sk(j,k,i,4)**2*rvol*(xmu+xmut(k)*sigma_star)
             DO m=1,6
                thalf(m) = 0.5*(turre(j,k-1,i,m)+turre(j,k,i,m))
                DO n=1,3  ! derivatives of tij
                   tijm(n,m) = diff(k,m)*sk(j,k,i,n)*sk(j,k,i,4)*rvol
                ENDDO
             ENDDO
             DO m=1,6
                ii=m2i(m)
                jj=m2j(m)
                ij = ij2m(ii,jj)
                IF(ij/=m) STOP "ij/=m"
                DO kk=1,3
                   cijk(kk,m) = 0.
                   jk = ij2m(jj,kk)
                   ik = ij2m(ii,kk)
                   DO mm=1,3
                      im = ij2m(ii,mm)
                      jm = ij2m(jj,mm)
                      km = ij2m(kk,mm)
                      cijk(kk,m) =  cijk(kk,m)+&
                           thalf(im)*tijm(mm,jk)+&
                           thalf(jm)*tijm(mm,ik)+&
                           thalf(km)*tijm(mm,ij)
                   ENDDO
                ENDDO
             ENDDO
             DO m=1,6
                flux(k,m) = flux(k,m)-sigma_star*xmut(k)/(tke_ave+1e-20)*&
                     (cijk(1,m)*sk(j,k,i,1)+cijk(2,m)*sk(j,k,i,2)+cijk(3,m)*sk(j,k,i,3))*sk(j,k,i,4)
             ENDDO

             !-----------------------------------------

             xmu = xmu_ave+xmut(k)*sigma
             flux(k,7) = diff(k,7)*sk(j,k,i,4)**2*rvol*xmu
             xcoef(k,2) = sk(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO k=1,kdim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(k+1,:)-flux(k,:))/(q(j,k,i,1)*vol(j,k,i))
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(k+1,2) + xcoef(k,2))/vol(j,k,i) 
             bl(2,j,k,i) = -xma_re*xcoef(k,2)/vol(j,k,i)
             br(2,j,k,i) = -xma_re*xcoef(k+1,2)/vol(j,k,i)
             
             bl(1,j,k,i) = -xma_re*xcoef(k,1)/vol(j,k,i)
             br(1,j,k,i) = -xma_re*xcoef(k+1,1)/vol(j,k,i)
             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(k,1)+xcoef(k+1,1))/vol(j,k,i)
             ENDDO
                          
          ENDDO
       ENDDO
    ENDDO

    IF(i2d/=1) THEN
    ! diffusion terms in the i-direction
    DO k=1,kdim-1
       DO j=1,jdim-1
          DO i=1,idim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             IF(i==1) THEN
                rho_ave = 0.5*(qi0(j,k,1,1)+q(j,k,i,1))
             ELSEIF(i==idim) THEN
                rho_ave  = 0.5*(qi0(j,k,1,3)+q(j,k,i-1,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k,i-1,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k,i-1))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k,i-1,7))
             xmut(i) = max(0.,rho_ave*tke_ave/ome_ave)

             diff(i,:) =turre(j,k,i,:) -  turre(j,k,i-1,:)
             IF(i==1) THEN
                vleft = vi0(j,k,1)
             ELSE
                vleft = vol(j,k,i-1)
             ENDIF
             IF(i==idim) THEN
                vright = vi0(j,k,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             xmu = xmu_ave+xmut(i)*sigma_star
             flux(i,1:6) = diff(i,1:6)*si(j,k,i,4)**2*rvol*xmu
             xcoef(i,1) = si(j,k,i,4)**2*rvol*xmu 
             xmu = xmu_ave+xmut(i)*sigma
             flux(i,7) = diff(i,7)*si(j,k,i,4)**2*rvol*xmu
             xcoef(i,2) = si(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO i=1,idim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(i+1,:)-flux(i,:))/(q(j,k,i,1)*vol(j,k,i))
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(i+1,2) + xcoef(i,2))/vol(j,k,i) 

             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(i,1)+xcoef(i+1,1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO    
    ENDIF
    
    RETURN
   END SUBROUTINE kws_cijk


   ! the new approach to implement the triple correlation terms
   ! using Launder-Reese-Rodi formulation 1. Eq(6.49) Wilcox Ed 3rd.
   SUBROUTINE kws_cijk0(jdim,kdim,idim,nummem,q,qj0,qk0,qi0,turre,tke,&
        sj,sk,si,vol,vj0,vk0,vi0,fmu,rhs,d,al,ar,bl,br)
     
     INTEGER,INTENT(in) :: jdim,idim,kdim,nummem
     REAL,INTENT(in) :: q(jdim,kdim,idim,5),&
         turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem),&
         qj0(kdim,idim-1,5,4),qk0(jdim,idim-1,5,4),qi0(jdim,kdim,5,4),&
         tke(0:jdim,0:kdim,0:idim), &
         sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
         si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1),&
         vj0(kdim,idim-1,4),vk0(jdim,idim-1,4),&
         vi0(jdim,kdim,4), &
         fmu(0:jdim,0:kdim,0:idim)
    REAL, INTENT(out) :: rhs(jdim-1,kdim-1,idim-1,nummem)
    REAL, INTENT(inout) :: d(nummem,nummem,jdim,kdim,idim)
    REAL, INTENT(out), DIMENSION(2,jdim,kdim,idim) :: al,ar,bl,br

    !common block needed
    COMMON /reyue/ reue,tinf
    REAL :: reue,tinf
    COMMON /info/ title(20),rkap(3),xmach
    REAL::title,rkap,xmach
    COMMON /twod/ i2d
    INTEGER :: i2d

    !local variables
    INTEGER :: i,j,k,m,n,ii,jj,kk,mm,ik,jk,im,jm,km,ij
  
    REAL:: diff(MAX(jdim,kdim,idim),nummem)
    REAL:: tijm(3,6),cijk(3,6),thalf(6)
    REAL:: xmut(MAX(jdim,kdim,idim))
    REAL:: xmu_ave,xmut_ave,vleft,vright,tke_ave,ome_ave,rho_ave
    REAL:: areax2,areay2,areaz2,area2,rvol,xma_re
    
    REAL, parameter :: one_third = 0.3333333333333333
    REAL, PARAMETER :: two_third = 2./3.
    REAL :: rv_l, rv_r,rv_c,xl,yl,zl, xr,yr,zr,xc,yc,zc
    REAL :: xmu,flux(max(jdim,kdim,idim),nummem)
    REAL :: xcoef(max(jdim,kdim,idim),nummem)
    INTEGER,parameter :: m2i(6)=(/1,2,3,1,2,1/),m2j(6)=(/1,2,3,2,3,3/)
    INTEGER,parameter :: ij2m(3,3)=RESHAPE((/1,4,6,4,2,5,6,5,3/),(/3,3/))

    xma_re= xmach/reue
    ! diffusion terms in the j-direction
    DO i=1,idim-1
       DO k=1,kdim-1
          DO j=1,jdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j-1,k,i))
             IF(j==1) THEN
                rho_ave = 0.5*(qj0(k,i,1,1)+q(j,k,i,1))
             ELSEif(j==jdim) then
                rho_ave  = 0.5*(qj0(k,i,1,3)+q(j-1,k,i,1))
             ELSE
                rho_ave  = 0.5*(q(j-1,k,i,1)+q(j,k,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j-1,k,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j-1,k,i,7))
             xmut(j) = MAX(0.0,rho_ave*tke_ave/ome_ave)
             diff(j,:) =turre(j,k,i,:) -  turre(j-1,k,i,:)
             IF(j==1) THEN
                vleft = vj0(k,i,1)
             ELSE
                vleft = vol(j-1,k,i)
             ENDIF
             IF(j==jdim) THEN
                vright = vj0(k,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             
             xmu = xmu_ave
             flux(j,1:6) = diff(j,1:6)*sj(j,k,i,4)**2*rvol*xmu  !laminar part diffusion
             xcoef(j,1) =  sj(j,k,i,4)**2*rvol*(xmu+xmut(j)*sigma_star)
             DO m=1,6
                thalf(m) = 0.5*(turre(j-1,k,i,m)+turre(j,k,i,m))
                DO n=1,3  ! derivatives of tij
                   tijm(n,m) = diff(j,m)*sj(j,k,i,n)*sj(j,k,i,4)*rvol
                ENDDO
             ENDDO
             cijk = 0
             DO m=1,6
                ii=m2i(m)
                jj=m2j(m)
                ij = ij2m(ii,jj)
                IF(ij/=m) STOP "ij/=m"
                DO kk=1,3
                   cijk(kk,m) = 0.
                   jk = ij2m(jj,kk)
                   ik = ij2m(ii,kk)
                   cijk(kk,m) = tijm(ii,jk)+tijm(jj,ik)+tijm(kk,ij)
                ENDDO
             ENDDO
             
             DO m=1,6
                flux(j,m) = flux(j,m)+sigma_star*xmut(j)*&
                     (cijk(1,m)*sj(j,k,i,1)+cijk(2,m)*sj(j,k,i,2)+cijk(3,m)*sj(j,k,i,3))*sj(j,k,i,4)
             ENDDO
             
             xmu = xmu_ave+xmut(j)*sigma
             flux(j,7) = diff(j,7)*sj(j,k,i,4)**2*rvol*xmu
             xcoef(j,2) = sj(j,k,i,4)**2*rvol*xmu 
          ENDDO
          
          !CALL cap_xmut(xmut,0,jdim) 
          DO j=1,jdim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(j+1,:)-flux(j,:))/(q(j,k,i,1)*vol(j,k,i))
             
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(j+1,2) + xcoef(j,2))/vol(j,k,i) 
             al(2,j,k,i) = -xma_re*xcoef(j,2)/vol(j,k,i)
             ar(2,j,k,i) = -xma_re*xcoef(j+1,2)/vol(j,k,i)
             
             al(1,j,k,i) = -xma_re*xcoef(j,1)/vol(j,k,i)
             ar(1,j,k,i) = -xma_re*xcoef(j+1,1)/vol(j,k,i)
             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(j,1)+xcoef(j+1,1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    
    ! diffusion terms in the k-direction
    DO i=1,idim-1
       DO j=1,jdim-1
          DO k=1,kdim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k-1,i))
             IF(k==1) THEN
                rho_ave = 0.5*(qk0(j,i,1,1)+q(j,k,i,1))
             ELSEIF(k==kdim) THEN
                rho_ave  = 0.5*(qk0(j,i,1,3)+q(j,k-1,i,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k-1,i,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k-1,i))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k-1,i,7))
             xmut(k) = max(0.0,rho_ave*tke_ave/ome_ave)
             diff(k,:) =turre(j,k,i,:) -  turre(j,k-1,i,:)
             IF(k==1) THEN
                vleft = vk0(j,i,1)
             ELSE
                vleft = vol(j,k-1,i)
             ENDIF
             IF(k==kdim) THEN
                vright = vk0(j,i,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             !-----------------------------------------
             xmu = xmu_ave
             flux(k,1:6) = diff(k,1:6)*sk(j,k,i,4)**2*rvol*xmu  !laminar part diffusion
             xcoef(k,1) =  sk(j,k,i,4)**2*rvol*(xmu+xmut(k)*sigma_star)
             DO m=1,6
                thalf(m) = 0.5*(turre(j,k-1,i,m)+turre(j,k,i,m))
                DO n=1,3  ! derivatives of tij
                   tijm(n,m) = diff(k,m)*sk(j,k,i,n)*sk(j,k,i,4)*rvol
                ENDDO
             ENDDO
             DO m=1,6
                ii=m2i(m)
                jj=m2j(m)
                ij = ij2m(ii,jj)
                IF(ij/=m) STOP "ij/=m"
                DO kk=1,3
                   cijk(kk,m) = 0.
                   jk = ij2m(jj,kk)
                   ik = ij2m(ii,kk)
                   cijk(kk,m) = tijm(ii,jk)+tijm(jj,ik)+tijm(kk,ij)
                ENDDO

             ENDDO
             DO m=1,6
                flux(k,m) = flux(k,m)+sigma_star*xmut(k)*&
                     (cijk(1,m)*sk(j,k,i,1)+cijk(2,m)*sk(j,k,i,2)+cijk(3,m)*sk(j,k,i,3))*sk(j,k,i,4)
             ENDDO

             !-----------------------------------------

             xmu = xmu_ave+xmut(k)*sigma
             flux(k,7) = diff(k,7)*sk(j,k,i,4)**2*rvol*xmu
             xcoef(k,2) = sk(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO k=1,kdim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(k+1,:)-flux(k,:))/(q(j,k,i,1)*vol(j,k,i))
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(k+1,2) + xcoef(k,2))/vol(j,k,i) 
             bl(2,j,k,i) = -xma_re*xcoef(k,2)/vol(j,k,i)
             br(2,j,k,i) = -xma_re*xcoef(k+1,2)/vol(j,k,i)
             
             bl(1,j,k,i) = -xma_re*xcoef(k,1)/vol(j,k,i)
             br(1,j,k,i) = -xma_re*xcoef(k+1,1)/vol(j,k,i)
             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(k,1)+xcoef(k+1,1))/vol(j,k,i)
             ENDDO
                          
          ENDDO
       ENDDO
    ENDDO

    IF(i2d/=1) THEN
    ! diffusion terms in the i-direction
    DO k=1,kdim-1
       DO j=1,jdim-1
          DO i=1,idim
             xmu_ave = 0.5*(fmu(j,k,i)+fmu(j,k,i-1))
             IF(i==1) THEN
                rho_ave = 0.5*(qi0(j,k,1,1)+q(j,k,i,1))
             ELSEIF(i==idim) THEN
                rho_ave  = 0.5*(qi0(j,k,1,3)+q(j,k,i-1,1))
             ELSE
                rho_ave  = 0.5*(q(j,k,i,1)+q(j,k,i-1,1))
             ENDIF
             tke_ave = 0.5*(tke(j,k,i)+tke(j,k,i-1))
             ome_ave = 0.5*(turre(j,k,i,7)+turre(j,k,i-1,7))
             xmut(i) = max(0.,rho_ave*tke_ave/ome_ave)

             diff(i,:) =turre(j,k,i,:) -  turre(j,k,i-1,:)
             IF(i==1) THEN
                vleft = vi0(j,k,1)
             ELSE
                vleft = vol(j,k,i-1)
             ENDIF
             IF(i==idim) THEN
                vright = vi0(j,k,3)
             ELSE
                vright = vol(j,k,i)
             ENDIF
             rvol = 1./(0.5*(vleft + vright))
             xmu = xmu_ave+xmut(i)*sigma_star
             flux(i,1:6) = diff(i,1:6)*si(j,k,i,4)**2*rvol*xmu
             xcoef(i,1) = si(j,k,i,4)**2*rvol*xmu 
             xmu = xmu_ave+xmut(i)*sigma
             flux(i,7) = diff(i,7)*si(j,k,i,4)**2*rvol*xmu
             xcoef(i,2) = si(j,k,i,4)**2*rvol*xmu 
          ENDDO

          !call cap_xmut(xmut,0,kdim) 
          
          DO i=1,idim-1
             rhs(j,k,i,:) = rhs(j,k,i,:) + xma_re*(flux(i+1,:)-flux(i,:))/(q(j,k,i,1)*vol(j,k,i))
             d(7,7,j,k,i) = d(7,7,j,k,i) - xma_re*(xcoef(i+1,2) + xcoef(i,2))/vol(j,k,i) 

             DO m=1,6
                d(m,m,j,k,i)= d(m,m,j,k,i) - xma_re*(xcoef(i,1)+xcoef(i+1,1))/vol(j,k,i)
             ENDDO
          ENDDO
       ENDDO
    ENDDO    
    ENDIF
    
    RETURN
   END SUBROUTINE kws_cijk0

!   The following subroutine - NO LONGER USED!!!!
   SUBROUTINE set_wallbc(jdim,kdim,idim,nummem,q,turb,smin,vist3d,tj0,tk0,ti0)
     INTEGER,INTENT(in) :: jdim,kdim,idim, nummem
     REAL,INTENT(in) :: turb(jdim,kdim,idim,nummem),smin(jdim-1,kdim-1,idim-1),&
          q(jdim,kdim,idim,5), vist3d(jdim,kdim,idim)
     REAL, INTENT(inout) :: tj0(kdim,idim-1,nummem,4),&
          tk0(jdim,idim-1,nummem,4),ti0(jdim,kdim,nummem,4)
     COMMON /twod/ i2d
     INTEGER :: i2d
     COMMON /fluid/ gamma,gm1,gp1,gm1g,gp1g,ggm1
     REAL::gamma,gm1,gp1,gm1g,gp1g,ggm1
     COMMON /fluid2/ pr,prt,cbar
     REAL ::pr,prt,cbar
     COMMON /reyue/reue,tinf,ivisc(3)
     REAL :: reue,tinf
     INTEGER::ivisc
     COMMON /info/ title(20),rkap(3),xmach
     REAL :: title,rkap,xmach

     INTEGER :: i,j,k,iv
     REAL :: tt,fnu,tke1,tke2,dkdy,c2b,c2bp,RE
     REAL :: dist

     c2b=cbar/tinf
     c2bp=c2b+1.0
     RE = REUE/XMACH
     ! apply the wall boundary condition for subiterations

     ! check jmin boundary
     j=1
     DO i=1,idim-1
        DO k=1,kdim-1
           IF(tj0(k,i,1,2)==HUGE(tj0(k,i,1,1))) THEN  ! "huge" is a wall bc indicator, c.f. bc2004.F
              DO iv=1,6
                 tj0(k,i,iv,1) = -turb(j,k,i,iv)
                 tj0(k,i,iv,2) = HUGE(tke1)
              ENDDO
              tt=gamma*q(j,k,i,5)/q(j,k,i,1)
              fnu=c2bp*tt*SQRT(tt)/(c2b+tt)
              dist=abs(smin(j,k,i))
              tj0(k,i,7,1) = &
                   2.*(60.*fnu/(re**2*q(j,k,i,1)*beta_0*dist**2))- &
                   turb(j,k,i,7)
              tj0(k,i,7,2) = 2*tj0(k,i,7,1) - turb(j,k,i,7) 
           ENDIF
        ENDDO
     ENDDO

     ! check jmin boundary
     j=jdim-1
     DO i=1,idim-1
        DO k=1,kdim-1
           IF(tj0(k,i,1,4)==HUGE(tj0(k,i,1,1))) THEN  ! "huge" is a wall bc indicator, c.f. bc2004.F
              DO iv=1,6
                 tj0(k,i,iv,3) = -turb(j,k,i,iv)
                 tj0(k,i,iv,4) = HUGE(tke1)
              ENDDO
              tt=gamma*q(k,j,i,5)/q(j,k,i,1)
              fnu=c2bp*tt*SQRT(tt)/(c2b+tt)
              dist=abs(smin(j,k,i))
              tj0(k,i,7,3) = &
                   2.*(60.*fnu/(re**2*q(j,k,i,1)*beta_0*dist**2))- &
                   turb(j,k,i,7)
              tj0(k,i,7,4) = 2*tj0(k,i,7,3) - turb(j,k,i,7)
           ENDIF
        ENDDO
     ENDDO

     !check for kmin boundary
     k=1
     DO i=1,idim-1
        DO j=1,jdim-1
           IF(tk0(j,i,1,2)==HUGE(tk0(j,i,1,1))) THEN
              DO iv=1,6
                 tk0(j,i,iv,1) = -turb(j,k,i,iv)
                 tk0(j,i,iv,2) = HUGE(tke1)
              ENDDO
              tt=gamma*q(j,k,i,5)/q(j,k,i,1)
              fnu=c2bp*tt*SQRT(tt)/(c2b+tt)
              dist=abs(smin(j,k,i))
              tk0(j,i,7,1) = &
                   2.*(60.*fnu/(re**2*q(j,k,i,1)*beta_0*dist**2))- &
                   turb(j,k,i,7)
              tk0(j,i,7,2) = 2*tk0(j,i,7,1) - turb(j,k,i,7)
           ENDIF
        ENDDO
     ENDDO

     !check for kmin boundary
     k=kdim-1
     DO i=1,idim-1
        DO j=1,jdim-1
           IF(tk0(j,i,1,4)==HUGE(tk0(j,i,1,1))) THEN
              DO iv=1,6
                 tk0(j,i,iv,3) = -turb(j,k,i,iv)
                 tk0(j,i,iv,4) = HUGE(tke1)
              ENDDO
              tt=gamma*q(j,k,i,5)/q(j,k,i,1)
              fnu=c2bp*tt*SQRT(tt)/(c2b+tt)
              dist=abs(smin(j,k,i))
              tk0(j,i,7,3) = &
                   2.*(60.*fnu/(re**2*q(j,k,i,1)*beta_0*dist**2))- &
                   turb(j,k,i,7)
              tk0(j,i,7,4) = 2*tk0(j,i,7,3)- turb(j,k,i,7)
              
           ENDIF
        ENDDO
     ENDDO
     
     IF(i2d==0) THEN
        STOP "to be implemented, set_wallbc"
     ENDIF
   END SUBROUTINE set_wallbc
 
   subroutine redefine_shear(jdim,kdim,idim,nummem,ux,vist3d,turb)
      integer, intent(in) :: jdim,kdim,idim,nummem
      real, intent(in) :: ux(jdim-1,kdim-1,idim-1,9),vist3d(jdim,kdim,idim)
      real, intent(inout) :: turb(jdim,kdim,idim, nummem)


    COMMON /reyue/ reue,tinf
    REAL::reue,tinf
    COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
         idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
         iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
    REAL:: title,rkap, xmach, alpha__,beta__, dt, fmax
    INTEGER ::nit,ntt, &
         idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
         iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv
    

    integer :: i,j,k
    real :: re
    return
    re = xmach/reue
    do i=1, idim
       do k=1,kdim
           do j=1,jdim
              turb(j,k,i,4) = re* vist3d(j,k,i)*(ux(j,k,i,2)+ux(j,k,i,4))
              turb(j,k,i,5) = re* vist3d(j,k,i)*(ux(j,k,i,6)+ux(j,k,i,8))
              turb(j,k,i,6) = re* vist3d(j,k,i)*(ux(j,k,i,3)+ux(j,k,i,7))
           enddo
       enddo
    enddo

   end subroutine 

!   The following subroutine - NO LONGER USED!!!!
   SUBROUTINE kws_init_from_2eq(lfrom2eq,nbl,jdim,kdim,idim,&
        nummem,vist3d,tursav)
     LOGICAL, intent(in) :: lfrom2eq
     INTEGER,INTENT(in)  :: nbl,jdim, kdim, idim,nummem
     REAL :: vist3d(jdim,kdim,idim)
     REAL :: tursav(jdim,kdim,idim,nummem)

     COMMON /mydist2/ nnodes,myhost,myid,mycomm
     INTEGER::nnodes,myhost,myid,mycomm

     character(len=30) :: str
     INTEGER :: j,k,i,m
     lkzstm_from_2eq = lfrom2eq
     IF(.NOT.lfrom2eq) RETURN
     IF(myid == myhost) THEN
        WRITE(11,*)"load 2eq solution ..."
     ENDIF
     WRITE(str,*)nbl
     OPEN(2001, file='2eq-'//TRIM(ADJUSTL(str))//".dat",form='unformatted',status='old')
     READ(2001) (((vist3d(j,k,i),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
     READ(2001) (((tursav(j,k,i,7),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
     READ(2001) (((tursav(j,k,i,1),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
     close(2001)
     tursav(:,:,:,1) = -tursav(:,:,:,1)*2./3.
   END SUBROUTINE kws_init_from_2eq

   
!   The following subroutine - NO LONGER USED!!!!
   SUBROUTINE kws_save_2eq(nbl,jdim,kdim,idim,nummem,vist3d, tursav)
     INTEGER,INTENT(in)  :: nbl,jdim, kdim, idim,nummem
     REAL :: vist3d(jdim,kdim,idim)
     REAL :: tursav(jdim,kdim,idim,nummem)
     
     INTEGER :: j,k,i
     character(len=30) :: str
     COMMON /mydist2/ nnodes,myhost,myid,mycomm
     INTEGER::nnodes,myhost,myid,mycomm

     
     IF(nummem/=2) return
     IF(myid ==myhost) THEN
        WRITE(11,*)"save 2-eq solution", nbl
     ENDIF
     WRITE(str,*)nbl
     OPEN(2001, file='2eq-'//TRIM(ADJUSTL(str))//".dat",form='unformatted',status='unknown')
     WRITE(2001) (((vist3d(j,k,i),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
     WRITE(2001) (((tursav(j,k,i,1),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
     WRITE(2001) (((tursav(j,k,i,2),j=1,jdim-1),k=1,kdim-1),i=1,idim-1)
     CLOSE(2001)
   END SUBROUTINE kws_save_2eq

   SUBROUTINE kws_init_stress_from_2eq(nbl,jdim,kdim,idim, nummem,ux,vist3d,rho,tursav)
     INTEGER,INTENT(in)  :: jdim, kdim, idim,nummem,nbl
     REAL,INTENT(in) :: vist3d(jdim,kdim,idim),ux(jdim-1,kdim-1,idim-1,9),rho(jdim,kdim,idim)
     REAL,intent(inout) :: tursav(jdim,kdim,idim,nummem)
     
     INTEGER :: j,k,i,m
     
     REAL :: sij(6),xnu,sii,tke23
     COMMON /reyue/ reue,tinf
     REAL::reue,tinf
     COMMON /info/ title(20),rkap(3),xmach,alpha__,beta__,dt,fmax,nit,ntt, &
          idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag, &
          iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
     REAL:: title,rkap, xmach, alpha__,beta__, dt, fmax
     INTEGER ::nit,ntt, &
          idiag,nitfo,iflagts,iflim,nres,levelb,mgflag, &
          iconsf,mseq,ncyc1,levelt,nitfo1,ngam,nsm,iipv

     real :: re

     re = xmach/reue
     call  kws_init_from_2eq(.true.,nbl,jdim,kdim,idim,&
        nummem,vist3d,tursav)

     DO i=1,idim-1
        DO k=1,kdim-1
           DO j=1,jdim-1
              sij(1) = ux(j,k,i,1)
              sij(2) = ux(j,k,i,5)
              sij(3) = ux(j,k,i,9)
              sij(4) = 0.5*(ux(j,k,i,2)+ux(j,k,i,4))
              sij(5) = 0.5*(ux(j,k,i,6)+ux(j,k,i,8))
              sij(6) = 0.5*(ux(j,k,i,3)+ux(j,k,i,7))
              sii = (sij(1) + sij(2) + sij(3))/3.
              tke23 = -tursav(j,k,i,1)
              xnu = vist3d(j,k,i)/rho(j,k,i)*re
              tursav(j,k,i,1) = 2.*xnu*(sij(1) - sii) - tke23
              tursav(j,k,i,2) = 2.*xnu*(sij(2) - sii) - tke23
              tursav(j,k,i,3) = 2.*xnu*(sij(3) - sii) - tke23
              tursav(j,k,i,4) = 2.*xnu*(sij(4))
              tursav(j,k,i,5) = 2.*xnu*(sij(5))
              tursav(j,k,i,6) = 2.*xnu*(sij(6))
           ENDDO
        ENDDO
     ENDDO
     RETURN
   END SUBROUTINE kws_init_stress_from_2eq

   SUBROUTINE kws_dbij_dx(jdim,kdim,idim,nummem, turre,tke,sj,sk,si,vol,dbijdx)
     INTEGER, INTENT(in) :: jdim,kdim,idim,nummem
     REAL, INTENT(in) :: turre(-1:jdim+1,-1:kdim+1,-1:idim+1,nummem)
     REAL, intent(in) :: tke(0:jdim,0:kdim,0:idim)
     REAL, INTENT(in) :: sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
     REAL,INTENT(out) :: dbijdx(jdim-1,kdim-1,idim-1,6,3)
     
     
     INTEGER :: j,k,i,m,n
     REAL :: dbdj(6),xj(3),dbdk(6),xk(3)

     !j-direction
     DO i=1,idim-1
        DO k=1,kdim-1
           DO j=1,jdim-1
              DO m=1,6
                 dbdj(m) = 0.5*(turre(j+1,k,i,m) - turre(j-1,k,i,m))
              ENDDO
              DO m=1,3
                 dbdj(m) = 1./3.*(tke(j+1,k,i) - tke(j-1,k,i)) + dbdj(m)
              ENDDO
              DO m=1,3
                 xj(m) = 0.5*(sj(j,k,i,m)*sj(j,k,i,4) + sj(j+1,k,i,m)*sj(j+1,k,i,4))/vol(j,k,i)
              ENDDO
              DO m = 1,6
                 DO n=1,3
                    dbijdx(j,k,i,m,n) = dbdj(m)*xj(n)
                 ENDDO
              ENDDO
              
           ENDDO
        ENDDO
     ENDDO
     
     !k-direction
     DO i=1,idim-1
        DO k=1,kdim-1
           DO j=1,jdim-1
              DO m=1,6
                 dbdk(m) = 0.5*(turre(j,k+1,i,m) - turre(j,k-1,i,m))
              ENDDO
              DO m=1,3
                 dbdk(m) = 1./3.*(tke(j,k+1,i) - tke(j,k-1,i)) + dbdj(m)
              ENDDO
              DO m=1,3
                 xk(m) = 0.5*(sk(j,k,i,m)*sk(j,k,i,4) + sk(j,k+1,i,m)*sk(j,k+1,i,4))/vol(j,k,i)
              ENDDO
              DO m = 1,6
                 DO n=1,3
                    dbijdx(j,k,i,m,n) = dbdk(m)*xk(n)+dbijdx(j,k,i,m,n)
                 ENDDO
              ENDDO         
           ENDDO
        ENDDO
     ENDDO
     RETURN
   END SUBROUTINE kws_dbij_dx

   SUBROUTINE kws_dbij_dxx(jdim,kdim,idim,nummem,dbijdx,sj,sk,si,vol,dbijdxx)
     INTEGER,INTENT(in)  :: jdim, kdim, idim,nummem
     REAL, INTENT(in) :: dbijdx(jdim-1,kdim-1,idim-1,6,3)
     REAL, INTENT(in) :: sj(jdim,kdim,idim-1,5),sk(jdim,kdim,idim-1,5),&
          si(jdim,kdim,idim,5),vol(jdim,kdim,idim-1)
     REAL,INTENT(out) :: dbijdxx(jdim-1,kdim-1,idim-1,6,3)

     
     INTEGER :: i,j,k,m,n,jp,jm,kp,km
     REAL ::dbdj(6,3),dbdk(6,3),xj(3),xk(3)
     !j-direction
     DO i=1,idim-1
        DO k=1,kdim-1
           DO j=1,jdim-1
              jp=1;jm=1
              IF(j==1) jm=0
              IF(j==jdim-1) jp=0
              DO n=1,3
                 DO m=1,6
                    dbdj(m,n)= (dbijdx(j+jp,k,i,m,n) - dbijdx(j-jm,k,i,m,n))/(jp+jm)
                 ENDDO
              ENDDO
              DO m=1,3
                 xj(m) = 0.5*(sj(j,k,i,m)*sj(j,k,i,4) + sj(j+1,k,i,m)*sj(j+1,k,i,4))/vol(j,k,i)
              ENDDO
              DO m = 1,6
                 DO n=1,3
                    dbijdxx(j,k,i,m,n) = dbdj(m,n)*xj(n)
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     
     !k-direction
     DO i=1,idim-1
        DO k=1,kdim-1
           kp=1;km=1
           IF(k==1) km=0
           IF(k==kdim-1) kp=0
           DO j=1,jdim-1
              DO n=1,3
                 DO m=1,6
                    dbdk(m,n)= (dbijdx(j,k+kp,i,m,n) - dbijdx(j,k-km,i,m,n))/(kp+km)
                 ENDDO
              ENDDO
              DO m=1,3
                 xk(m) = 0.5*(sk(j,k,i,m)*sk(j,k,i,4) + sk(j,k+1,i,m)*sk(j,k+1,i,4))/vol(j,k,i)
              ENDDO
              DO m = 1,6
                 DO n=1,3
                    dbijdxx(j,k,i,m,n) = dbdk(m,n)*xk(n)
                 ENDDO
              ENDDO
              dbijdxx(j,k,i,:,:) = 0.5*dbijdxx(j,k,i,:,:)  ! due to k * bij  is defined as - k*(Cij + 2./3 tke delta ij)/2k
           ENDDO
        ENDDO
     ENDDO

     
   END SUBROUTINE kws_dbij_dxx
END MODULE module_kwstm
