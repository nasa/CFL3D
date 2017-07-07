c  ---------------------------------------------------------------------------
c  CFL3D is a structured-grid, cell-centered, upwind-biased, Reynolds-averaged
c  Navier-Stokes (RANS) code. It can be run in parallel on multiple grid zones
c  with point-matched, patched, overset, or embedded connectivities. Both
c  multigrid and mesh sequencing are available in time-accurate or
c  steady-state modes.
c
c  Copyright 2001 United States Government as represented by the Administrator
c  of the National Aeronautics and Space Administration. All Rights Reserved.
c 
c  The CFL3D platform is licensed under the Apache License, Version 2.0 
c  (the "License"); you may not use this file except in compliance with the 
c  License. You may obtain a copy of the License at 
c  http://www.apache.org/licenses/LICENSE-2.0. 
c 
c  Unless required by applicable law or agreed to in writing, software 
c  distributed under the License is distributed on an "AS IS" BASIS, WITHOUT 
c  WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the 
c  License for the specific language governing permissions and limitations 
c  under the License.
c  ---------------------------------------------------------------------------
c
      program cgns_to_cfl3dinput
c
c***********************************************************************
c   Purpose: Reads in CGNS file and writes out 3d plot3d-MG grid file,
c   plus an initial-guess CFL3D-type input file.
c   This program is not capable of reading or writing iblanking info.
c   It also assumes that the dimensions (or nondimensionalizations)
c   are already in the CGNS file exactly as they are desired to come
c   out in the PLOT3D file.
c***********************************************************************
c
#ifdef CGNS
c
#     include "cgnslib_f.h"
c
      dimension isize(9)
c
      character*80 file1,file2
      character*32 zonename
c
c  Read CGNS file
      write(6,'('' Note:  currently, you must be linked to CGNS'',
     +  '' V2.5.2 or later'',/)')
      write(6,'('' input name of CGNS file to read'')')
      read(5,'(a80)') file2
      write(6,'('' input 1 if you want to key BCs off of boconame'',
     +  '' instead of the preferred bctype'')')
      write(6,'(''     (looks for body, farfield, and symmetry)'')')
      write(6,'(''     0 is default:'')')
      read(5,*) ibcname_read
      call cg_open_f(file2,CG_MODE_READ,iccg,ier)
      if (ier .ne. 0) call cg_error_exit_f
      call cg_nbases_f(iccg,nbases,ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nbases .ne. 1) then
        write(6,'('' Error nbases ='',i5,'' (should be 1)'')') nbases
        write(6,'('' stopping'')')
        stop
      end if
      ibase=1
c   Get number of zones
      call cg_nzones_f(iccg,ibase,nblock,ier)
      if (ier .ne. 0) call cg_error_exit_f
      write(6,'('' number of zones='',i9)') nblock
c   Get and write out the grid
      call gridtop(iccg,ibase,nblock,file1)
c   Now get 1-to-1 connectivity info
      call cg_n1to1_global_f(iccg,ibase,num1to1,ier)
      if (ier .ne. 0) call cg_error_exit_f
      write(6,'('' total number of 1-to-1 interfaces present='',i9)')
     +  num1to1
      write(6,'(''   (this reflects total number divided by 2)'')')
c
      call reader(iccg,ibase,nblock,num1to1,file1,ibcname_read)
c
c   Close data base
      call cg_close_f(iccg,ier)
      if (ier .ne. 0) call cg_error_exit_f
c
      stop
      end
c
c*******************************************************************
c
      subroutine gridtop(iccg,ibase,nblock,file1)
c
      integer stats
      dimension isize(9)
      character*80 file1
      character*32 zonename
c
      allocatable :: idim(:)
      allocatable :: jdim(:)
      allocatable :: kdim(:)
c   allocate memory
      memuse=0
      allocate( idim(nblock), stat=stats )
      call umalloc_r(nblock,0,'idim',memuse,stats)
      allocate( jdim(nblock), stat=stats )
      call umalloc_r(nblock,0,'jdim',memuse,stats)
      allocate( kdim(nblock), stat=stats )
      call umalloc_r(nblock,0,'kdim',memuse,stats)
c
      do n=1,nblock
c   Read general zone information
        call cg_zone_read_f(iccg, ibase, n, zonename, isize, ier)
        if (ier .ne. 0) call cg_error_exit_f
        idim(n)=isize(1)
        jdim(n)=isize(2)
        kdim(n)=isize(3)
      enddo
c
c  Write out plot3d grid and q files:
      write(6,'('' do you want 0=formatted or 1=unformatted'',
     +  '' grid file output?'')')
      read(5,*) iform
      if (iform .eq. 1) then
c  unformatted-type
      write(6,'('' input name of unformatted grid filename to write'')')
      read(5,'(a80)') file1
      open(2,file=file1,form='unformatted',status='unknown')
        write(2) nblock
        write(2) (idim(n),jdim(n),kdim(n),n=1,nblock)
      write(6,'('' unformatted files written in whatever'',
     + '' precision (single or double)'',/,'' that you compiled'',
     + '' this translator program in'')')
      else
c  formatted-type
      write(6,'('' input name of formatted grid filename to write'')')
      read(5,'(a80)') file1
      open(2,file=file1,form='formatted',status='unknown')
        write(2,*) nblock
        write(2,*) (idim(n),jdim(n),kdim(n),n=1,nblock)
      end if
c
      do n=1,nblock
        call gridbot(n,iccg,ibase,idim(n),jdim(n),kdim(n),iform)
      enddo
c
c   deallocate memory
      deallocate(idim)
      deallocate(jdim)
      deallocate(kdim)
      close(2)
c
      return
      end
c
c*******************************************************************
c
      subroutine gridbot(n,iccg,ibase,idm,jdm,kdm,iform)
c
      integer stats
      dimension isize(9)
      allocatable :: x(:,:,:)
      allocatable :: y(:,:,:)
      allocatable :: z(:,:,:)
c   allocate memory
      memuse=0
      allocate( x(idm,jdm,kdm), stat=stats )
      call umalloc_r(idm*jdm*kdm,0,'x',memuse,stats)
      allocate( y(idm,jdm,kdm), stat=stats )
      call umalloc_r(idm*jdm*kdm,0,'y',memuse,stats)
      allocate( z(idm,jdm,kdm), stat=stats )
      call umalloc_r(idm*jdm*kdm,0,'z',memuse,stats)
c
c  Read CGNS file
      nblock=1
      write(6,'('' for zone '',i9,'', idim,jdim,kdim='',3i5)')
     +    n,idm,jdm,kdm
      call getgrd(iccg,ibase,n,idm,jdm,kdm,x,y,z)
c
c  Write out plot3d grid and q files:
      if (iform .eq. 1) then
c  unformatted-type
        write(6,'('' block#'',i9,'': id,jd,kd='',3i5)') n,idm,jdm,kdm
        write(2) (((x(i,j,k),i=1,idm),j=1,jdm),k=1,kdm),
     +           (((y(i,j,k),i=1,idm),j=1,jdm),k=1,kdm),
     +           (((z(i,j,k),i=1,idm),j=1,jdm),k=1,kdm)
      else
        write(6,'('' block#'',i9,'': id,jd,kd='',3i5)') n,idm,jdm,kdm
        write(2,*) (((x(i,j,k),i=1,idm),j=1,jdm),k=1,kdm),
     +             (((y(i,j,k),i=1,idm),j=1,jdm),k=1,kdm),
     +             (((z(i,j,k),i=1,idm),j=1,jdm),k=1,kdm)
      end if
c
c   deallocate memory
      deallocate(x)
      deallocate(y)
      deallocate(z)
c
      return
      end
c*******************************************************************
c
      subroutine reader(iccg,ibase,nblock,num1to1,file1,ibcname_read)
c
      parameter(numbcperzone=100,nfammax=100)
c
      integer stats
      dimension isize(9)
      allocatable :: idim(:)
      allocatable :: jdim(:)
      allocatable :: kdim(:)
      allocatable :: irange(:,:)
      allocatable :: idonor_range(:,:)
      allocatable :: itransform(:,:)
      allocatable :: nbc(:)
      allocatable :: ibctype(:,:)
      allocatable :: ips(:,:,:)
c
      character*80 file1
      character*32 zonename(nblock)
      character*32 connectname(num1to1),znname(num1to1),
     +  donorname(num1to1)
c
#     include "cgnslib_f.h"
c
      character*32 boconame,familyname
      character*32 famname(nfammax),dumfam
      character*3 boconame_short
      dimension normalindex(3),ipnts(6)
      dimension data_double(6)
      dimension nfambc(nfammax),ibctypef(nfammax)
c
c   allocate memory
      memuse=0
      allocate( idim(nblock), stat=stats )
      call umalloc_r(nblock,0,'idim',memuse,stats)
      allocate( jdim(nblock), stat=stats )
      call umalloc_r(nblock,0,'jdim',memuse,stats)
      allocate( kdim(nblock), stat=stats )
      call umalloc_r(nblock,0,'kdim',memuse,stats)
      allocate( irange(6,num1to1), stat=stats )
      call umalloc_r(6*num1to1,0,'irange',memuse,stats)
      allocate( idonor_range(6,num1to1), stat=stats )
      call umalloc_r(6*num1to1,0,'idonor_range',memuse,stats)
      allocate( itransform(3,num1to1), stat=stats )
      call umalloc_r(3*num1to1,0,'itransform',memuse,stats)
      allocate( nbc(nblock), stat=stats )
      call umalloc_r(nblock,0,'nbc',memuse,stats)
      allocate( ibctype(nblock,numbcperzone), stat=stats )
      call umalloc_r(nblock*numbcperzone,0,'ibctype',memuse,stats)
      allocate( ips(nblock,numbcperzone,6), stat=stats )
      call umalloc_r(nblock*numbcperzone*6,0,'ips',memuse,stats)
c
c  First, see if there are any Family_t nodes
      call cg_nfamilies_f(iccg, ibase, nfamilies, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (nfamilies .gt. nfammax) then
        write(6,'('' Error, need to increase nfammax to at least '',
     +     i5)') nfamilies
        stop
      end if
      if (nfamilies .gt. 0) then
      write(6,'('' Found '',i5,'' families:'')') nfamilies
      do n=1,nfamilies
        call cg_family_read_f(iccg, ibase, n, famname(n), nfambc(n), 
     +     ngeo, ier)
        call cg_fambc_read_f(iccg, ibase, n, 1, dumfam,
     +     ibctypef(n), ier)
        write(6,'('' n,famname(n),type='',i2,3x,a32,3x,a32)') 
     +     n,famname(n),BCTypeName(ibctypef(n))
      enddo
      end if
c
c  Read CGNS file
      nbl=nblock
c   Do loop over the zones
      do n=1,nbl
c   Read general zone information
        call cg_zone_read_f(iccg, ibase, n, zonename(n), isize, ier)
        if (ier .ne. 0) call cg_error_exit_f
        idim(n)=isize(1)
        jdim(n)=isize(2)
        kdim(n)=isize(3)
      enddo
c
c   Now get 1-to-1 connectivity info
      n1to1_global=num1to1
      call cg_1to1_read_global_f(iccg,ibase,connectname,znname,
     +  donorname,irange,idonor_range,itransform,ier)
c  Need the following commented out in case there are zero 1-to-1 connections:
c     if (ier .ne. 0) call cg_error_exit_f
      if (n1to1_global .gt. 0) then
      do n=1,n1to1_global
        write(6,'('' n='',i9,'':'')') n
        write(6,'(''    connectname ='',a32)') connectname(n)
        write(6,'(''    znname      ='',a32)') znname(n)
        write(6,'(''    donorname   ='',a32)') donorname(n)
        write(6,'(''    irange      ='',6i7)') irange(1,n),irange(2,n),
     +   irange(3,n),irange(4,n),irange(5,n),irange(6,n)
        write(6,'(''    idonor_range='',6i7)') idonor_range(1,n),
     +   idonor_range(2,n),idonor_range(3,n),idonor_range(4,n),
     +   idonor_range(5,n),idonor_range(6,n)
        write(6,'(''    itransform  ='',3i7)') itransform(1,n),
     +   itransform(2,n),itransform(3,n)
      enddo
      end if
c
c   See if there is any BC info
      do n=1,nbl
        nbc(n)=0
        call cg_nbocos_f(iccg, ibase, n, nbocos, ier)
        if (ier .ne. 0) call cg_error_exit_f
        write(6,'('' for zone '',i9,'', nbocos='',i9)') n,nbocos
        if (nbocos .gt. 0) then
          do m=1,nbocos
            call cg_boco_info_f(iccg, ibase, n, m, boconame, ibocotype, 
     +        iptset_type, npnts, normalindex, normalflag, 
     +        idatatype, ndataset, ier)
            if (ier .ne. 0) call cg_error_exit_f
            if (npnts .ne. 2) then
              write(6,'('' ERROR, cannot use BCs... npnts must be 2'')')
              goto 333
            end if
            if (iptset_type .ne. PointRange) then
              write(6,'('' ERROR, cannot use BCs... iptset_type must'',
     +         '' be PointRange'')')
              goto 333
            end if
            if (ibcname_read .eq. 1) then
              boconame_short(1:3)=boconame(1:3)
              if (boconame_short .ne. "bod" .and.
     +            boconame_short .ne. "BOD" .and.
     +            boconame_short .ne. "Bod" .and.
     +            boconame_short .ne. "sym" .and.
     +            boconame_short .ne. "SYM" .and.
     +            boconame_short .ne. "Sym" .and.
     +            boconame_short .ne. "far" .and.
     +            boconame_short .ne. "FAR" .and.
     +            boconame_short .ne. "Far") then
                write(6,'('' ERROR:'',a32,'' not recognized'')') 
     +            boconame
                goto 333
              end if
            else
              if (ibocotype .ne. BCFarfield .and.
     +            ibocotype .ne. BCInflow .and.
     +            ibocotype .ne. BCSymmetryPlane .and.
     +            ibocotype .ne. BCWallViscous .and.
     +            ibocotype .ne. BCExtrapolate .and.
     +            ibocotype .ne. BCWallViscousHeatFlux .and.
     +            ibocotype .ne. BCWallViscousIsothermal .and.
     +            ibocotype .ne. BCWallInviscid .and.
     +            ibocotype .ne. BCInflowSupersonic .and.
     +            ibocotype .ne. BCTunnelInflow .and.
     +            ibocotype .ne. BCSymmetryPolar .and.
     +            ibocotype .ne. BCWall .and.
     +            ibocotype .ne. FamilySpecified .and.
     +            ibocotype .ne. BCDegenerateLine) then
                write(6,'('' ERROR: bocotype '',a32,
     +            '' not recognized'')') BCTypeName(ibocotype)
                goto 333
              end if
            end if
            call cg_boco_read_f(iccg, ibase, n, m, ipnts,
     +        data_double, ier)
            if (ier .ne. 0) call cg_error_exit_f
            if ((ipnts(1).eq.ipnts(4)).and.(ipnts(2).eq.ipnts(5)).or.
     +          (ipnts(1).eq.ipnts(4)).and.(ipnts(3).eq.ipnts(6)).or.
     +          (ipnts(2).eq.ipnts(5)).and.(ipnts(3).eq.ipnts(6))) then
              write(6,'('' ERROR: BC points screwed up'')')
              goto 333
            end if
c           write(6,'(''... boconame='',a32)') boconame
c           write(6,'(''... ibocotype='',a32)') BCTypeName(ibocotype)
c           write(6,'(''... iptset_type='',a32)') 
c    +       PointSetTypeName(iptset_type)
c           write(6,'(''... npnts='',i5)') npnts
c           write(6,'(''... normalindex='',3i6)') normalindex(1),
c    +       normalindex(2),normalindex(3)
c           write(6,'(''... normalflag='',i6)') normalflag
c           write(6,'(''... idatatype='',a32)') 
c    +       DataTypeName(idatatype)
c           write(6,'(''... ndataset='',i6)') ndataset
c           write(6,'(''... ipnts='',6i12)') ipnts(1),ipnts(2),ipnts(3),
c    +       ipnts(4),ipnts(5),ipnts(6)
c   the relevant info that I can use is ibocotype and ipnts(1-6)
            nbc(n)=nbc(n)+1
            if (ibcname_read .eq. 1) then
              if (boconame_short .eq. "bod" .or.
     +            boconame_short .eq. "BOD" .or.
     +            boconame_short .eq. "Bod") then
                ibctype(n,nbc(n))=22
              else if (boconame_short .eq. "sym" .or.
     +            boconame_short .eq. "SYM" .or.
     +            boconame_short .eq. "Sym") then
                ibctype(n,nbc(n))=16
              else if (boconame_short .eq. "far" .or.
     +            boconame_short .eq. "FAR" .or.
     +            boconame_short .eq. "Far") then
                ibctype(n,nbc(n))=7
              end if
            else
              if (ibocotype .eq. FamilySpecified) then
                call cg_goto_f(iccg, ibase, ier, 'Zone_t', n, 
     +            'ZoneBC_t', 1, 'BC_t', m, 'end')
                if (ier .ne. 0) call cg_error_exit_f
                call cg_famname_read_f(familyname, ier)
                if (ier .ne. 0) call cg_error_exit_f
                write(6,'('' familyname='',a32)') familyname
                do mmm=1,nfamilies
                  if(familyname .eq. famname(mmm)) then
                    ibctype(n,nbc(n))=ibctypef(mmm)
                  end if
                enddo
              else
                ibctype(n,nbc(n))=ibocotype
              end if
            end if
            do ii=1,6
              ips(n,nbc(n),ii)=ipnts(ii)
            enddo
 333        continue
          enddo
        end if
        write(6,'('' nbc in zone '',i9,'' ='',i7)') n,nbc(n)
        do nn=1,nbc(n)
          write(6,'('' zone='',i9,'', ibctype='',a32,6i7)') n,
     +      BCTypeName(ibctype(n,nn)),ips(n,nn,1),ips(n,nn,2),
     +      ips(n,nn,3),ips(n,nn,4),ips(n,nn,5),ips(n,nn,6)
        enddo
        if (nbc(n) .gt. numbcperzone) then
          write(6,'('' stopping... need to increase numbcperzone'',
     +      '' to '',i7)') nbc(n)
          stop
        end if
      enddo
c
c   Write out an initial-guess CFL3D-type input file, including
c   all 1-to-1 information and an attempt at BC information
      open(7,file='cfl3d_guess.inp',form='formatted',status='unknown')
c
      ncgmin=1000
      do n=1,nbl
        iddd=idim(n)
        jddd=jdim(n)
        kddd=kdim(n)
        ncg=0
        do m=1,4
          if (idim(n) .ne. 2) then
          if (float((iddd+1)/2) .ne. float(iddd+1)/2.) goto 100
          iddd=(iddd+1)/2
          end if
          if (float((jddd+1)/2) .ne. float(jddd+1)/2.) goto 100
          jddd=(jddd+1)/2
          if (float((kddd+1)/2) .ne. float(kddd+1)/2.) goto 100
          kddd=(kddd+1)/2
          ncg=ncg+1
        enddo
 100    continue
        ncgmin=min(ncgmin,ncg)
      enddo
c
      write(7,'(''FILES:'')')
      write(7,'(a80)') file1
      write(7,'(''plot3dg.bin'')')
      write(7,'(''plot3dq.bin'')')
      write(7,'(''cfl3d.out'')')
      write(7,'(''cfl3d.res'')')
      write(7,'(''cfl3d.turres'')')
      write(7,'(''cfl3d.blomax'')')
      write(7,'(''cfl3d.out15'')')
      write(7,'(''cfl3d.prout'')')
      write(7,'(''cfl3d.out20'')')
      write(7,'(''ovrlp.bin'')')
      write(7,'(''patch.bin'')')
      write(7,'(''restart.bin'')')
      write(7,'(''    Title'')')
      write(7,'(''     XMACH     ALPHA      BETA  REUE,MIL   TINF,DR'',
     +  ''     IALPH     IHIST'')')
      write(7,'(''    x.xxxx    x.xxxx    0.0000    x.xxxx  520.0000'',
     +  ''         x         0'')')
      write(7,'(''      SREF      CREF      BREF       XMC       YMC'',
     +  ''       ZMC'')')
      write(7,'(''    x.xxxx    x.xxxx    x.xxxx    0.0000    0.0000'',
     +  ''    0.0000'')')
      write(7,'(''        DT     IREST   IFLAGTS      FMAX     IUNST'',
     +  ''    CFLTAU'')')
      write(7,'(''   -5.0000         0         0    5.0000         0'',
     +  ''   10.0000'')')
      write(7,'(''     NGRID   NPLOT3D    NPRINT    NWREST      ICHK'',
     +  ''       I2D    NTSTEP       ITA'')')
      if (idim(1) .eq. 2) then
      write(7,'(2i10,''         0      9900         0         1'',
     +  ''         1        -2'')') -nbl,nbl
      else
      write(7,'(2i10,''         0      9900         0         0'',
     +  ''         1        -2'')') -nbl,nbl
      end if
      write(7,'(''       NCG       IEM  IADVANCE    IFORCE  IVISC(I)'',
     +  ''  IVISC(J)  IVISC(K)'')')
      do n=1,nbl
      write(7,'(i10,''         0         0       333         5'',
     +  ''         5         5'')') ncgmin
      enddo
      write(7,'(''      IDIM      JDIM      KDIM'')')
      do n=1,nbl
      write(7,'(3i10)') idim(n),jdim(n),kdim(n)
      enddo
      write(7,'(''    ILAMLO    ILAMHI    JLAMLO    JLAMHI    KLAMLO'',
     +  ''    KLAMHI'')')
      do n=1,nbl
      write(7,'(''         0         0         0         0         0'',
     +  ''         0'')')
      enddo
      write(7,'(''     INEWG    IGRIDC        IS        JS        KS'',
     +  ''        IE        JE        KE'')')
      do n=1,nbl
      write(7,'(''         0         0         0         0         0'',
     +  ''         0         0         0'')')
      enddo
      write(7,'(''  IDIAG(I)  IDIAG(J)  IDIAG(K)  IFLIM(I)  IFLIM(J)'',
     +  ''  IFLIM(K)'')')
      do n=1,nbl
      write(7,'(''         1         1         1         4         4'',
     +  ''         4'')')
      enddo
      write(7,'(''   IFDS(I)   IFDS(J)   IFDS(K)  RKAP0(I)  RKAP0(J)'',
     +  ''  RKAP0(K)'')')
      do n=1,nbl
      write(7,'(''         1         1         1   0.33333   0.33333'',
     +  ''   0.33333'')')
      enddo
      write(7,'(''      GRID     NBCI0   NBCIDIM     NBCJ0   NBCJDIM'',
     +  ''     NBCK0   NBCKDIM    IOVRLP'')')
      do n=1,nbl
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(1,m) .eq. irange(4,m) .and.
     +        (irange(1,m) .eq. 1)) then
              iseg=iseg+1
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(1,m) .eq. idonor_range(4,m) .and.
     +        (idonor_range(1,m) .eq. 1)) then
              iseg=iseg+1
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,1).eq.ips(n,nn,4)) .and.
     +          (ips(n,nn,1).eq.1)) then
              iseg=iseg+1
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
          iseg=1
        else
        if (isum .ne. (jdim(n)-1)*(kdim(n)-1)) iseg=iseg+1
        end if
        iloseg=iseg
c
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(1,m) .eq. irange(4,m) .and.
     +        (irange(1,m) .eq. idim(n))) then
              iseg=iseg+1
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(1,m) .eq. idonor_range(4,m) .and.
     +        (idonor_range(1,m) .eq. idim(n))) then
              iseg=iseg+1
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,1).eq.ips(n,nn,4)) .and.
     +          (ips(n,nn,1).ne.1)) then
              iseg=iseg+1
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
          iseg=1
        else
        if (isum .ne. (jdim(n)-1)*(kdim(n)-1)) iseg=iseg+1
        end if
        ihiseg=iseg
c
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(2,m) .eq. irange(5,m) .and.
     +        (irange(2,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(2,m) .eq. idonor_range(5,m) .and.
     +        (idonor_range(2,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,2).eq.ips(n,nn,5)) .and.
     +          (ips(n,nn,2).eq.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
          iseg=1
        else
        if (isum .ne. (idim(n)-1)*(kdim(n)-1)) iseg=iseg+1
        end if
        jloseg=iseg
c
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(2,m) .eq. irange(5,m) .and.
     +        (irange(2,m) .eq. jdim(n))) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(2,m) .eq. idonor_range(5,m) .and.
     +        (idonor_range(2,m) .eq. jdim(n))) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,2).eq.ips(n,nn,5)) .and.
     +          (ips(n,nn,2).ne.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
          iseg=1
        else
        if (isum .ne. (idim(n)-1)*(kdim(n)-1)) iseg=iseg+1
        end if
        jhiseg=iseg
c
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(3,m) .eq. irange(6,m) .and.
     +        (irange(3,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(3,m) .eq. idonor_range(6,m) .and.
     +        (idonor_range(3,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,3).eq.ips(n,nn,6)) .and.
     +          (ips(n,nn,3).eq.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
          iseg=1
        else
        if (isum .ne. (idim(n)-1)*(jdim(n)-1)) iseg=iseg+1
        end if
        kloseg=iseg
c
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(3,m) .eq. irange(6,m) .and.
     +        (irange(3,m) .eq. kdim(n))) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(3,m) .eq. idonor_range(6,m) .and.
     +        (idonor_range(3,m) .eq. kdim(n))) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,3).eq.ips(n,nn,6)) .and.
     +          (ips(n,nn,3).ne.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
          iseg=1
        else
        if (isum .ne. (idim(n)-1)*(jdim(n)-1)) iseg=iseg+1
        end if
        khiseg=iseg
c
      write(7,'(7i10,''         0'')') n,iloseg,ihiseg,jloseg,jhiseg,
     +  kloseg,khiseg
      enddo
      write(7,'(''I0:   GRID   SEGMENT    BCTYPE      JSTA      JEND'',
     +  ''      KSTA      KEND     NDATA'')')
      do n=1,nbl
        ibctyp=0
        ndata=0
c       search thru all the 1-to-1 connections
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(1,m) .eq. irange(4,m) .and. 
     +        (irange(1,m) .eq. 1)) then
              iseg=iseg+1
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,kend,ndata
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(1,m) .eq. idonor_range(4,m) .and.
     +        (idonor_range(1,m) .eq. 1)) then
              iseg=iseg+1
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,kend,ndata
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,1).eq.ips(n,nn,4)) .and. 
     +          (ips(n,nn,1).eq.1)) then
              iseg=iseg+1
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              if (ibctype(n,nn) .eq. BCFarfield .or.
     +            ibctype(n,nn) .eq. BCInflow) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPlane) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCExtrapolate) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallInviscid) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCInflowSupersonic) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCTunnelInflow) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPolar) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCDegenerateLine) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallViscousIsothermal) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                 -1.        0.'')')
              else if (ibctype(n,nn) .eq. BCWallViscous .or.
     +         ibctype(n,nn) .eq. BCWallViscousHeatFlux .or.
     +         ibctype(n,nn) .eq. BCWall) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                  0.        0.'')')
              end if
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
        write(7,'(i10,''         1      xxxx         0         0'',
     +    ''         0         0         0'')') n
        else
        if (isum .ne. (jdim(n)-1)*(kdim(n)-1)) then
        write(7,'(2i10,''      xxxx        xx        xx'',
     +    ''        xx        xx         0'')') n,iseg+1
        end if
        end if
      enddo
      write(7,'(''IDIM: GRID   SEGMENT    BCTYPE      JSTA      JEND'',
     +  ''      KSTA      KEND     NDATA'')')
      do n=1,nbl
        ibctyp=0
        ndata=0
c       search thru all the 1-to-1 connections
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(1,m) .eq. irange(4,m) .and.
     +        (irange(1,m) .eq. idim(n))) then
              iseg=iseg+1
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,kend,ndata
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(1,m) .eq. idonor_range(4,m) .and.
     +        (idonor_range(1,m) .eq. idim(n))) then
              iseg=iseg+1
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,kend,ndata
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,1).eq.ips(n,nn,4)) .and.
     +          (ips(n,nn,1).ne.1)) then
              iseg=iseg+1
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              if (ibctype(n,nn) .eq. BCFarfield .or.
     +            ibctype(n,nn) .eq. BCInflow) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPlane) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCExtrapolate) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallInviscid) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCInflowSupersonic) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCTunnelInflow) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPolar) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCDegenerateLine) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallViscousIsothermal) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                 -1.        0.'')')
              else if (ibctype(n,nn) .eq. BCWallViscous .or.
     +         ibctype(n,nn) .eq. BCWallViscousHeatFlux .or.
     +         ibctype(n,nn) .eq. BCWall) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                  0.        0.'')')
              end if
              isum=isum+((jend-jsta)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
        write(7,'(i10,''         1      xxxx         0         0'',
     +    ''         0         0         0'')') n
        else
        if (isum .ne. (jdim(n)-1)*(kdim(n)-1)) then
        write(7,'(2i10,''      xxxx        xx        xx'',
     +    ''        xx        xx         0'')') n,iseg+1
        end if
        end if
      enddo
      write(7,'(''J0:   GRID   SEGMENT    BCTYPE      ISTA      IEND'',
     +  ''      KSTA      KEND     NDATA'')')
      do n=1,nbl
        ibctyp=0
        ndata=0
c       search thru all the 1-to-1 connections
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(2,m) .eq. irange(5,m) .and.
     +        (irange(2,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,kend,ndata
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(2,m) .eq. idonor_range(5,m) .and.
     +        (idonor_range(2,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,kend,ndata
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,2).eq.ips(n,nn,5)) .and.
     +          (ips(n,nn,2).eq.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              if (ibctype(n,nn) .eq. BCFarfield .or.
     +            ibctype(n,nn) .eq. BCInflow) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPlane) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCExtrapolate) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallInviscid) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCInflowSupersonic) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCTunnelInflow) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPolar) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCDegenerateLine) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallViscousIsothermal) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                 -1.        0.'')')
              else if (ibctype(n,nn) .eq. BCWallViscous .or.
     +         ibctype(n,nn) .eq. BCWallViscousHeatFlux .or.
     +         ibctype(n,nn) .eq. BCWall) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                  0.        0.'')')
              end if
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
        write(7,'(i10,''         1      xxxx         0         0'',
     +    ''         0         0         0'')') n
        else
        if (isum .ne. (idim(n)-1)*(kdim(n)-1)) then
        write(7,'(2i10,''      xxxx        xx        xx'',
     +    ''        xx        xx         0'')') n,iseg+1
        end if
        end if
      enddo
      write(7,'(''JDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND'',
     +  ''      KSTA      KEND     NDATA'')')
      do n=1,nbl
        ibctyp=0
        ndata=0
c       search thru all the 1-to-1 connections
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(2,m) .eq. irange(5,m) .and.
     +        (irange(2,m) .eq. jdim(n))) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              ksta=min(irange(3,m),irange(6,m))
              kend=max(irange(3,m),irange(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,kend,ndata
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(2,m) .eq. idonor_range(5,m) .and.
     +        (idonor_range(2,m) .eq. jdim(n))) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              ksta=min(idonor_range(3,m),idonor_range(6,m))
              kend=max(idonor_range(3,m),idonor_range(6,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,kend,ndata
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,2).eq.ips(n,nn,5)) .and.
     +          (ips(n,nn,2).ne.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              ksta=min(ips(n,nn,3),ips(n,nn,6))
              kend=max(ips(n,nn,3),ips(n,nn,6))
              if (ibctype(n,nn) .eq. BCFarfield .or.
     +            ibctype(n,nn) .eq. BCInflow) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPlane) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCExtrapolate) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallInviscid) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCInflowSupersonic) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCTunnelInflow) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPolar) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCDegenerateLine) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. BCWallViscousIsothermal) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                 -1.        0.'')')
              else if (ibctype(n,nn) .eq. BCWallViscous .or.
     +         ibctype(n,nn) .eq. BCWallViscousHeatFlux .or.
     +         ibctype(n,nn) .eq. BCWall) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                  0.        0.'')')
              end if
              isum=isum+((iend-ista)*(kend-ksta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
        write(7,'(i10,''         1      xxxx         0         0'',
     +    ''         0         0         0'')') n
        else
        if (isum .ne. (idim(n)-1)*(kdim(n)-1)) then
        write(7,'(2i10,''      xxxx        xx        xx'',
     +    ''        xx        xx         0'')') n,iseg+1
        end if
        end if
      enddo
      write(7,'(''K0:   GRID   SEGMENT    BCTYPE      ISTA      IEND'',
     +  ''      JSTA      JEND     NDATA'')')
      do n=1,nbl
        ibctyp=0
        ndata=0
c       search thru all the 1-to-1 connections
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(3,m) .eq. irange(6,m) .and.
     +        (irange(3,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,jend,ndata
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(3,m) .eq. idonor_range(6,m) .and.
     +        (idonor_range(3,m) .eq. 1)) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,jend,ndata
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,3).eq.ips(n,nn,6)) .and.
     +          (ips(n,nn,3).eq.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              if (ibctype(n,nn) .eq. BCFarfield .or.
     +            ibctype(n,nn) .eq. BCInflow) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPlane) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCExtrapolate) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCWallInviscid) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCInflowSupersonic) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCTunnelInflow) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPolar) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCDegenerateLine) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCWallViscousIsothermal) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                 -1.        0.'')')
              else if (ibctype(n,nn) .eq. BCWallViscous .or.
     +         ibctype(n,nn) .eq. BCWallViscousHeatFlux .or.
     +         ibctype(n,nn) .eq. BCWall) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                  0.        0.'')')
              end if
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
        write(7,'(i10,''         1      xxxx         0         0'',
     +    ''         0         0         0'')') n
        else
        if (isum .ne. (idim(n)-1)*(jdim(n)-1)) then
        write(7,'(2i10,''      xxxx        xx        xx'',
     +    ''        xx        xx         0'')') n,iseg+1
        end if
        end if
      enddo
      write(7,'(''KDIM: GRID   SEGMENT    BCTYPE      ISTA      IEND'',
     +  ''      JSTA      JEND     NDATA'')')
      do n=1,nbl
        ibctyp=0
        ndata=0
c       search thru all the 1-to-1 connections
        iseg=0
        isum=0
        do m=1,n1to1_global
          if (znname(m) .eq. zonename(n)) then
            if (irange(3,m) .eq. irange(6,m) .and.
     +        (irange(3,m) .eq. kdim(n))) then
              iseg=iseg+1
              ista=min(irange(1,m),irange(4,m))
              iend=max(irange(1,m),irange(4,m))
              jsta=min(irange(2,m),irange(5,m))
              jend=max(irange(2,m),irange(5,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,jend,ndata
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
          if (donorname(m) .eq. zonename(n)) then
            if (idonor_range(3,m) .eq. idonor_range(6,m) .and.
     +        (idonor_range(3,m) .eq. kdim(n))) then
              iseg=iseg+1
              ista=min(idonor_range(1,m),idonor_range(4,m))
              iend=max(idonor_range(1,m),idonor_range(4,m))
              jsta=min(idonor_range(2,m),idonor_range(5,m))
              jend=max(idonor_range(2,m),idonor_range(5,m))
              write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,jend,ndata
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          end if
        enddo
        if (nbc(n) .gt. 0) then
          do nn=1,nbc(n)
            if ((ips(n,nn,3).eq.ips(n,nn,6)) .and.
     +          (ips(n,nn,3).ne.1)) then
              iseg=iseg+1
              ista=min(ips(n,nn,1),ips(n,nn,4))
              iend=max(ips(n,nn,1),ips(n,nn,4))
              jsta=min(ips(n,nn,2),ips(n,nn,5))
              jend=max(ips(n,nn,2),ips(n,nn,5))
              if (ibctype(n,nn) .eq. BCFarfield .or.
     +            ibctype(n,nn) .eq. BCInflow) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPlane) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCExtrapolate) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCWallInviscid) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCInflowSupersonic) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCTunnelInflow) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCSymmetryPolar) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCDegenerateLine) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. BCWallViscousIsothermal) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                 -1.        0.'')')
              else if (ibctype(n,nn) .eq. BCWallViscous .or.
     +         ibctype(n,nn) .eq. BCWallViscousHeatFlux .or.
     +         ibctype(n,nn) .eq. BCWall) then
                ibctyp=2004
                ndata=2
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
                write(7,'(''              TWTYPE        CQ'')')
                write(7,'(''                  0.        0.'')')
              end if
              isum=isum+((iend-ista)*(jend-jsta))
            end if
          enddo
        end if
        if (iseg .eq. 0) then
        write(7,'(i10,''         1      xxxx         0         0'',
     +    ''         0         0         0'')') n
        else
        if (isum .ne. (idim(n)-1)*(jdim(n)-1)) then
        write(7,'(2i10,''      xxxx        xx        xx'',
     +    ''        xx        xx         0'')') n,iseg+1
        end if
        end if
      enddo
      write(7,'(''      MSEQ    MGFLAG    ICONSF       MTT'',
     +  ''      NGAM'')')
      if (ncgmin .eq. 0) then
      write(7,'(''         1         0         0         0'',
     +  ''         2'')')
      else
      write(7,'(''         1         1         0         0'',
     +  ''         2'')')
      end if
      write(7,'(''      ISSC EPSSSC(1) EPSSSC(2) EPSSSC(3)      ISSR'',
     +  '' EPSSSR(1) EPSSSR(2) EPSSSR(3)'')')
      write(7,'(''         0    0.3000    0.3000    0.3000         0'',
     +  ''    0.3000    0.3000    0.3000'')')
      write(7,'(''      NCYC    MGLEVG     NEMGL     NITFO'')')
      write(7,'(''      xxxx'',i10,''         0         0'')') 
     +  min(ncgmin+1,3)
      write(7,'(''      MIT1      MIT2      MIT3      MIT4      MIT5'',
     +  ''      MIT6      MIT7     MIT8'')')
      write(7,'(''         1         1         1         1         1'',
     +  ''         1         1        1'')')
      write(7,'(''   1-1 BLOCKING DATA:'')')
      write(7,'(''      NBLI'')')
      write(7,'(i10)') n1to1_global
      write(7,'('' NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND'',
     +  ''   JEND   KEND  ISVA1  ISVA2'')')
      do m=1,n1to1_global
        igrid=0
        do n=1,nbl
          if (znname(m) .eq. zonename(n)) igrid=n
        enddo
        if (igrid .eq. 0) then
          write(6,'('' error... zone not found'')')
          stop
        end if
        if (irange(1,m) .eq. irange(4,m)) then
          isva1=2
          isva2=3
        else if (irange(2,m) .eq. irange(5,m)) then
          isva1=1
          isva2=3
        else if (irange(3,m) .eq. irange(6,m)) then
          isva1=1
          isva2=2
        end if
        write(7,'(2i7,7x,8i7)') m,igrid,irange(1,m),irange(2,m),
     +   irange(3,m),irange(4,m),irange(5,m),irange(6,m),isva1,
     +   isva2
      enddo
      write(7,'('' NUMBER   GRID     :    ISTA   JSTA   KSTA   IEND'',
     +  ''   JEND   KEND  ISVA1  ISVA2'')')
      do m=1,n1to1_global
        igrid=0
        do n=1,nbl
          if (donorname(m) .eq. zonename(n)) igrid=n
        enddo
        if (igrid .eq. 0) then
          write(6,'('' error... zone not found'')')
          stop
        end if
        if (idonor_range(1,m) .eq. idonor_range(4,m)) then
          if (abs(itransform(1,m)) .eq. 1) then
            if (abs(itransform(2,m)).lt.abs(itransform(3,m))) then
            isva1=2
            isva2=3
            else
            isva1=3
            isva2=2
            end if
          else if (abs(itransform(2,m)) .eq. 1) then
            if (abs(itransform(1,m)).lt.abs(itransform(3,m))) then
            isva1=2
            isva2=3
            else
            isva1=3
            isva2=2
            end if
          else if (abs(itransform(3,m)) .eq. 1) then
            if (abs(itransform(1,m)).lt.abs(itransform(2,m))) then
            isva1=2
            isva2=3
            else
            isva1=3
            isva2=2
            end if
          end if
        else if (idonor_range(2,m) .eq. idonor_range(5,m)) then
          if (abs(itransform(1,m)) .eq. 2) then
            if (abs(itransform(2,m)).lt.abs(itransform(3,m))) then
            isva1=1
            isva2=3
            else
            isva1=3
            isva2=1
            end if
          else if (abs(itransform(2,m)) .eq. 2) then
            if (abs(itransform(1,m)).lt.abs(itransform(3,m))) then
            isva1=1
            isva2=3
            else
            isva1=3
            isva2=1
            end if
          else if (abs(itransform(3,m)) .eq. 2) then
            if (abs(itransform(1,m)).lt.abs(itransform(2,m))) then
            isva1=1
            isva2=3
            else
            isva1=3
            isva2=1
            end if
          end if
        else if (idonor_range(3,m) .eq. idonor_range(6,m)) then
          if (abs(itransform(1,m)) .eq. 3) then
            if (abs(itransform(2,m)).lt.abs(itransform(3,m))) then
            isva1=1
            isva2=2
            else
            isva1=2
            isva2=1
            end if
          else if (abs(itransform(2,m)) .eq. 3) then
            if (abs(itransform(1,m)).lt.abs(itransform(3,m))) then
            isva1=1
            isva2=2
            else
            isva1=2
            isva2=1
            end if
          else if (abs(itransform(3,m)) .eq. 3) then
            if (abs(itransform(1,m)).lt.abs(itransform(2,m))) then
            isva1=1
            isva2=2
            else
            isva1=2
            isva2=1
            end if
          end if
        end if
        write(7,'(2i7,7x,8i7)') m,igrid,idonor_range(1,m),
     +   idonor_range(2,m),idonor_range(3,m),idonor_range(4,m),
     +   idonor_range(5,m),idonor_range(6,m),isva1,isva2
      enddo
      write(7,'(''  PATCH SURFACE DATA:'')')
      write(7,'(''    NINTER'')')
      write(7,'(''         0'')')
      write(7,'(''  PLOT3D OUTPUT:'')')
      write(7,'(''  BLOCK IPTYPE ISTART   IEND   IINC JSTART   JEND'',
     +  ''   JINC KSTART   KEND   KINC'')')
      do n=1,nbl
      write(7,'(i7,''      0      0      0      0      0      0'',
     +  ''      0      0      0      0'')') n
      enddo
      write(7,'(''  MOVIE'')')
      write(7,'(''      0'')')
      write(7,'(''  PRINT OUT:'')')
      write(7,'(''  BLOCK IPTYPE ISTART   IEND   IINC JSTART   JEND'',
     +  ''   JINC KSTART   KEND   KINC'')')
      write(7,'(''  CONTROL SURFACE:'')')
      write(7,'(''  NCS'')')
      write(7,'(''    0'')')
      write(7,'(''   GRID ISTART   IEND   JSTART   JEND   KSTART'',
     +  ''   KEND  IWALL  INORM'')')
c
      write(6,'(/,'' PLOT3D file is 3-D, MG type (no iblank)'')')
      write(6,'(''     name='',a80)') file1
      write(6,'('' Initial-guess CFL3D-type input file output'',
     +  '' to cfl3d_guess.inp'')')
      write(6,'(''     (includes all 1-to-1 info from CGNS file'')')
      write(6,'(''      and an attempt at the BCs, if available)'')')
      write(6,'(''     Edit cfl3d_guess.inp and replace all'',
     +  '' occurrences of xxxx and xx'')')
c
c   deallocate memory
      deallocate(idim)
      deallocate(jdim)
      deallocate(kdim)
      deallocate(irange)
      deallocate(idonor_range)
      deallocate(itransform)
c
      return
      end
c
c*******************************************************************
c
      subroutine getgrd(iccg,ibase,igrid,idim,jdim,kdim,x,y,z)
c
c   Gets x,y, and z from CGNS data base for zone number "igrid".
c   The CGNS file must already be opened.
c   If there is a grid-size inconsistency with what is expected (idim,jdim,kdim)
c   or if the 3 coordinates CoordinateX, CoordinateY, and CoordinateZ 
c   do not exist, the routine stops execution.
c
c   INPUTS:
c      iccg.............CGNS file index number (determined outside this routine) (integer)
c      ibase............CGNS base index number (determined outside this routine) (integer)
c      igrid............zone (or grid) number (integer)
c      idim,jdim,kdim...expected dimensions of this zone (zone "igrid") (integers)
c   OUTPUTS
c      x,y,z............coordinates, output in (i,j,k) order (real)
c
#     include "cgnslib_f.h"
c
      dimension x(idim,jdim,kdim),y(idim,jdim,kdim),
     +          z(idim,jdim,kdim)
      dimension isize(3*3),irmin(3),irmax(3)
c
      character*32 zonename,coordname,testname(3)
c
c   Determine if single or double precision is being used:
c     idouble=0
c#ifdef DOUBLE
      idouble=1
c#endif
c
c   Coordinate names that we are looking for:
        testname(1)='CoordinateX'
        testname(2)='CoordinateY'
        testname(3)='CoordinateZ'
c
c   Read general zone information
      call cg_zone_read_f(iccg, ibase, igrid, zonename, isize, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if(isize(1) .ne. idim .or. isize(2) .ne. jdim .or.
     +   isize(3) .ne. kdim) then
        write(6,'('' Grid index inconsistencies:  isize='',3i5,
     +   ''idim,jdim,kdim='',3i5)') isize(1),isize(2),isize(3),
     +   idim,jdim,kdim
        write(6,'('' Be sure to order the zones alphabetically'',
     +   '' in the input file!'')')
        stop
      end if
c   Find out how many grid coordinates exist
      call cg_ncoords_f(iccg, ibase, igrid, ncoords, ier)
      if (ier .ne. 0) call cg_error_exit_f
      if (ncoords .ne. 3) then
        write(6,'('' ncoords='',i5,''.  Expecting 3.'')') ncoords
        stop
      end if
c   Check coordinate names in data base
      do icoord=1,3
        call cg_coord_info_f(iccg, ibase, igrid, icoord, itype,
     +    coordname,ier)
        if (ier .ne. 0) call cg_error_exit_f
        if(coordname .eq. testname(1) .or.
     +     coordname .eq. testname(2) .or.
     +     coordname .eq. testname(3)) then
          continue
        else
          write(6,'('' coordname of '',a32,'' unrecognized.'')') 
     +     coordname
          write(6,'('' Looking for CoordinateX, CoordinateY, and'',
     +     '' CoordinateZ'')')
          stop
        end if
      enddo
c   Set up array bounds:
      irmin(1)=1
      irmin(2)=1
      irmin(3)=1
      irmax(1)=idim
      irmax(2)=jdim
      irmax(3)=kdim
c   Read x,y,z
      if (idouble .eq. 1) then
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateX', 
     &                     RealDouble, irmin, irmax, x, ier)
      else
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateX', 
     &                     RealSingle, irmin, irmax, x, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateY', 
     &                     RealDouble, irmin, irmax, y, ier)
      else
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateY', 
     &                     RealSingle, irmin, irmax, y, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
      if (idouble .eq. 1) then
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateZ', 
     &                     RealDouble, irmin, irmax, z, ier)
      else
      call cg_coord_read_f(iccg, ibase, igrid, 'CoordinateZ', 
     &                     RealSingle, irmin, irmax, z, ier)
      end if
      if (ier .ne. 0) call cg_error_exit_f
c
      return
      end
#else
c     this is now just a dummy routine since CGNS libs have not been
c     used
c
      write(6,*)'This code is non-functional since the installation'
      write(6,*)'of cfl3d was done without cgns libraries'
c
      stop
      end
#endif
