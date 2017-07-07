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
      program nmf_to_cfl3dinput
c
c***********************************************************************
c   Purpose: Reads in nmf file 
c   (.nmf=Neutral Map File) - see http://geolab.larc.nasa.gov/Volume/Doc/nmf.htm
c   and writes out an initial-guess CFL3D-type input file.
c***********************************************************************
c
      parameter(nn=1000000)
c
      character*80 file1
      character*10 namex
c
c  Read nmf file
      write(6,'('' input name of nmf file to read'')')
      write(6,'('' NOTE: names in .nmf file cannot be longer than 20'',
     +  '' characters (incl quotes)'')')
      write(6,'('' and numbers MUST start after the 20th column!'')')
      read(5,'(a80)') file1
      open(3,file=file1,form='formatted',status='old')
c
c   Get number of zones
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*) nblock
      read(3,*)
c
      do n=1,nblock
        read(3,*) nnn,idum,jdum,kdum
      enddo
c
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
c
      num1to1=0
      numbocos=0
      do n=1,nn
        read(3,'(a10)',end=1010) namex
        if (namex .eq. "ONE_TO_ONE" .or. namex .eq. "one_to_one") then
          num1to1=num1to1+1
        else
          numbocos=numbocos+1
        end if
      enddo
      write(6,'('' Error, need to increase nn'')')
      stop
 1010 continue
c
      write(6,'('' number of zones='',i6)') nblock
      write(6,'('' number of BC patches (not including 1-to-1:'',i6)')
     +  numbocos
      write(6,'('' total number of 1-to-1 interfaces present='',i6)')
     +  num1to1
      write(6,'(''   (this reflects total number divided by 2)'')')
c
      call reader(nblock,num1to1,numbocos)
c
      write(6,'('' successful completion'')')
      stop
      end
c
c*******************************************************************
c
      subroutine reader(nblock,num1to1,numbocos)
c
      parameter(numbcperzone=100)
c
      integer stats
      allocatable :: idim(:)
      allocatable :: jdim(:)
      allocatable :: kdim(:)
      allocatable :: irange(:,:)
      allocatable :: idonor_range(:,:)
      allocatable :: nbc(:)
      allocatable :: ibctype(:,:)
      allocatable :: ips(:,:,:)
c
      dimension i_one_to(num1to1),i_one_from(num1to1)
c
      character*20 nameo,nameb
      character*10 namex
      character*5 rev(num1to1)
c
      dimension normalindex(3),ipnts(6)
      dimension data_double(6)
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
      allocate( nbc(nblock), stat=stats )
      call umalloc_r(nblock,0,'nbc',memuse,stats)
      allocate( ibctype(nblock,numbcperzone), stat=stats )
      call umalloc_r(nblock*numbcperzone,0,'ibctype',memuse,stats)
      allocate( ips(nblock,numbcperzone,6), stat=stats )
      call umalloc_r(nblock*numbcperzone*6,0,'ips',memuse,stats)
c
      nbc(:)=0
c
c  Read nmf file again
      rewind 3
c
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*) nbl
      read(3,*)
      do n=1,nbl
        read(3,*) nnn,idim(n),jdim(n),kdim(n)
      enddo
      read(3,*)
      read(3,*)
      read(3,*)
      read(3,*)
c
      n_o=0
      n_b=0
      do n=1,1000000
        read(3,'(a10)',end=1010) namex
c       write(6,'('' namex='',a10)') namex
        backspace 3
        if (namex .eq. "ONE_TO_ONE" .or. 
     +      namex .eq. "one_to_one") then
          n_o=n_o+1
          read(3,*) nameo,i_one_to(n_o),if1,ilo1,ihi1,
     +      jlo1,jhi1,
     +      i_one_from(n_o),if2,ilo2,ihi2,jlo2,
     +      jhi2,rev(n_o)
          if (if1 .eq. 1) then
            irange(1,n_o)=ilo1
            irange(2,n_o)=jlo1
            irange(3,n_o)=1
            irange(4,n_o)=ihi1
            irange(5,n_o)=jhi1
            irange(6,n_o)=1
          else if (if1 .eq. 2) then
            irange(1,n_o)=ilo1
            irange(2,n_o)=jlo1
            irange(3,n_o)=kdim(i_one_to(n_o))
            irange(4,n_o)=ihi1
            irange(5,n_o)=jhi1
            irange(6,n_o)=kdim(i_one_to(n_o))
          else if (if1 .eq. 3) then
            irange(1,n_o)=1
            irange(2,n_o)=ilo1
            irange(3,n_o)=jlo1
            irange(4,n_o)=1
            irange(5,n_o)=ihi1
            irange(6,n_o)=jhi1
          else if (if1 .eq. 4) then
            irange(1,n_o)=idim(i_one_to(n_o))
            irange(2,n_o)=ilo1
            irange(3,n_o)=jlo1
            irange(4,n_o)=idim(i_one_to(n_o))
            irange(5,n_o)=ihi1
            irange(6,n_o)=jhi1
          else if (if1 .eq. 5) then
            irange(1,n_o)=jlo1
            irange(2,n_o)=1
            irange(3,n_o)=ilo1
            irange(4,n_o)=jhi1
            irange(5,n_o)=1
            irange(6,n_o)=ihi1
          else
            irange(1,n_o)=jlo1
            irange(2,n_o)=jdim(i_one_to(n_o))
            irange(3,n_o)=ilo1
            irange(4,n_o)=jhi1
            irange(5,n_o)=jdim(i_one_to(n_o))
            irange(6,n_o)=ihi1
          end if
          if (if2 .eq. 1) then
            idonor_range(1,n_o)=ilo2
            idonor_range(2,n_o)=jlo2
            idonor_range(3,n_o)=1
            idonor_range(4,n_o)=ihi2
            idonor_range(5,n_o)=jhi2
            idonor_range(6,n_o)=1
          else if (if2 .eq. 2) then
            idonor_range(1,n_o)=ilo2
            idonor_range(2,n_o)=jlo2
            idonor_range(3,n_o)=kdim(i_one_from(n_o))
            idonor_range(4,n_o)=ihi2
            idonor_range(5,n_o)=jhi2
            idonor_range(6,n_o)=kdim(i_one_from(n_o))
          else if (if2 .eq. 3) then
            idonor_range(1,n_o)=1
            idonor_range(2,n_o)=ilo2
            idonor_range(3,n_o)=jlo2
            idonor_range(4,n_o)=1
            idonor_range(5,n_o)=ihi2
            idonor_range(6,n_o)=jhi2
          else if (if2 .eq. 4) then
            idonor_range(1,n_o)=idim(i_one_from(n_o))
            idonor_range(2,n_o)=ilo2
            idonor_range(3,n_o)=jlo2
            idonor_range(4,n_o)=idim(i_one_from(n_o))
            idonor_range(5,n_o)=ihi2
            idonor_range(6,n_o)=jhi2
          else if (if2 .eq. 5) then
            idonor_range(1,n_o)=jlo2
            idonor_range(2,n_o)=1
            idonor_range(3,n_o)=ilo2
            idonor_range(4,n_o)=jhi2
            idonor_range(5,n_o)=1
            idonor_range(6,n_o)=ihi2
          else
            idonor_range(1,n_o)=jlo2
            idonor_range(2,n_o)=jdim(i_one_from(n_o))
            idonor_range(3,n_o)=ilo2
            idonor_range(4,n_o)=jhi2
            idonor_range(5,n_o)=jdim(i_one_from(n_o))
            idonor_range(6,n_o)=ihi2
          end if
        else
          n_b=n_b+1
          read(3,*) nameb,ig,if,ilo,ihi,jlo,jhi
          nbc(ig)=nbc(ig)+1
          if (nbc(ig) .gt. numbcperzone) then
            write(6,'('' Error. numbcperzone not big enough.'')')
            stop
          end if
          if (nameb .eq. "WALL" .or. nameb .eq. "wall" .or.
     +        nameb .eq. "VISCOUS_SOLID" .or. 
     +        nameb .eq. "viscous_solid") then
            ibctype(ig,nbc(ig))=2004
          else if (nameb .eq. "INVISCID" .or. 
     +             nameb .eq. "inviscid" ) then 
            ibctype(ig,nbc(ig))=1005
          else if (nameb .eq. "FARFIELD" .or. 
     +             nameb .eq. "farfield" .or.
     +             nameb .eq. "INFLOW"   .or.
     +             nameb .eq. "OUTFLOW" ) then
            ibctype(ig,nbc(ig))=1003
          else if (nameb .eq. "SYMMETRY" .or.
     +             nameb .eq. "symmetry" .or.
     +             nameb .eq. "SYMMETRY-X" .or.
     +             nameb .eq. "SYMMETRY_X" .or.
     +             nameb .eq. "SYMMETRY_X_STRONG" .or.
     +             nameb .eq. "symmetry-x" .or.
     +             nameb .eq. "symmetry_x" .or.
     +             nameb .eq. "symmetry_x_strong" .or.
     +             nameb .eq. "SYMMETRY-Y" .or.
     +             nameb .eq. "SYMMETRY_Y" .or.
     +             nameb .eq. "SYMMETRY_Y_STRONG" .or.
     +             nameb .eq. "symmetry-y" .or.
     +             nameb .eq. "symmetry_y" .or.
     +             nameb .eq. "symmetry_y_strong" .or.
     +             nameb .eq. "SYMMETRY-Z" .or.
     +             nameb .eq. "SYMMETRY_Z" .or.
     +             nameb .eq. "SYMMETRY_Z_STRONG" .or.
     +             nameb .eq. "symmetry-z" .or.
     +             nameb .eq. "symmetry_z" .or.
     +             nameb .eq. "symmetry_z_strong") then
            ibctype(ig,nbc(ig))=1001
          else if (nameb .eq. "BACK_PRESSURE" .or.
     +             nameb .eq. "back_pressure") then
            ibctype(ig,nbc(ig))=1002
          else if (nameb .eq. "SUBSONIC_INFLOW_PT" .or.
     +             nameb .eq. "subsonic_inflow_pt") then
            ibctype(ig,nbc(ig))=2010
          else if (nameb .eq. "BC72" ) then
            ibctype(ig,nbc(ig))=1011
          else if (nameb .eq. "BC71") then
            ibctype(ig,nbc(ig))=1013
          else if (nameb .eq. "POLE" .or.
     +             nameb .eq. "pole") then
            write(6,'('' Pole detected.  Input:'')')
            write(6,'(''    1011 = half-plane'')')
            write(6,'(''    1012 = full-plane'')')
            write(6,'(''    1013 = partial-plane:'')')
            read(5,*) i_pole
            ibctype(ig,nbc(ig))=i_pole
          else
            write(6,'('' BC type not known'')')
            write(6,'(a16)') nameb
            stop
          end if
          if (if .eq. 1) then
            ips(ig,nbc(ig),1) = ilo
            ips(ig,nbc(ig),2) = jlo
            ips(ig,nbc(ig),3) = 1
            ips(ig,nbc(ig),4) = ihi
            ips(ig,nbc(ig),5) = jhi
            ips(ig,nbc(ig),6) = 1
          else if (if .eq. 2) then
            ips(ig,nbc(ig),1) = ilo
            ips(ig,nbc(ig),2) = jlo
            ips(ig,nbc(ig),3) = kdim(ig)
            ips(ig,nbc(ig),4) = ihi
            ips(ig,nbc(ig),5) = jhi
            ips(ig,nbc(ig),6) = kdim(ig)
          else if (if .eq. 3) then
            ips(ig,nbc(ig),1) = 1
            ips(ig,nbc(ig),2) = ilo
            ips(ig,nbc(ig),3) = jlo
            ips(ig,nbc(ig),4) = 1
            ips(ig,nbc(ig),5) = ihi
            ips(ig,nbc(ig),6) = jhi
          else if (if .eq. 4) then
            ips(ig,nbc(ig),1) = idim(ig)
            ips(ig,nbc(ig),2) = ilo
            ips(ig,nbc(ig),3) = jlo
            ips(ig,nbc(ig),4) = idim(ig)
            ips(ig,nbc(ig),5) = ihi
            ips(ig,nbc(ig),6) = jhi
          else if (if .eq. 5) then
            ips(ig,nbc(ig),1) = jlo
            ips(ig,nbc(ig),2) = 1
            ips(ig,nbc(ig),3) = ilo
            ips(ig,nbc(ig),4) = jhi
            ips(ig,nbc(ig),5) = 1
            ips(ig,nbc(ig),6) = ihi
          else
            ips(ig,nbc(ig),1) = jlo
            ips(ig,nbc(ig),2) = jdim(ig)
            ips(ig,nbc(ig),3) = ilo
            ips(ig,nbc(ig),4) = jhi
            ips(ig,nbc(ig),5) = jdim(ig)
            ips(ig,nbc(ig),6) = ihi
          end if
        end if
      enddo
 1010 continue
c
      n1to1_global=num1to1
      nbocos=n_b
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
      write(7,'(''xxxxxxxxxxxxxx'')')
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
              if (ibctype(n,nn) .eq. 1003) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1001) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1002) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1005) then
                ibctyp=1005
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1006) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1000) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1008) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1011) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1012) then
                ibctyp=1012
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1013) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 2004) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
              if (ibctype(n,nn) .eq. 1003) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1001) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1002) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1005) then
                ibctyp=1005
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1006) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1000) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1008) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1011) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1012) then
                ibctyp=1012
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1013) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,jsta,jend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 2004) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
              if (ibctype(n,nn) .eq. 1003) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1001) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1002) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1005) then
                ibctyp=1005
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1006) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1000) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1008) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1011) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1012) then
                ibctyp=1012
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1013) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 2004) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
              if (ibctype(n,nn) .eq. 1003) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1001) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1002) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1005) then
                ibctyp=1005
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1006) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1000) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1008) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1011) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1012) then
                ibctyp=1012
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 1013) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,ksta,
     +            kend,ndata
              else if (ibctype(n,nn) .eq. 2004) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
              if (ibctype(n,nn) .eq. 1003) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1001) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1002) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1005) then
                ibctyp=1005
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1006) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1000) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1008) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1011) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1012) then
                ibctyp=1012
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1013) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 2004) then
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
          if (i_one_to(m) .eq. n) then
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
          if (i_one_from(m) .eq. n) then
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
              if (ibctype(n,nn) .eq. 1003) then
                ibctyp=1003
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1001) then
                ibctyp=1001
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1002) then
                ibctyp=1002
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1005) then
                ibctyp=1005
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1006) then
                ibctyp=1006
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1000) then
                ibctyp=1000
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1008) then
                ibctyp=1008
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1011) then
                ibctyp=1011
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1012) then
                ibctyp=1012
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 1013) then
                ibctyp=1013
                ndata=0
                write(7,'(8i10)') n,iseg,ibctyp,ista,iend,jsta,
     +            jend,ndata
              else if (ibctype(n,nn) .eq. 2004) then
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
          if (i_one_to(m) .eq. n) igrid=n
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
          if (i_one_from(m) .eq. n) igrid=n
        enddo
        if (igrid .eq. 0) then
          write(6,'('' error... zone not found'')')
          stop
        end if
        if (idonor_range(1,m) .eq. idonor_range(4,m)) then
          if(irange(2,m) .eq. irange(5,m)) then
            if (rev(m) .eq. "TRUE") then
              isva1=2
              isva2=3
            else
              isva1=3
              isva2=2
            end if
          else
            if (rev(m) .eq. "FALSE") then
              isva1=2
              isva2=3
            else
              isva1=3
              isva2=2
            end if
          end if
        else if (idonor_range(2,m) .eq. idonor_range(5,m)) then
          if(irange(2,m) .eq. irange(5,m)) then
            if (rev(m) .eq. "FALSE") then
              isva1=1
              isva2=3
            else
              isva1=3
              isva2=1
            end if
          else
            if (rev(m) .eq. "TRUE") then
              isva1=1
              isva2=3
            else
              isva1=3
              isva2=1
            end if
          end if
        else if (idonor_range(3,m) .eq. idonor_range(6,m)) then
          if(irange(2,m) .eq. irange(5,m)) then
            if (rev(m) .eq. "TRUE") then
              isva1=1
              isva2=2
            else
              isva1=2
              isva2=1
            end if
          else
            if (rev(m) .eq. "FALSE") then
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
c
      return
      end
