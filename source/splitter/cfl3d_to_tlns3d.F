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
      subroutine cfl3d_to_tlns3d(iver,ipatch,iptyp,ipar,nnodes,isd,
     .                           mbloc,msegn,msegt,mxbli,il,jl,kl,
     .                           rkap0g,levelg,igridg,iflimg,ifdsg,
     .                           iviscg,jdimg,kdimg,idimg,idiagg,
     .                           nblcg,idegg,jsg,ksg,jeg,keg,
     .                           ieg,mit,ilamlog,ilamhig,jlamlog,
     .                           jlamhig,klamlog,klamhig,
     .                           iredundant,nseg1,nseg,nsilo,nsihi,
     .                           nsjlo,nsjhi,nsklo,nskhi,nb1,ne1,
     .                           nb2,ne2,ibct,bcval,iovrlp,ifoflg,
     .                           ndat,xmap,imap,ivisb,twotref,itrb1,
     .                           itrb2,jtrb1,jtrb2,ktrb1,ktrb2,iturbb,
     .                           ibif1,ibif2,nbi1,nei1,nbj1,nej1,nbk1,
     .                           nek1,nbi2,nei2,nbj2,nej2,nbk2,nek2,
     .                           nd11,nd21,nd12,nd22,ifsor,nb1s,ne1s,
     .                           nb2s,ne2s,nrotat,ntrans,tlref,rlref,
     .                           ioflag,nou,bou,nbuf,ibufdim,cflout)
c
c     $Id$
c
c***********************************************************************
c
c     Purpose: convert a cfl3d input file (and ronnie input file, if
c              a patched case) to a tlns3d map file 
c
c     initial coding by m.d.sanetrik (6/95)
c
c     modications by r.t.biedron (3/97):
c
c       1) add ability to read zero-range shortcut in cfl3d input file 
c       2) change read order for input bc indices i,k on j-const. faces
c       3) add ability to read non-zero ndata input in cfl3d input file
c       4) add ability to handle a non-patched case
c       5) add option for version 4 or version 5 cfl3d input file
c       6) add description of "type 1 or type 2 patch input file" 
c          during prompt
c       7) added ability to handle a mesh-sequenced case
c       8) changed '(a40)' read of input file names to be '(a60)'
c          to be consistant with cfl3d standard
c       9) removed hardwired ivisb(1-3,mbloc) - now reads ivisc from
c          cfl3d input file to set ivisb(1-3,mbloc)
c      10) changed the order of writing the source and target data
c          when patching to preserve the target data from the ronnie file
c      11) added checks for parameters msegn and mbloc
c      12) set tlns3d bc type variable typ = 0 for an inner cut that 
c          maps to same block, in accordance with the tlns3d standard
c      13) eliminated double counting of 1-1 interfaces 
c      14) for baldwin-lomax cases, set non-zero search ranges
c
c     modications by r.t.biedron (1/98):
c
c       1) changed dimensions of ibif1,ibif2,nbi1...nd22 from msegn
c          to mxbli, and added the necessary parameter statement to
c          set mxbli; also changed dimensions of ifsor...ne2s from
c          msegn to mbloc
c       2) in setting the 1-1 indicies on the source side in the imap
c          array, moved the check for whether imap(8,..,..) should be
c          negatve to a spot ahead of where the check to reverse the
c          source indicies is made; otherwise, the wrong pair of source
c          indicies get swapped if imap(8,..,..) is negative
c
c     NOTES (known conversion limitations):
c
c       1) there seems to be no universal way to go from the cfl3d 
c          "ilamlo, ilamhi..." specification of a laminar region
c          in a field-equation turb. model to the  tlns3d "itrb1,itrb2..."
c          specification of a turbulent region. Thus, this code
c          outputs the "itrb1,itrb2..." as corresponding to the
c          entire index range, and uses the tlns3d convention of
c          zeros for both the starting and ending range to do that 
c       2) when using the type 1 ronnie input, the index range of
c          the source zone(s) is set to the entire range, since any
c          other info is not available from the ronnie input file;
c          this can lead to some discrepancy if converting from a
c          tlns3d map file to a cfl3d/ronnie input file and back 
c          again. this *should* not be a problem, since tlns3d allows
c          for overspecification of the source range, but....
c       3) in tlns3d there is provision for specifying search parameters
c          for the baldwin-lomax turb. model; the cfl3d input file does
c          not have the corresponding information, so if balwin-lomax
c          is used, the the entire segment is assume turbulent, and
c          the search range is set from 3 to max dim. 
c
c**********************************************************************
c
      character*80 gridin,plt3dg,plt3dq,output,resid,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             gridout,sdgridin,sdgridout,dovrlap,dpatch,dresid,
     .             rout
      character*80 gridindum,sdgridindum,cflout
      character*120 bou(ibufdim,nbuf)
c
      dimension nou(nbuf)
      dimension iredundant(msegn,mbloc),nseg1(mbloc)
      dimension nseg(mbloc)
      dimension nsilo(mbloc),nsihi(mbloc),nsjlo(mbloc),
     .          nsjhi(mbloc),nsklo(mbloc),nskhi(mbloc)
      dimension nb1(mbloc,6,msegn),ne1(mbloc,6,msegn),
     .          nb2(mbloc,6,msegn),ne2(mbloc,6,msegn),
     .          ibct(mbloc,6,msegn)
      dimension bcval(mbloc,6,msegn,7)
      dimension iovrlp(mbloc),ifoflg(mbloc,6,msegn),
     .          ndat(mbloc,6,msegn),xmap(msegt,msegn,mbloc)
      dimension imap(msegt,msegn,mbloc)
      dimension ivisb(msegt,mbloc),twotref(msegn,mbloc)
      dimension itrb1(mbloc),itrb2(mbloc),jtrb1(mbloc),jtrb2(mbloc),
     .          ktrb1(mbloc),ktrb2(mbloc),iturbb(mbloc)
c
      dimension ibif1(mxbli),ibif2(mxbli)
      dimension nbi1(mxbli),nei1(mxbli),nbj1(mxbli),
     .          nej1(mxbli),nbk1(mxbli),nek1(mxbli),
     .          nbi2(mxbli),nei2(mxbli),nbj2(mxbli),
     .          nej2(mxbli),nbk2(mxbli),nek2(mxbli),
     .          nd11(mxbli),nd21(mxbli),nd12(mxbli),
     .          nd22(mxbli)
      dimension ifsor(mbloc),nb1s(mbloc),ne1s(mbloc),
     .          nb2s(mbloc),ne2s(mbloc)
      dimension il(mbloc),jl(mbloc),kl(mbloc),rkap0g(mbloc,3),
     .          levelg(mbloc),igridg(mbloc),iflimg(mbloc,3),
     .          ifdsg(mbloc,3),iviscg(mbloc,3),jdimg(mbloc),
     .          kdimg(mbloc),idimg(mbloc),idiagg(mbloc,3),
     .          nblcg(mbloc),idegg(mbloc,3),jsg(mbloc),ksg(mbloc),
     .          jeg(mbloc),keg(mbloc),ieg(mbloc),mit(5,mbloc),
     .          ilamlog(mbloc),ilamhig(mbloc),jlamlog(mbloc),
     .          jlamhig(mbloc),klamlog(mbloc),klamhig(mbloc)
c
      common /unit5/ iunit5
      common /grdinfo/ nbloc
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /info2/ tinf,reue,c2spe,bref,cref,sref,xmc,ymc,zmc,
     .               rfreq,alphau,cloc,cfltau
      common /info3/ ncg,isnd,ialph,irest,iunst,ntstep,
     .               ita,ihstry,nplot3d,nprint,nwrest,ichk,i2d,
     .               mglev(5),nem(5),mitl(5,5),mtt,ndv,isdform,
     .               ip3dgrd,iplt3dtyp
      common /mgv/ epsssc(3),epsssr(3),issc,issr
      common /ron/ titleron(20)
      common /alphait/ cltarg,resupdt,cltol,dalim,dcldal,alphalast,
     .                 ialphit,nttlast,icycupdt
      common /precon1/ cprec,uref,avn,iprecon
      common /cflfiles/gridin,plt3dg,plt3dq,output,resid,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 gridout,sdgridin,sdgridout,dovrlap,dpatch,dresid
      common /ronfiles/rout
c
c*************************************************************************
c--------------------  mapping function description ---------------------
c
c      imap  : mapping function containing topological information
c      msegt : maximum no. of types of operations/boundary conditons
c              (currently set to 50; should be large enough for any case)
c      msegn : maximum no. of segments permitted on all faces of a block
c
c      imap(1 , ) : specifies boundary/operation type
c      imap(2 , ) : specifies face number (1-6)
c      imap(3 , ) : specifies beginning of direction 1
c      imap(4 , ) : specifies end       of direction 1
c      imap(5 , ) : specifies beginning of direction 2
c      imap(6 , ) : specifies end       of direction 2
c
c      imap(8,  ) : if the boundary type is symmetry
c                   takes the value 1,2 or 3 to indicate symmetry 
c                   about x, y or z = constant planes, respectively
c 
c  if the boundary type is a cut/interface
c  additional information described below is required 
c
c      imap(7 , ) : specifies block no. of source segment
c     |imap(8 , )|: specifies face no. of source segment
c      imap(8 , ) > 0 corresponds to direction 1 of source segment
c                      matching with direction 1 of target segment
c                      and same  for direction 2
c      imap(8 , ) < 0 corresponds to direction 2 of source segment
c                      matching with direction 1 of target segment
c                      and vice-versa
c
c      imap(9 , ) : specifies beginning of direction 1 of source segment
c      imap(10, ) : specifies end       of direction 1 of source segment
c      imap(11, ) : specifies beginning of direction 2 of source segment
c      imap(12, ) : specifies end       of direction 2 of source segment
c
c      turbulence/transitional related information
c    
c      imap(13, ) : specifies if there is turbulent flow on this segment
c                   1  corresponds to turbulent flow
c                   0  corresponds to laminar   flow
c      imap(14, ) : begining index in direction 1 for turbulent flow
c      imap(15, ) : ending   index in direction 1 for turbulent flow
c      imap(16, ) : begining index in direction 2 for turbulent flow
c      imap(17, ) : ending   index in direction 2 for turbulent flow
c      imap(18, ) : begining index in normal direction for fmax 
c      imap(19, ) : ending   index in normal direction for fmax 
c      imap(20, ) : ending   index in normal direction for turb. flow 
c
c      the following items are added to facilitate the splitting of 
c      the cfl3d input files, and are not part of the standard tlns3d
c      imap array
c
c      the following additions apply to face segments
c      imap(21, ) : number of additional data for 2000 series bc's
c      imap(22, ) : force flag for this segment
c      imap(23, ) : overlap flag for this block (same for all segments)
c      xmap( 1, ) : 1st additional data entry for 2000 series bc's
c      xmap( 2, ) : 2nd additional data entry for 2000 series bc's
c      xmap( 3, ) : 3rd additional data entry for 2000 series bc's
c      xmap( 4, ) : 4th additional data entry for 2000 series bc's
c      xmap( 5, ) : 5th additional data entry for 2000 series bc's
c      xmap( 6, ) : 6th additional data entry for 2000 series bc's
c      xmap( 7, ) : 7th additional data entry for 2000 series bc's
c      xmap( 8, ) : itrans
c      xmap( 9, ) : rfreqt
c      xmap(10, ) : xmag
c      xmap(11, ) : ymag
c      xmap(12, ) : zmag
c      xmap(13, ) : dxmax
c      xmap(14, ) : dymax
c      xmap(15, ) : dzmax
c      xmap(16, ) : irotat
c      xmap(17, ) : rfreqr
c      xmap(18, ) : thxmag
c      xmap(19, ) : thymag
c      xmap(20, ) : thzmag
c      xmap(21, ) : xorig
c      xmap(22, ) : yorig
c      xmap(23, ) : zorig
c      xmap(24, ) : thxmax
c      xmap(25, ) : thymax
c      xmap(26, ) : thzmax
c
c      the following additions apply to blocks; ivisb(1-3, ) pertain 
c      to the turb model in the i,j,k directions, as before. 
c      ivisb( 4, ) : ncg
c      ivisb( 5, ) : iem
c      ivisb( 6, ) : iadvance
c      ivisb( 7, ) : iforce
c      ivisb( 8, ) : ilamlo
c      ivisb( 9, ) : ilamhi
c      ivisb(10, ) : jlamlo
c      ivisb(11, ) : jlamhi
c      ivisb(12, ) : klamlo
c      ivisb(13, ) : klamhi
c      ivisb(14, ) : inewg
c      ivisb(15, ) : igridc
c      ivisb(16, ) : is
c      ivisb(17, ) : js
c      ivisb(18, ) : ks
c      ivisb(19, ) : ie
c      ivisb(20, ) : je
c      ivisb(21, ) : ke
c      ivisb(22, ) : idiag(i)
c      ivisb(23, ) : idiag(j)
c      ivisb(24, ) : idiag(k)
c      ivisb(25, ) : iflim(i)
c      ivisb(26, ) : iflim(j)
c      ivisb(27, ) : iflim(k)
c      ivisb(28, ) : ifdsg(i)
c      ivisb(29, ) : ifdsg(j)
c      ivisb(30, ) : ifdsg(k)
c      ivisb(31, ) : rkap0g(i) (integer used to represent a real value)
c      ivisb(32, ) : rkap0g(j)                    "
c      ivisb(33, ) : rkap0g(k)                    "
c      ivisb(34, ) : iovrlp
c
c--------------------  boundary/operation type description ---------------------
c
c      nbctype    = imap(1 , )
c                 = 0  corresponds to an inner cut that maps to same block
c                 = 1  corresponds to a cut that maps to another block
c                 = 2  corresponds to a slip (inviscid) wall
c                 = 3  corresponds to a noslip (viscous) wall
c                 = 4  symmetry condition (imap(8) tells about which plane)
c                 = 5  downstream/extrapolation condition
c                 = 6  far-field condition (Riemann invariants)
c
c***************************************************************************
c
c     initial some parameters that need to be set if input version differs
c     from output version
c
      ialph  = 1
      isnd   = 0
      ihstry = 0
      cfltau = 5.0
      c2spe  = 0.0
      nnodes = 1
c
c     zero out xmap and ivisb arrays
c
      do nn=1,msegt
         do ib=1,mbloc
            ivisb(nn,ib) = 0
            do ns=1,msegn
               xmap(nn,ns,ib) = 0
            end do
         end do
      end do
c
      read(10,*)
      read(10,'(a60)') gridindum
      read(10,'(a60)') plt3dg
      read(10,'(a60)') plt3dq
      read(10,'(a60)') output
      read(10,'(a60)') resid
      read(10,'(a60)') turbres
      read(10,'(a60)') blomx
      read(10,'(a60)') output2
      read(10,'(a60)') printout
      read(10,'(a60)') pplunge
      read(10,'(a60)') ovrlap
      read(10,'(a60)') patch
      read(10,'(a60)') restrt
c
c     check for keyword input
c
      iunit5sav = iunit5
      iunit5    = 10
      call readkey(ititr,myid,ibufdim,nbuf,bou,nou,-99,-1)
      iunit5    = iunit5sav
c
      if (ititr.eq.0) then
         read(10,10)(title(i),i=1,20)
      end if
   10 format(20a4)
c
      read(10,*)
      if(iver.eq.4) then 
         if(ipar.eq.0) then
         read(10,*) xmach,alpha,beta,reue,tinf,isnd,c2spe
         else
         read(10,*) xmach,alpha,beta,reue,tinf,isnd,c2spe,nnodes
         end if
      else if(iver.ge.5) then
         read(10,*) xmach,alpha,beta,reue,tinf,ialph,ihstry
      end if
c
c     versions 4.1 and 6/5 have different ways of triggering the
c     Cl-specified option...v4.1 uses the sign on Mach no. as a 
c     trigger, whereas v6 uses the sign on Sref (in v6, the sign
c     on Mach no. is used to trigger the preconditioning option;
c     there is no preconditioning option in v4.1). The location
c     of the extra data for the specification of Cl also differs
c     between the two versions.
c
      ialphit = 0
      iprecon = 0
      if (xmach .lt. 0.) then
         if (iver .eq. 4) then
               ialphit = 1
         else
               iprecon = 1
         end if
         xmach = abs(xmach)
      end if
c
      if (iver .eq. 4) then
         read(10,*)
         if (ialphit .eq. 0) then
            read (10,*) sref,cref,bref,xmc,ymc,zmc
         else
            read (10,*) sref,cref,bref,xmc,ymc,zmc,cltarg,resupdt
            resupdt = -abs(resupdt)
         end if
         read(10,*)
         read(10,*) dt,irest,iflagts,fmax,iunst,rfreq,alphau,cloc
      else
         read(10,*)
         read (10,*) sref,cref,bref,xmc,ymc,zmc
         read(10,*)
         read(10,*) dt,irest,iflagts,fmax,iunst,cfltau
      end if
c
c     read in number of blocs
c
      read(10,*)
      read(10,*) nbloc,nplot3d,nprint,nwrest,ichk,i2d,ntstep,ita
      ip3dgrd = 0
      if (nbloc .lt. 0) ip3dgrd = 1
      nbloc = abs(nbloc)
      if (nbloc .gt. mbloc) then
         write(6,*) 'stopping: parameter mbloc must be ',
     .              'at least ',nbloc
         call termn8(0,-3,ibufdim,nbuf,bou,nou)
      end if
c
      read(10,*)
      do 100 ibloc=1,nbloc
      read(10,*) ncg,iem,iadv,ifor,ivisb(1,ibloc),ivisb(2,ibloc),
     .           ivisb(3,ibloc)
      ivisb( 4,ibloc) = ncg
      ivisb( 5,ibloc) = iem
      ivisb( 6,ibloc) = iadv
      ivisb( 7,ibloc) = ifor
c
      itr = 0
      do iv =1,3
         itr = max(itr,abs(ivisb(iv,ibloc)))
      end do
      iturbb(ibloc) = 0
      if (itr.ge.2) then
         iturbb(ibloc) = 1
      end if
c
c     itrb1, etc for tlns3d only - use for placeholder in map file
c
      itrb1(  ibloc)  = 0
      itrb2(  ibloc)  = 0
      jtrb1(  ibloc)  = 0
      jtrb2(  ibloc)  = 0
      ktrb1(  ibloc)  = 0
      ktrb2(  ibloc)  = 0
  100 continue
c
c     read in size of each block
c
      read(10,*)
      do 110 ibloc=1,nbloc
      read(10,*) il(ibloc),jl(ibloc),kl(ibloc)
  110 continue
c
c
c     ilamlo, etc. NOT PRESERVED during spliiting...issue a warning
c     if any of these laminar regions are active (i.e. non-zero range)
c
      ilamflg = 0
      isum    = 0
c
      read(10,*)
      do 120 ibloc=1,nbloc
      read(10,*) ilamlo,ilamhi,jlamlo,jlamhi,klamlo,klamhi
      ivisb( 8,ibloc) = ilamlo
      ivisb( 9,ibloc) = ilamhi
      ivisb(10,ibloc) = jlamlo
      ivisb(11,ibloc) = jlamhi
      ivisb(12,ibloc) = klamlo
      ivisb(13,ibloc) = klamhi 
      isum = isum + ilamlo+ilamhi+jlamlo+jlamhi+klamlo+klamhi 
  120 continue
      if (isum .gt.0 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: laminar regions specified...these '',
     .             ''are NOT preserved during splitting'')')
         write(6,*)
      end if
c
c     embedded grid data NOT PRESERVED during spliiting...issue a
c     warning if any of the embedded-grid data is non-zero
c
      iembflg = 0
      isum    = 0
c
      read(10,*)
      do 121 ibloc=1,nbloc
      read(10,*) inewg,igridc,is,js,ks,ie,je,ke
      ivisb(14,ibloc) = inewg
      ivisb(15,ibloc) = igridc
      ivisb(16,ibloc) = is
      ivisb(17,ibloc) = js
      ivisb(18,ibloc) = ks
      ivisb(19,ibloc) = ie
      ivisb(20,ibloc) = je
      ivisb(21,ibloc) = ke
      sum = isum + inewg+igridc+is+js+ks+ie+je+ke
  121 continue
      if (isum .gt.0 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: embeded grid data specified...these '',
     .             ''are NOT preserved during splitting'')')
         write(6,*)
      end if
c
      read(10,*)
      do 122 ibloc=1,nbloc
      read(10,*) idi,idj,idk,ifi,ifj,ifk
      ivisb(22,ibloc) = idi
      ivisb(23,ibloc) = idj
      ivisb(24,ibloc) = idk
      ivisb(25,ibloc) = ifi
      ivisb(26,ibloc) = ifj
      ivisb(27,ibloc) = ifk
  122 continue
c
      read(10,*)
      do 123 ibloc=1,nbloc
      read(10,*) ifdsi,ifdsj,ifdsk,rkapi,rkapj,rkapk
      ivisb(28,ibloc) = ifdsi
      ivisb(29,ibloc) = ifdsj
      ivisb(30,ibloc) = ifdsk
c     integer values to represent real input values for rkap
      if(rkapi.lt.-0.9  .and. rkapi.gt.-1.1)  irkapi = -1
      if(rkapi.lt.0.1   .and. rkapi.gt.-0.1)  irkapi =  0
      if(rkapi.lt.1.1   .and. rkapi.gt. 0.9)  irkapi =  1
      if(rkapi.lt.0.40  .and. rkapi.gt. 0.25) irkapi =  3
      if(rkapj.lt.-0.9  .and. rkapj.gt.-1.1)  irkapj = -1
      if(rkapj.lt.0.1   .and. rkapj.gt.-0.1)  irkapj =  0
      if(rkapj.lt.1.1   .and. rkapj.gt. 0.9)  irkapj =  1
      if(rkapj.lt.0.40  .and. rkapj.gt. 0.25) irkapj =  3
      if(rkapk.lt.-0.9  .and. rkapk.gt.-1.1)  irkapk = -1
      if(rkapk.lt.0.1   .and. rkapk.gt.-0.1)  irkapk =  0
      if(rkapk.lt.1.1   .and. rkapk.gt. 0.9)  irkapk =  1
      if(rkapk.lt.0.40  .and. rkapk.gt. 0.25) irkapk =  3
      ivisb(31,ibloc) = irkapi
      ivisb(32,ibloc) = irkapj
      ivisb(33,ibloc) = irkapk
  123 continue
c
c     initialize bcval and ifoflg array
c
      do ibloc=1,mbloc
         do iface = 1,6
            do iseg  = 1,msegn
               do ll =1,7
                  bcval(ibloc,iface,iseg,ll) = 0.
               end do
               ifoflg(ibloc,iface,iseg) = 1
            end do
         end do
      end do
c
c     initialize bc data file counter
c
      ibcfilen = 0
c
c     read in number of segments on each face of each block
c
      read(10,*)
      do 130 ibloc=1,nbloc
      read(10,*) ifoo,nsilo(ibloc),nsihi(ibloc),nsjlo(ibloc),
     &                nsjhi(ibloc),nsklo(ibloc),nskhi(ibloc),
     &                iovrlp(ibloc)
      ivisb(34,ibloc) = iovrlp(ibloc)
  130 continue
c
c     set up non-interface segments for i=imin faces
c
      read(10,*)
      do 140 ibloc=1,nbloc
      do 140 iseg=1,nsilo(ibloc)
      iface = 1
      read(10,*) ib,is,ibct(ibloc,iface,iseg),
     &           nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),ndata
      call shortinp(nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           il(ibloc),jl(ibloc),kl(ibloc),iface)
      if (ibct(ibloc,iface,iseg).eq.2005 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: periodic bc data specified...these '',
     .             ''are NOT preserved during splitting'')')
      else if (ibct(ibloc,iface,iseg).eq.2006 
     .         .and. cflout .ne. 'null') then
         write(6,'(''WARNING: radial equilibrium bc data specified'',
     .             ''...these MAY NOT be preserved during splitting'')')
      end if
      ndat(ibloc,iface,iseg) = ndata
      if (abs(ndata).gt.0) then
         read(10,*)
         if (ndata.gt.0) then
            read(10,*) (bcval(ibloc,iface,iseg,ll),ll=1,ndata)
         else
c           can't actually deal with file name yet, so just assign
c           a number to the file
            read(10,*)
            ibcfilen = ibcfilen + 1
            bcval(ibloc,iface,iseg,1) = ibcfilen
         end if
      end if
      if (is.lt.0) ifoflg(ibloc,iface,iseg) = 0
  140 continue
c
c     set up non-interface segments for i=imax faces
c
      read(10,*)
      do 145 ibloc=1,nbloc
      do 145 iseg=1,nsihi(ibloc)
      iface = 2
      read(10,*) ib,is,ibct(ibloc,iface,iseg),
     &           nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),ndata
      call shortinp(nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           il(ibloc),jl(ibloc),kl(ibloc),iface)
      if (ibct(ibloc,iface,iseg).eq.2005 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: periodic bc data specified...these '',
     .             ''are NOT preserved during splitting'')')
      else if (ibct(ibloc,iface,iseg).eq.2006
     .         .and. cflout .ne. 'null') then
         write(6,'(''WARNING: radial equilibrium bc data specified'',
     .             ''...these MAY NOT be preserved during splitting'')')
      end if
      ndat(ibloc,iface,iseg) = ndata
      if (abs(ndata).gt.0) then
         read(10,*)
         if (ndata.gt.0) then
            read(10,*) (bcval(ibloc,iface,iseg,ll),ll=1,ndata)
         else
c           can't actually deal with file name yet, so just assign
c           a number to the file
            read(10,*)
            ibcfilen = ibcfilen + 1
            bcval(ibloc,iface,iseg,1) = ibcfilen
         end if
      end if
      if (is.lt.0) ifoflg(ibloc,iface,iseg) = 0
  145 continue
c
c     set up non-interface segments for j=jmin faces
c
      read(10,*)
      do 150 ibloc=1,nbloc
      do 150 iseg=1,nsjlo(ibloc)
      iface = 3
      read(10,*) ib,is,ibct(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),ndata
      call shortinp(nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           il(ibloc),jl(ibloc),kl(ibloc),iface)
      if (ibct(ibloc,iface,iseg).eq.2005 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: periodic bc data specified...these '',
     .             ''are NOT preserved during splitting'')')
      else if (ibct(ibloc,iface,iseg).eq.2006
     .         .and. cflout .ne. 'null') then
         write(6,'(''WARNING: radial equilibrium bc data specified'',
     .             ''...these MAY NOT be preserved during splitting'')')
      end if
      ndat(ibloc,iface,iseg) = ndata
      if (abs(ndata).gt.0) then
         read(10,*)
         if (ndata.gt.0) then
            read(10,*) (bcval(ibloc,iface,iseg,ll),ll=1,ndata)
         else
c           can't actually deal with file name yet, so just assign
c           a number to the file
            ibcfilen = ibcfilen + 1
            read(10,*)
            bcval(ibloc,iface,iseg,1) = ibcfilen
         end if
      end if
      if (is.lt.0) ifoflg(ibloc,iface,iseg) = 0
  150 continue
c
c     set up non-interface segments for j=jmax faces
c
      read(10,*)
      do 155 ibloc=1,nbloc
      do 155 iseg=1,nsjhi(ibloc)
      iface = 4
      read(10,*) ib,is,ibct(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),ndata
      call shortinp(nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           il(ibloc),jl(ibloc),kl(ibloc),iface)
      if (ibct(ibloc,iface,iseg).eq.2005 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: periodic bc data specified...these '',
     .             ''are NOT preserved during splitting'')')
      else if (ibct(ibloc,iface,iseg).eq.2006
     .         .and. cflout .ne. 'null') then
         write(6,'(''WARNING: radial equilibrium bc data specified'',
     .             ''...these MAY NOT be preserved during splitting'')')
      end if
      ndat(ibloc,iface,iseg) = ndata
      if (abs(ndata).gt.0) then
         read(10,*)
         if (ndata.gt.0) then
            read(10,*) (bcval(ibloc,iface,iseg,ll),ll=1,ndata)
         else
c           can't actually deal with file name yet, so just assign
c           a number to the file
            read(10,*)
            ibcfilen = ibcfilen + 1
            bcval(ibloc,iface,iseg,1) = ibcfilen
         end if
      end if
      if (is.lt.0) ifoflg(ibloc,iface,iseg) = 0
  155 continue
c
c     set up non-interface segments for k=kmin faces
c
      read(10,*)
      do 160 ibloc=1,nbloc
      do 160 iseg=1,nsklo(ibloc)
      iface = 5
      read(10,*) ib,is,ibct(ibloc,iface,iseg),
     &           nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),ndata
      call shortinp(nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           il(ibloc),jl(ibloc),kl(ibloc),iface)
      if (ibct(ibloc,iface,iseg).eq.2005 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: periodic bc data specified...these '',
     .             ''are NOT preserved during splitting'')')
      else if (ibct(ibloc,iface,iseg).eq.2006
     .         .and. cflout .ne. 'null') then
         write(6,'(''WARNING: radial equilibrium bc data specified'',
     .             ''...these MAY NOT be preserved during splitting'')')
      end if
      ndat(ibloc,iface,iseg) = ndata
      if (abs(ndata).gt.0) then
         read(10,*)
         if (ndata.gt.0) then
            read(10,*) (bcval(ibloc,iface,iseg,ll),ll=1,ndata)
         else
c           can't actually deal with file name yet, so just assign
c           a number to the file
            read(10,*)
            ibcfilen = ibcfilen + 1
            bcval(ibloc,iface,iseg,1) = ibcfilen
         end if
      end if
      if (is.lt.0) ifoflg(ibloc,iface,iseg) = 0
  160 continue
c
c     set up non-interface segments for k=kmax faces
c
      read(10,*)
      do 165 ibloc=1,nbloc
      do 165 iseg=1,nskhi(ibloc)
      iface = 6
      read(10,*) ib,is,ibct(ibloc,iface,iseg),
     &           nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),ndata
      call shortinp(nb1(ibloc,iface,iseg),ne1(ibloc,iface,iseg),
     &           nb2(ibloc,iface,iseg),ne2(ibloc,iface,iseg),
     &           il(ibloc),jl(ibloc),kl(ibloc),iface)
      if (ibct(ibloc,iface,iseg).eq.2005 .and. cflout .ne. 'null') then
         write(6,'(''WARNING: periodic bc data specified...these '',
     .             ''are NOT preserved during splitting'')')
      else if (ibct(ibloc,iface,iseg).eq.2006
     .         .and. cflout .ne. 'null') then
         write(6,'(''WARNING: radial equilibrium bc data specified'',
     .             ''...these MAY NOT be preserved during splitting'')')
      end if
      ndat(ibloc,iface,iseg) = ndata
      if (abs(ndata).gt.0) then
         read(10,*)
         if (ndata.gt.0) then
            read(10,*) (bcval(ibloc,iface,iseg,ll),ll=1,ndata)
         else
c           can't actually deal with file name yet, so just assign
c           a number to the file
            read(10,*)
            ibcfilen = ibcfilen + 1
            bcval(ibloc,iface,iseg,1) = ibcfilen
         end if
      end if
      if (is.lt.0) ifoflg(ibloc,iface,iseg) = 0
  165 continue
c
c     iniitialize imap array
c
      do iii = 1,msegt
         do nnsg = 1,msegn
            do ibl = 1,mbloc
               imap(iii,nnsg,ibl) = 0
            end do
         end do
      end do
c
c     set up imap array using the info that we have so far
c
      do 200 ibl = 1,nbloc
      nnsg = 0
      do 201 ifc = 1,6
      if(ifc .eq. 1) nsf = nsilo(ibl)
      if(ifc .eq. 2) nsf = nsihi(ibl)
      if(ifc .eq. 3) nsf = nsjlo(ibl)
      if(ifc .eq. 4) nsf = nsjhi(ibl)
      if(ifc .eq. 5) nsf = nsklo(ibl)
      if(ifc .eq. 6) nsf = nskhi(ibl)
      do 202 isg=1,nsf
      if(ibct(ibl,ifc,isg) .ne. 0) then
         nnsg      = nnsg +1
         nseg(ibl) = nnsg
         if (nseg(ibl) .gt. msegn) then
            write(6,*) 'stopping: parameter msegn must be ',
     .                 'at least ',nseg(ibl)
            call termn8(0,-7,ibufdim,nbuf,bou,nou)
         end if
c
         imap( 1,nnsg,ibl) = ibct(ibl,ifc,isg)
         imap( 2,nnsg,ibl) = ifc
         imap( 3,nnsg,ibl) = nb1(ibl,ifc,isg)
         imap( 4,nnsg,ibl) = ne1(ibl,ifc,isg)
         imap( 5,nnsg,ibl) = nb2(ibl,ifc,isg)
         imap( 6,nnsg,ibl) = ne2(ibl,ifc,isg)
         imap( 7,nnsg,ibl) = 0
         imap( 8,nnsg,ibl) = 0
         if(imap(1,nnsg,ibl) .eq. 4) then
            idir = 2
            if(ifc .le. 2) idir = 1
            if(ifc .ge. 5) idir = 3
            imap( 8,nnsg,ibl) = idir
         end if
         imap( 9,nnsg,ibl) = 0
         imap(10,nnsg,ibl) = 0
         imap(11,nnsg,ibl) = 0
         imap(12,nnsg,ibl) = 0
c
         imap(21,nnsg,ibl) = ndat(ibl,ifc,isg)
         imap(22,nnsg,ibl) = ifoflg(ibl,ifc,isg)
         imap(23,nnsg,ibl) = iovrlp(ibl)
         xmap( 1,nnsg,ibl) = bcval(ibl,ifc,isg,1)
         xmap( 2,nnsg,ibl) = bcval(ibl,ifc,isg,2)
         xmap( 3,nnsg,ibl) = bcval(ibl,ifc,isg,3)
         xmap( 4,nnsg,ibl) = bcval(ibl,ifc,isg,4)
         xmap( 5,nnsg,ibl) = bcval(ibl,ifc,isg,5)
         xmap( 6,nnsg,ibl) = bcval(ibl,ifc,isg,6)
         xmap( 7,nnsg,ibl) = bcval(ibl,ifc,isg,7)
c
      end if
  202 continue
  201 continue
  200 continue
c
      read(10,*)
      if (iver .eq. 4) then
         read(10,*) mseq,mgflag,iconsf,mtt,ngam
      else
         if (iprecon .eq. 0) then
            read(10,*) mseq,mgflag,iconsf,mtt,ngam
         else
            read(10,*) mseq,mgflag,iconsf,mtt,ngam,cprec,uref,avn
         end if
      end if
      read(10,*)
      read(10,*) issc,(epsssc(i),i=1,3),issr,(epsssr(i),i=1,3)
      read(10,*)
      do m=1,mseq
        read(10,*) ncyc1(m),mglev(m),nem(m),nitfo1(m)
      enddo
      read(10,*)
      do m=1,mseq
        read(10,*) (mitl(i,m),i=1,(mglev(m)+nem(m)))
      enddo
      read(10,*)
      read(10,*)
c
c     set up point-to-point interface segments
c
      read(10,*) nptpif
c
      if (iver.eq.4) then
         read(10,*)
         do 210 iptpif=1,nptpif
         read(10,*) ifoo,ibif1(iptpif),ibif2(iptpif)
  210    continue
c
         read(10,*) 
         do 220 iptpif=1,nptpif
         read(10,*) nbi1(iptpif),nbj1(iptpif),nbk1(iptpif),
     &              nei1(iptpif),nej1(iptpif),nek1(iptpif),
     &              nd11(iptpif),nd21(iptpif)
  220    continue
c
         read(10,*) 
         do 225 iptpif=1,nptpif
         read(10,*) nbi2(iptpif),nbj2(iptpif),nbk2(iptpif),
     &              nei2(iptpif),nej2(iptpif),nek2(iptpif),
     &              nd12(iptpif),nd22(iptpif)
  225    continue
      else
         read(10,*)
         do 222 iptpif=1,nptpif
         read(10,*) ifoo,ibif1(iptpif),
     &              nbi1(iptpif),nbj1(iptpif),nbk1(iptpif),
     &              nei1(iptpif),nej1(iptpif),nek1(iptpif),
     &              nd11(iptpif),nd21(iptpif)
  222    continue
         read(10,*)
         do 223 iptpif=1,nptpif
         read(10,*) ifoo,ibif2(iptpif),
     &              nbi2(iptpif),nbj2(iptpif),nbk2(iptpif),
     &              nei2(iptpif),nej2(iptpif),nek2(iptpif),
     &              nd12(iptpif),nd22(iptpif)
  223    continue
      end if
c
      do 230 iptpif=1,nptpif
c
c     determine face type for target and source interfaces
c
      if(nbi1(iptpif) .eq. nei1(iptpif)) then
         nface1 = 1
         if(nbi1(iptpif) .ne. 1) nface1 = 2
      else if(nbj1(iptpif) .eq. nej1(iptpif)) then
         nface1 = 3
         if(nbj1(iptpif) .ne. 1) nface1 = 4
      else
         nface1 = 5
         if(nbk1(iptpif) .ne. 1) nface1 = 6
      end if
c
      if(nbi2(iptpif) .eq. nei2(iptpif)) then
         nface2 = 1
         if(nbi2(iptpif) .ne. 1) nface2 = 2
      else if(nbj2(iptpif) .eq. nej2(iptpif)) then
         nface2 = 3
         if(nbj2(iptpif) .ne. 1) nface2 = 4
      else
         nface2 = 5
         if(nbk2(iptpif) .ne. 1) nface2 = 6
      end if
c
c     block 1 side (1=target, 2=source)
c
c     set up everything but the index ranges for the segment
c
      ibl       = ibif1(iptpif)
      nseg(ibl) = nseg(ibl) +1
      if (nseg(ibl) .gt. msegn) then
         write(6,*) 'stopping: parameter msegn must be ',
     .              'at least ',nseg(ibl)
         call termn8(0,-7,ibufdim,nbuf,bou,nou)
      end if
      isg       = nseg(ibl)
c
      imap( 1,isg,ibl) = 1
      imap( 2,isg,ibl) = nface1
      imap( 7,isg,ibl) = ibif2(iptpif)
      imap( 8,isg,ibl) = nface2
      imap(22,isg,ibl) = 1
      imap(23,isg,ibl) = iovrlp(ibl)
c
c     determine whether nfaces should be negative
c
      icyc1 = 1
      icyc2 = 1
c
      if(nface1 .le. 2) then
         if(nd11(iptpif) .ne. 2) icyc1 = -1
      else if(nface1 .ge. 5) then
         if(nd11(iptpif) .ne. 1) icyc1 = -1
      else
         if(nd11(iptpif) .ne. 3) icyc1 = -1
      end if
c
      if(nface2 .le. 2) then
         if(nd12(iptpif) .ne. 2) icyc2 = -1
      else if(nface2 .ge. 5) then
         if(nd12(iptpif) .ne. 1) icyc2 = -1
      else
         if(nd12(iptpif) .ne. 3) icyc2 = -1
      end if
c
      if(icyc1*icyc2 .lt. 0) imap(8,isg,ibl) = -imap(8,isg,ibl)
c
      irev1 = 0
      irev2 = 0
c
c     i-face (target)
c
      if(nface1 .le. 2) then
         imap( 3,isg,ibl) = nbj1(iptpif)
         imap( 4,isg,ibl) = nej1(iptpif)
         imap( 5,isg,ibl) = nbk1(iptpif)
         imap( 6,isg,ibl) = nek1(iptpif)
         if(nej1(iptpif) .lt. nbj1(iptpif)) then
            irev1 = 1
            imap( 3,isg,ibl) = nej1(iptpif)
            imap( 4,isg,ibl) = nbj1(iptpif)
         end if
         if(nek1(iptpif) .lt. nbk1(iptpif)) then
            irev2 = 1
            imap( 5,isg,ibl) = nek1(iptpif)
            imap( 6,isg,ibl) = nbk1(iptpif)
         end if
c
c     k-face (target)
c
      else if(nface1 .ge. 5) then
         imap( 3,isg,ibl) = nbi1(iptpif)
         imap( 4,isg,ibl) = nei1(iptpif)
         imap( 5,isg,ibl) = nbj1(iptpif)
         imap( 6,isg,ibl) = nej1(iptpif)
         if(nei1(iptpif) .lt. nbi1(iptpif)) then
            irev1 = 1
            imap( 3,isg,ibl) = nei1(iptpif)
            imap( 4,isg,ibl) = nbi1(iptpif)
         end if
         if(nej1(iptpif) .lt. nbj1(iptpif)) then
            irev2 = 1
            imap( 5,isg,ibl) = nej1(iptpif)
            imap( 6,isg,ibl) = nbj1(iptpif)
         end if
c
c     j-face (target)
c
      else
         imap( 3,isg,ibl) = nbk1(iptpif)
         imap( 4,isg,ibl) = nek1(iptpif)
         imap( 5,isg,ibl) = nbi1(iptpif)
         imap( 6,isg,ibl) = nei1(iptpif)
         if(nek1(iptpif) .lt. nbk1(iptpif)) then
            irev1 = 1
            imap( 3,isg,ibl) = nek1(iptpif)
            imap( 4,isg,ibl) = nbk1(iptpif)
         end if
         if(nei1(iptpif) .lt. nbi1(iptpif)) then
            irev2 = 1
            imap( 5,isg,ibl) = nei1(iptpif)
            imap( 6,isg,ibl) = nbi1(iptpif)
         end if
      end if
c
c     i-face (source)
c
      if(nface2 .le. 2) then
         imap( 9,isg,ibl) = nbj2(iptpif)
         imap(10,isg,ibl) = nej2(iptpif)
         imap(11,isg,ibl) = nbk2(iptpif)
         imap(12,isg,ibl) = nek2(iptpif)
         if(irev1 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap( 9,isg,ibl) = nej2(iptpif)
            imap(10,isg,ibl) = nbj2(iptpif)
            else
            imap(11,isg,ibl) = nek2(iptpif)
            imap(12,isg,ibl) = nbk2(iptpif)
            end if
         end if
         if(irev2 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap(11,isg,ibl) = nek2(iptpif)
            imap(12,isg,ibl) = nbk2(iptpif)
            else
            imap( 9,isg,ibl) = nej2(iptpif)
            imap(10,isg,ibl) = nbj2(iptpif)
            end if
         end if
c
c     k-face (source)
c
      else if(nface2 .ge. 5) then
         imap( 9,isg,ibl) = nbi2(iptpif)
         imap(10,isg,ibl) = nei2(iptpif)
         imap(11,isg,ibl) = nbj2(iptpif)
         imap(12,isg,ibl) = nej2(iptpif)
         if(irev1 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap( 9,isg,ibl) = nei2(iptpif)
            imap(10,isg,ibl) = nbi2(iptpif)
            else
            imap(11,isg,ibl) = nej2(iptpif)
            imap(12,isg,ibl) = nbj2(iptpif)
            end if
         end if
         if(irev2 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap(11,isg,ibl) = nej2(iptpif)
            imap(12,isg,ibl) = nbj2(iptpif)
            else
            imap( 9,isg,ibl) = nei2(iptpif)
            imap(10,isg,ibl) = nbi2(iptpif)
            end if
         end if
c
c     j-face (source)
c
      else
         imap( 9,isg,ibl) = nbk2(iptpif)
         imap(10,isg,ibl) = nek2(iptpif)
         imap(11,isg,ibl) = nbi2(iptpif)
         imap(12,isg,ibl) = nei2(iptpif)
         if(irev1 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap( 9,isg,ibl) = nek2(iptpif)
            imap(10,isg,ibl) = nbk2(iptpif)
            else
            imap(11,isg,ibl) = nei2(iptpif)
            imap(12,isg,ibl) = nbi2(iptpif)
            end if
         end if
         if(irev2 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap(11,isg,ibl) = nei2(iptpif)
            imap(12,isg,ibl) = nbi2(iptpif)
            else
            imap( 9,isg,ibl) = nek2(iptpif)
            imap(10,isg,ibl) = nbk2(iptpif)
            end if
         end if
      end if
c
c     set imap( 1,isg,ibl) = 0 for 1-1 segments that map to themselves
      if (ibif1(iptpif) .eq. ibif2(iptpif)) then
         ift1 = abs(imap( 2,isg,ibl))
         ift2 = abs(imap( 8,isg,ibl))
         idel11 = abs(imap( 3,isg,ibl) - imap( 4,isg,ibl))
         idel12 = abs(imap( 5,isg,ibl) - imap( 6,isg,ibl))
         idel21 = abs(imap( 9,isg,ibl) - imap(10,isg,ibl))
         idel22 = abs(imap(11,isg,ibl) - imap(12,isg,ibl))
         if ((ift1 .eq. ift2) .and. (idel11 .eq. idel21) .and.
     .   (idel12 .eq. idel22)) imap( 1,isg,ibl) = 0
      end if
c
c     block 2 side (1=source, 2=target)
c
c     set up everything but the index ranges for the segment
c
      ibl       = ibif2(iptpif)
      nseg(ibl) = nseg(ibl) +1
      if (nseg(ibl) .gt. msegn) then
         write(6,*) 'stopping: parameter msegn must be ',
     .              'at least ',nseg(ibl)
         call termn8(0,-7,ibufdim,nbuf,bou,nou)
      end if
      isg       = nseg(ibl)
c
      imap( 1,isg,ibl) = 1
      imap( 2,isg,ibl) = nface2
      imap( 7,isg,ibl) = ibif1(iptpif)
      imap( 8,isg,ibl) = nface1
      imap(22,isg,ibl) = 1
      imap(23,isg,ibl) = iovrlp(ibl)
c
c     determine whether nfaces should be negative
c
      icyc1 = 1
      icyc2 = 1
c
      if(nface2 .le. 2) then
         if(nd12(iptpif) .ne. 2) icyc2 = -1
      else if(nface2 .ge. 5) then
         if(nd12(iptpif) .ne. 1) icyc2 = -1
      else
         if(nd12(iptpif) .ne. 3) icyc2 = -1
      end if
c
      if(nface1 .le. 2) then
         if(nd11(iptpif) .ne. 2) icyc1 = -1
      else if(nface1 .ge. 5) then
         if(nd11(iptpif) .ne. 1) icyc1 = -1
      else
         if(nd11(iptpif) .ne. 3) icyc1 = -1
      end if
c
      if(icyc1*icyc2 .lt. 0) imap(8,isg,ibl) = -imap(8,isg,ibl)
c
      irev1 = 0
      irev2 = 0
c
c     i-face (target)
c
      if(nface2 .le. 2) then
         imap( 3,isg,ibl) = nbj2(iptpif)
         imap( 4,isg,ibl) = nej2(iptpif)
         imap( 5,isg,ibl) = nbk2(iptpif)
         imap( 6,isg,ibl) = nek2(iptpif)
         if(nej2(iptpif) .lt. nbj2(iptpif)) then
            irev1 = 1
            imap( 3,isg,ibl) = nej2(iptpif)
            imap( 4,isg,ibl) = nbj2(iptpif)
         end if
         if(nek2(iptpif) .lt. nbk2(iptpif)) then
            irev2 = 1
            imap( 5,isg,ibl) = nek2(iptpif)
            imap( 6,isg,ibl) = nbk2(iptpif)
         end if
c
c     k-face (target)
c
      else if(nface2 .ge. 5) then
         imap( 3,isg,ibl) = nbi2(iptpif)
         imap( 4,isg,ibl) = nei2(iptpif)
         imap( 5,isg,ibl) = nbj2(iptpif)
         imap( 6,isg,ibl) = nej2(iptpif)
         if(nei2(iptpif) .lt. nbi2(iptpif)) then
            irev1 = 1
            imap( 3,isg,ibl) = nei2(iptpif)
            imap( 4,isg,ibl) = nbi2(iptpif)
         end if
         if(nej2(iptpif) .lt. nbj2(iptpif)) then
            irev2 = 1
            imap( 5,isg,ibl) = nej2(iptpif)
            imap( 6,isg,ibl) = nbj2(iptpif)
         end if
c
c     j-face (target)
c
      else
         imap( 3,isg,ibl) = nbk2(iptpif)
         imap( 4,isg,ibl) = nek2(iptpif)
         imap( 5,isg,ibl) = nbi2(iptpif)
         imap( 6,isg,ibl) = nei2(iptpif)
         if(nek2(iptpif) .lt. nbk2(iptpif)) then
            irev1 = 1
            imap( 3,isg,ibl) = nek2(iptpif)
            imap( 4,isg,ibl) = nbk2(iptpif)
         end if
         if(nei2(iptpif) .lt. nbi2(iptpif)) then
            irev2 = 1
            imap( 5,isg,ibl) = nei2(iptpif)
            imap( 6,isg,ibl) = nbi2(iptpif)
         end if
      end if
c
c     i-face (source)
c
      if(nface1 .le. 2) then
         imap( 9,isg,ibl) = nbj1(iptpif)
         imap(10,isg,ibl) = nej1(iptpif)
         imap(11,isg,ibl) = nbk1(iptpif)
         imap(12,isg,ibl) = nek1(iptpif)
         if(irev1 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap( 9,isg,ibl) = nej1(iptpif)
            imap(10,isg,ibl) = nbj1(iptpif)
            else
            imap(11,isg,ibl) = nek1(iptpif)
            imap(12,isg,ibl) = nbk1(iptpif)
            end if
         end if
         if(irev2 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap(11,isg,ibl) = nek1(iptpif)
            imap(12,isg,ibl) = nbk1(iptpif)
            else
            imap( 9,isg,ibl) = nej1(iptpif)
            imap(10,isg,ibl) = nbj1(iptpif)
            end if
         end if
c
c     k-face (source)
c
      else if(nface1 .ge. 5) then
         imap( 9,isg,ibl) = nbi1(iptpif)
         imap(10,isg,ibl) = nei1(iptpif)
         imap(11,isg,ibl) = nbj1(iptpif)
         imap(12,isg,ibl) = nej1(iptpif)
         if(irev1 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap( 9,isg,ibl) = nei1(iptpif)
            imap(10,isg,ibl) = nbi1(iptpif)
            else
            imap(11,isg,ibl) = nej1(iptpif)
            imap(12,isg,ibl) = nbj1(iptpif)
            end if
         end if
         if(irev2 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap(11,isg,ibl) = nej1(iptpif)
            imap(12,isg,ibl) = nbj1(iptpif)
            else
            imap( 9,isg,ibl) = nei1(iptpif)
            imap(10,isg,ibl) = nbi1(iptpif)
            end if
         end if
c
c     j-face (source)
c
      else
         imap( 9,isg,ibl) = nbk1(iptpif)
         imap(10,isg,ibl) = nek1(iptpif)
         imap(11,isg,ibl) = nbi1(iptpif)
         imap(12,isg,ibl) = nei1(iptpif)
         if(irev1 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap( 9,isg,ibl) = nek1(iptpif)
            imap(10,isg,ibl) = nbk1(iptpif)
            else
            imap(11,isg,ibl) = nei1(iptpif)
            imap(12,isg,ibl) = nbi1(iptpif)
            end if
         end if
         if(irev2 .eq. 1) then
            if(imap(8,isg,ibl) .gt. 0) then
            imap(11,isg,ibl) = nei1(iptpif)
            imap(12,isg,ibl) = nbi1(iptpif)
            else
            imap( 9,isg,ibl) = nek1(iptpif)
            imap(10,isg,ibl) = nbk1(iptpif)
            end if
         end if
      end if
c
c     set imap( 1,isg,ibl) = 0 for 1-1 segments that map to themselves
      if (ibif1(iptpif) .eq. ibif2(iptpif)) then
         ift1 = abs(imap( 2,isg,ibl))
         ift2 = abs(imap( 8,isg,ibl))
         idel11 = abs(imap( 3,isg,ibl) - imap( 4,isg,ibl))
         idel12 = abs(imap( 5,isg,ibl) - imap( 6,isg,ibl))
         idel21 = abs(imap( 9,isg,ibl) - imap(10,isg,ibl))
         idel22 = abs(imap(11,isg,ibl) - imap(12,isg,ibl))
         if ((ift1 .eq. ift2) .and. (idel11 .eq. idel21) .and.
     .   (idel12 .eq. idel22)) imap( 1,isg,ibl) = 0
      end if
c
  230 continue
c
c     eliminate double-counted 1-1 interfaces
c
      do ibloc = 1,nbloc
c
c        identify those segments that are redundant, and
c        flag them with iredundant() = 1
c        idebug > 0 for debugging purposes to write out 1-1 
c        data before elimination of iredundant 
c
         idebug  = 0
         do is = 1,msegn
            iredundant(is,ibloc) = 0
         end do
         do is1 = 1,nseg(ibloc)
            im11  = imap( 1,is1,ibloc)
            im12  = imap( 2,is1,ibloc)
            im13  = imap( 3,is1,ibloc)
            im14  = imap( 4,is1,ibloc)
            im15  = imap( 5,is1,ibloc)
            im16  = imap( 6,is1,ibloc)
            im17  = imap( 7,is1,ibloc)
            im18  = imap( 8,is1,ibloc)
            im19  = imap( 9,is1,ibloc)
            im110 = imap(10,is1,ibloc)
            im111 = imap(11,is1,ibloc)
            im112 = imap(12,is1,ibloc)
            do is2 = 1,nseg(ibloc)
               if (iredundant(is1,ibloc) .eq.0) then
                  im21  = imap( 1,is2,ibloc)
                  im22  = imap( 2,is2,ibloc)
                  im23  = imap( 3,is2,ibloc)
                  im24  = imap( 4,is2,ibloc)
                  im25  = imap( 5,is2,ibloc)
                  im26  = imap( 6,is2,ibloc)
                  im27  = imap( 7,is2,ibloc)
                  im28  = imap( 8,is2,ibloc)
                  im29  = imap( 9,is2,ibloc)
                  im210 = imap(10,is2,ibloc)
                  im211 = imap(11,is2,ibloc)
                  im212 = imap(12,is2,ibloc)
                  if (is2 .ne. is1) then
                     if ((im11.eq.im21)   .and. (im12.eq.im22)   .and. 
     .                   (im12.eq.im22)   .and. (im13.eq.im23)   .and.
     .                   (im13.eq.im23)   .and. (im14.eq.im24)   .and.
     .                   (im15.eq.im25)   .and. (im16.eq.im26)   .and.
     .                   (im17.eq.im27)   .and. (im18.eq.im28)   .and.
     .                   (im19.eq.im29)   .and. (im110.eq.im210) .and.
     .                   (im111.eq.im211) .and. (im112.eq.im212) ) then
                          iredundant(is2,ibloc) = 1
                          if (idebug.gt.0) then
                          write(6,*)'block ',ibloc,' seg ',is1,
     .                    ' is redundant to ',is2
                          write(6,*)'iredundant(is1) = ',
     .                    iredundant(is1,ibloc),' redundant(is2) = ',
     .                    iredundant(is2,ibloc)
                          end if
                     end if
                  end if
               end if
            end do
         end do
      end do 
c
c     read remainder of input file
c
      if (iver.eq.4) then
c
c        patch data
         read(10,*)
         read(10,*)
         read(10,*) ninter
c
c        plot3d data
         read(10,*)
         read(10,*)
         if (abs(nplot3d).gt.0) then
            do mm=1,abs(nplot3d)
            read(10,*) idum,iplt3dtyp
            end do
         end if
c
c        movie data
         read(10,*)
         if (abs(nplot3d).gt.0) then
            read(10,*) idum
         end if
c
c        print out data
         read(10,*)
         read(10,*)
         if (abs(nprint).gt.0) then
            do mm=1,abs(nprint)
            read(10,*)
            end do
         end if
c
c        sensitivity data
         read(10,*,end=8989)
         read(10,*,end=8989)
         read(10,*,end=8989) ndv,isdform
         read(10,*)
         read(10,'(a60)') sdgridindum
         read(10,'(a60)') dovrlap
         read(10,'(a60)') dpatch
         read(10,'(a60)') dresid
 8989    continue
c
      else 
c
c        patch data
         read(10,*)
         read(10,*)
         read(10,*) ninter
c
c        plot3d data
         read(10,*)
         read(10,*)
         if (abs(nplot3d).gt.0) then
            nplot3d = abs(nplot3d)
            do mm=1,nplot3d
            read(10,*)
            end do
         end if
c
c        movie data
         read(10,*)
         read(10,*) idum
c
c        print out data
         read(10,*)
         read(10,*)
         if (abs(nprint).gt.0) then
            nprint = abs(nprint)
            do mm=1,nprint
            read(10,*)
            end do
         end if
c
c        control surfaces
c
         read(10,*)
         read(10,*)
         read(10,*) ncs
         read(10,*)
         if (ncs .gt. 0) then
            do mm=1,ncs
               read(10,*)
            end do
         end if
c
c        dynamic mesh data
c
         read(10,*,end=999)
         read(10,*)
         read(10,*) ntrans
         read(10,*)
         if (ntrans.gt.0) read(10,*) tlref
         read(10,*)
         if (ntrans.gt.0) then
            do n=1,ntrans
               read(10,*) igd,itrans,rfreqt,xmag,ymag,zmag
               xmap( 8,1,igd) = itrans
               xmap( 9,1,igd) = rfreqt
               xmap(10,1,igd) = xmag
               xmap(11,1,igd) = ymag
               xmap(12,1,igd) = zmag
            end do
         end if
         read(10,*)
         if (ntrans.gt.0) then
            do n=1,ntrans
               read(10,*) igd,dxmax,dymax,dzmax
               xmap(13,1,igd) = dxmax
               xmap(14,1,igd) = dymax
               xmap(15,1,igd) = dzmax
            end do
         end if
c
         read(10,*)
         read(10,*)
         read(10,*) nrotat
         read(10,*)
         if (nrotat.gt.0) read(10,*) rlref
         read(10,*)
         if (nrotat.gt.0) then
            do n=1,nrotat
               read(10,*) igd,irotat,rfreqr,thxmag,thymag,thzmag,
     .                    xorig,yorig,zorig
               do iseg=1,msegn
                  xmap(16,iseg,igd) = irotat 
                  xmap(17,iseg,igd) = rfreqr
                  xmap(18,iseg,igd) = thxmag
                  xmap(19,iseg,igd) = thymag
                  xmap(20,iseg,igd) = thzmag
                  xmap(21,iseg,igd) = xorig
                  xmap(22,iseg,igd) = yorig
                  xmap(23,iseg,igd) = zorig
               end do
            end do
         end if
         read(10,*)
         if (nrotat.gt.0) then
            do n=1,nrotat
               read(10,*) igd,thxmax,thymax,thzmax
               do iseg=1,msegn
                  xmap(24,iseg,igd) = thxmax
                  xmap(25,iseg,igd) = thymax
                  xmap(26,iseg,igd) = thzmax
               end do
            end do
         end if
c
c        ignore dynamic mesh data for now!
c
  999    continue
c
      end if
c
c     read in patched grid info from ronnie input file, if applicable
c
      if (ipatch .ne. 0) then
c
      read(15,*)
      read(15,*) 
      read(15,'(a60)') rout
      read(15,*)
      if(iptyp .gt.0) then
        read(15,*)
        read(15,*) ioflag,itrace
      else
        ioflag = 0
      end if
      read(15,10) (titleron(i),i=1,20)
      read(15,*)
      read(15,*) 
      read(15,*) 
c
      do 300 ibloc=1,nbloc
      read(15,*) 
  300 continue
c
      read(15,*)
      read(15,*) npatif
      read(15,*)
c
      do 310 intr=1,npatif
      read(15,*)
  310 continue
c
      read(15,*)
      if(ioflag.ne.0) read(15,*)
c
      do 350 ipif=1,npatif
      if (ioflag .eq. 0) then
         read(15,*) iif,iftar,nb1t,ne1t,nb2t,ne2t,nfsor,
     &                            (ifsor(is),is=1,nfsor)
c
c        parse target face value & set up index ranges
c
         iblt   = iftar/100
         iftar  = iftar -iblt*100
         iijkt  = iftar/10
         imnmxt = iftar -iijkt*10
         if(iijkt .eq. 1) then
            ift = 1
            if(imnmxt .eq. 2) ift = 2
         else if(iijkt .eq. 2) then
            ift = 3
            if(imnmxt .eq. 2) ift = 4
         else 
            ift = 5
            if(imnmxt .eq. 2) ift = 6
         end if
c
         if(ift .le. 2) then
c
c        i-face
c
            nbeg1t = nb1t
            nend1t = ne1t
            nbeg2t = nb2t
            nend2t = ne2t
            if(nbeg1t .eq. 0) nbeg1t = 1
            if(nend1t .eq. 0) nend1t = jl(iblt)
            if(nbeg2t .eq. 0) nbeg2t = 1
            if(nend2t .eq. 0) nend2t = kl(iblt)
c
         else if(ift .ge. 5) then
c
c        k-face
c
            nbeg1t = nb2t
            nend1t = ne2t
            nbeg2t = nb1t
            nend2t = ne1t
            if(nbeg1t .eq. 0) nbeg1t = 1
            if(nend1t .eq. 0) nend1t = il(iblt)
            if(nbeg2t .eq. 0) nbeg2t = 1
            if(nend2t .eq. 0) nend2t = jl(iblt)
c
         else
c
c        j-face
c
            nbeg1t = nb1t
            nend1t = ne1t
            nbeg2t = nb2t
            nend2t = ne2t
            if(nbeg1t .eq. 0) nbeg1t = 1
            if(nend1t .eq. 0) nend1t = kl(iblt)
            if(nbeg2t .eq. 0) nbeg2t = 1
            if(nend2t .eq. 0) nend2t = il(iblt)
c
         end if
c
         if(nend1t .lt. nbeg1t) then
            ifoo   = nend1t
            nend1t = nbeg1t
            nbeg1t = ifoo
         end if
         if(nend2t .lt. nbeg2t) then
            ifoo   = nend2t
            nend2t = nbeg2t
            nbeg2t = ifoo
         end if
c
c        parse source face values and set up index ranges
c
         do 320 i=1,nfsor
c
         ibls     = ifsor(i)/100
         ifsor(i) = ifsor(i) -ibls*100
         iijks    = ifsor(i)/10
         imnmxs   = ifsor(i) -iijks*10
         if(iijks .eq. 1) then
            ifs = 1
            if(imnmxs .eq. 2) ifs = 2
         else if(iijks .eq. 2) then
            ifs = 3
            if(imnmxs .eq. 2) ifs = 4
         else 
            ifs = 5
            if(imnmxs .eq. 2) ifs = 6
         end if
c
         if(ifs .le. 2) then
c
c        i-face
c
            nbeg1s = 1
            nend1s = jl(ibls)
            nbeg2s = 1
            nend2s = kl(ibls)
c
         else if(ifs .ge. 5) then
c
c        k-face
c
            nbeg1s = 1
            nend1s = il(ibls)
            nbeg2s = 1
            nend2s = jl(ibls)
c
         else
c
c        j-face
c
            nbeg1s = 1
            nend1s = kl(ibls)
            nbeg2s = 1
            nend2s = il(ibls)
c
         end if
c
         nseg(iblt)      = nseg(iblt) +1
         it              = nseg(iblt)
         if (nseg(iblt) .gt. msegn) then
            write(6,*) 'stopping: parameter msegn must be ',
     .                 'at least ',nseg(iblt)
            call termn8(0,-7,ibufdim,nbuf,bou,nou)
         end if
         imap( 1,it,iblt) = -1
         imap( 2,it,iblt) = ift
         imap( 3,it,iblt) = nbeg1t
         imap( 4,it,iblt) = nend1t
         imap( 5,it,iblt) = nbeg2t
         imap( 6,it,iblt) = nend2t
         imap( 7,it,iblt) = ibls
         imap( 8,it,iblt) = ifs
         imap( 9,it,iblt) = nbeg1s
         imap(10,it,iblt) = nend1s
         imap(11,it,iblt) = nbeg2s
         imap(12,it,iblt) = nend2s
         imap(21,it,iblt) = 0
         imap(22,it,iblt) = 1
         imap(23,it,iblt) = iovrlp(iblt)
c
  320    continue
c
c
      else
         read(15,*) iif,iftar,nb1t,ne1t,nb2t,ne2t,nfsor
         do 330 is=1,nfsor
         read(15,*) ifsor(is),nb1s(is),ne1s(is),nb2s(is),ne2s(is)
  330    continue
c
c        parse target face value & set up index ranges
c
         iblt   = iftar/100
         iftar  = iftar -iblt*100
         iijkt  = iftar/10
         imnmxt = iftar -iijkt*10
         if(iijkt .eq. 1) then
            ift = 1
            if(imnmxt .eq. 2) ift = 2
         else if(iijkt .eq. 2) then
            ift = 3
            if(imnmxt .eq. 2) ift = 4
         else 
            ift = 5
            if(imnmxt .eq. 2) ift = 6
         end if
c
         if(ift .le. 2) then
c
c        i-face
c
            nbeg1t = nb1t
            nend1t = ne1t
            nbeg2t = nb2t
            nend2t = ne2t
            if(nbeg1t .eq. 0) nbeg1t = 1
            if(nend1t .eq. 0) nend1t = jl(iblt)
            if(nbeg2t .eq. 0) nbeg2t = 1
            if(nend2t .eq. 0) nend2t = kl(iblt)
c
         else if(ift .ge. 5) then
c
c        k-face
c
            nbeg1t = nb2t
            nend1t = ne2t
            nbeg2t = nb1t
            nend2t = ne1t
            if(nbeg1t .eq. 0) nbeg1t = 1
            if(nend1t .eq. 0) nend1t = il(iblt)
            if(nbeg2t .eq. 0) nbeg2t = 1
            if(nend2t .eq. 0) nend2t = jl(iblt)
c
         else
c
c        j-face
c
            nbeg1t = nb1t
            nend1t = ne1t
            nbeg2t = nb2t
            nend2t = ne2t
            if(nbeg1t .eq. 0) nbeg1t = 1
            if(nend1t .eq. 0) nend1t = kl(iblt)
            if(nbeg2t .eq. 0) nbeg2t = 1
            if(nend2t .eq. 0) nend2t = il(iblt)
c
         end if
c
         if(nend1t .lt. nbeg1t) then
            ifoo   = nend1t
            nend1t = nbeg1t
            nbeg1t = ifoo
         end if
         if(nend2t .lt. nbeg2t) then
            ifoo   = nend2t
            nend2t = nbeg2t
            nbeg2t = ifoo
         end if
c
c        parse source face values and set up index ranges
c
         do 340 i=1,nfsor
c
         ibls     = ifsor(i)/100
         ifsor(i) = ifsor(i) -ibls*100
         iijks    = ifsor(i)/10
         imnmxs   = ifsor(i) -iijks*10
         if(iijks .eq. 1) then
            ifs = 1
            if(imnmxs .eq. 2) ifs = 2
         else if(iijks .eq. 2) then
            ifs = 3
            if(imnmxs .eq. 2) ifs = 4
         else 
            ifs = 5
            if(imnmxs .eq. 2) ifs = 6
         end if
c
         if(ifs .le. 2) then
c
c        i-face
c
            nbeg1s = nb1s(i)
            nend1s = ne1s(i)
            nbeg2s = nb2s(i)
            nend2s = ne2s(i)
            if(nbeg1s .eq. 0) nbeg1s = 1
            if(nend1s .eq. 0) nend1s = jl(ibls)
            if(nbeg2s .eq. 0) nbeg2s = 1
            if(nend2s .eq. 0) nend2s = kl(ibls)
c
         else if(ifs .ge. 5) then
c
c        k-face
c
            nbeg1s = nb2s(i)
            nend1s = ne2s(i)
            nbeg2s = nb1s(i)
            nend2s = ne1s(i)
            if(nbeg1s .eq. 0) nbeg1s = 1
            if(nend1s .eq. 0) nend1s = il(ibls)
            if(nbeg2s .eq. 0) nbeg2s = 1
            if(nend2s .eq. 0) nend2s = jl(ibls)
c
         else
c
c        j-face
c
            nbeg1s = nb1s(i)
            nend1s = ne1s(i)
            nbeg2s = nb2s(i)
            nend2s = ne2s(i)
            if(nbeg1s .eq. 0) nbeg1s = 1
            if(nbeg1s .eq. 0) nend1s = kl(ibls)
            if(nend2s .eq. 0) nbeg2s = 1
            if(nend2s .eq. 0) nend2s = il(ibls)
c
         end if
c
         if(nend1s .lt. nbeg1s) then
            ifoo   = nend1s
            nend1s = nbeg1s
            nbeg1s = ifoo
         end if
         if(nend2s .lt. nbeg2s) then
            ifoo   = nend2s
            nend2s = nbeg2s
            nbeg2s = ifoo
         end if
c
         nseg(iblt)      = nseg(iblt) +1
         if (nseg(iblt) .gt. msegn) then
            write(6,*) 'stopping: parameter msegn must be ',
     .                 'at least ',nseg(iblt)
            call termn8(0,-7,ibufdim,nbuf,bou,nou)
         end if
         it              = nseg(iblt)
         imap( 1,it,iblt) = -1
         imap( 2,it,iblt) = ift
         imap( 3,it,iblt) = nbeg1t
         imap( 4,it,iblt) = nend1t
         imap( 5,it,iblt) = nbeg2t
         imap( 6,it,iblt) = nend2t
         imap( 7,it,iblt) = ibls
         imap( 8,it,iblt) = ifs
         imap( 9,it,iblt) = nbeg1s
         imap(10,it,iblt) = nend1s
         imap(11,it,iblt) = nbeg2s
         imap(12,it,iblt) = nend2s
         imap(21,it,iblt) = 0
         imap(22,it,iblt) = 1
         imap(23,it,iblt) = iovrlp(iblt)
c
  340    continue
c
      end if
c
  350 continue
c
      end if
c
c     set parameters for baldwin-lomax model: imap(13-20,it,iblt)
c     (these are set non-zero only if the wall is viscous)
c     set these to be the entire range if the block is blomax
c     (ivisb(idir,iblt) is the cfl3d turb. model number)
c
      do iblt = 1,nbloc
         do it = 1,nseg(iblt)
c
            imap(13,it,iblt) = 0
            imap(14,it,iblt) = 0
            imap(15,it,iblt) = 0
            imap(16,it,iblt) = 0
            imap(17,it,iblt) = 0
            imap(18,it,iblt) = 0
            imap(19,it,iblt) = 0
            imap(20,it,iblt) = 0
c
            idir = 2
            if(imap(2,it,iblt) .le. 2) idir = 1
            if(imap(2,it,iblt) .ge. 5) idir = 3
c
            if (imap(1,it,iblt).eq.3 .and. 
     .         (abs(ivisb(idir,iblt)).ge.2 .or.
     .         abs(ivisb(idir,iblt)).le.3)) then
               imap(13,it,iblt) = 1
               imap(14,it,iblt) = imap( 3,it,iblt)
               imap(15,it,iblt) = imap( 4,it,iblt)
               imap(16,it,iblt) = imap( 5,it,iblt)
               imap(17,it,iblt) = imap( 6,it,iblt)
c              ranges for normal direction; search range for fmax taken
c              as entire range, minus the 2 points closest to the wall
c              the ending index for the turb direction is taken so
c              that the entire direction is turbulent
               if (imap( 2,it,iblt) .eq. 1) then
                  imap(18,it,iblt) = 3
                  imap(19,it,iblt) = il(iblt)
                  imap(20,it,iblt) = il(iblt)
               else if (imap( 2,it,iblt) .eq. 2) then
                  imap(18,it,iblt) = il(iblt)-2
                  imap(19,it,iblt) = 1
                  imap(20,it,iblt) = 1
               else if (imap( 2,it,iblt) .eq. 3) then
                  imap(18,it,iblt) = 3
                  imap(19,it,iblt) = jl(iblt)
                  imap(20,it,iblt) = jl(iblt)
               else if (imap( 2,it,iblt) .eq. 4) then
                  imap(18,it,iblt) = jl(iblt)-2
                  imap(19,it,iblt) = 1
                  imap(20,it,iblt) = 1
               else if (imap( 2,it,iblt) .eq. 5) then
                  imap(18,it,iblt) = 3
                  imap(19,it,iblt) = kl(iblt)
                  imap(20,it,iblt) = kl(iblt)
               else if (imap( 2,it,iblt) .eq. 6) then
                  imap(18,it,iblt) = kl(iblt)-2
                  imap(19,it,iblt) = 1
                  imap(20,it,iblt) = 1
               end if
            end if
c
c           default to single tref
c
            twotref(it,iblt) = 0.0
c
         end do
      end do
c
c    check for block faces that have not been set; note that in
c    most cases this would be an error in the cfl3d or ronnie input 
c    files. however, if this translator is being used solely to set
c    up input to paul pao's patcher, then there may be certain cases
c    where not all patch boundaries are set in the ronnie input deck
c    used to set up pao's patcher input. If not all block faces are
c    explicitly set in the tlns3d map file, pao's patcher seems to 
c    run into problems. so, query if the user wants to fill in the 
c    missing faces with extrapolation, which don't affect the patch
c    data. 
c
      do ibloc = 1,nbloc
c        ifchkn .ne. if at least one segment on face n has been set
         ifchk1 = 0
         ifchk2 = 0
         ifchk3 = 0
         ifchk4 = 0
         ifchk5 = 0
         ifchk6 = 0
         do iseg  = 1,nseg(ibloc)
            if (imap(2,iseg,ibloc) .eq. 1) ifchk1 = 1
            if (imap(2,iseg,ibloc) .eq. 2) ifchk2 = 1
            if (imap(2,iseg,ibloc) .eq. 3) ifchk3 = 1
            if (imap(2,iseg,ibloc) .eq. 4) ifchk4 = 1
            if (imap(2,iseg,ibloc) .eq. 5) ifchk5 = 1
            if (imap(2,iseg,ibloc) .eq. 6) ifchk6 = 1
         enddo 
         if (ifchk1 .eq. 0) 
     .      write(6,*)'WARNING: block ',ibloc,
     .                ' face imin, has no bc data'
         if (ifchk2 .eq. 0) 
     .      write(6,*)'WARNING: block ',ibloc,
     .                ' face imax, has no bc data'
         if (ifchk3 .eq. 0) 
     .      write(6,*)'WARNING: block ',ibloc,
     .                ' face jmin, has no bc data'
         if (ifchk4 .eq. 0) 
     .      write(6,*)'WARNING: block ',ibloc,
     .                ' face jmax, has no bc data'
         if (ifchk5 .eq. 0) 
     .      write(6,*)'WARNING: block ',ibloc,
     .                ' face kmin, has no bc data'
         if (ifchk6 .eq. 0) 
     .      write(6,*)'WARNING: block ',ibloc,
     .                ' face kmax, has no bc data'
c
         ifill = 0
         if (ifchk1.eq.0 .or. ifchk2.eq.0 .or. ifchk3.eq.0 .or.
     .       ifchk4.eq.0 .or. ifchk5.eq.0 .or. ifchk6.eq.0) then
             write(6,*)' will fill this missing data',
     .                 ' with extrapolation'
             ifill = 1
         end if
c
c        fill in missing faces with extrapolation
c
         if (ifchk1.eq.0 .and. ifill.eq.1) then
            nseg(ibloc) = nseg(ibloc) + 1
            imap(1,nseg(ibloc),ibloc) = 1002
            imap(2,nseg(ibloc),ibloc) = 1
            imap(3,nseg(ibloc),ibloc) = 1
            imap(4,nseg(ibloc),ibloc) = jl(ibloc)
            imap(5,nseg(ibloc),ibloc) = 1
            imap(6,nseg(ibloc),ibloc) = kl(ibloc)
            do ii = 7,20
              imap(ii,nseg(ibloc),ibloc) = 0
            enddo
         end if
         if (ifchk2.eq.0 .and. ifill.eq.1) then
            nseg(ibloc) = nseg(ibloc) + 1
            imap(1,nseg(ibloc),ibloc) = 1002
            imap(2,nseg(ibloc),ibloc) = 2
            imap(3,nseg(ibloc),ibloc) = 1
            imap(4,nseg(ibloc),ibloc) = jl(ibloc)
            imap(5,nseg(ibloc),ibloc) = 1
            imap(6,nseg(ibloc),ibloc) = kl(ibloc)
            do ii = 7,20
              imap(ii,nseg(ibloc),ibloc) = 0
            enddo
         end if
         if (ifchk3.eq.0 .and. ifill.eq.1) then
            nseg(ibloc) = nseg(ibloc) + 1
            imap(1,nseg(ibloc),ibloc) = 1002
            imap(2,nseg(ibloc),ibloc) = 3
            imap(3,nseg(ibloc),ibloc) = 1
            imap(4,nseg(ibloc),ibloc) = kl(ibloc)
            imap(5,nseg(ibloc),ibloc) = 1
            imap(6,nseg(ibloc),ibloc) = il(ibloc)
            do ii = 7,20
              imap(ii,nseg(ibloc),ibloc) = 0
            enddo
         end if
         if (ifchk4.eq.0 .and. ifill.eq.1) then
            nseg(ibloc) = nseg(ibloc) + 1
            imap(1,nseg(ibloc),ibloc) = 1002
            imap(2,nseg(ibloc),ibloc) = 4
            imap(3,nseg(ibloc),ibloc) = 1
            imap(4,nseg(ibloc),ibloc) = kl(ibloc)
            imap(5,nseg(ibloc),ibloc) = 1
            imap(6,nseg(ibloc),ibloc) = il(ibloc)
            do ii = 7,20
              imap(ii,nseg(ibloc),ibloc) = 0
            enddo
         end if
         if (ifchk5.eq.0 .and. ifill.eq.1) then
            nseg(ibloc) = nseg(ibloc) + 1
            imap(1,nseg(ibloc),ibloc) = 1002
            imap(2,nseg(ibloc),ibloc) = 5
            imap(3,nseg(ibloc),ibloc) = 1
            imap(4,nseg(ibloc),ibloc) = il(ibloc)
            imap(5,nseg(ibloc),ibloc) = 1
            imap(6,nseg(ibloc),ibloc) = jl(ibloc)
            do ii = 7,20
              imap(ii,nseg(ibloc),ibloc) = 0
            enddo
         end if
         if (ifchk6.eq.0 .and. ifill.eq.1) then
            nseg(ibloc) = nseg(ibloc) + 1
            imap(1,nseg(ibloc),ibloc) = 1002
            imap(2,nseg(ibloc),ibloc) = 6
            imap(3,nseg(ibloc),ibloc) = 1
            imap(4,nseg(ibloc),ibloc) = il(ibloc)
            imap(5,nseg(ibloc),ibloc) = 1
            imap(6,nseg(ibloc),ibloc) = jl(ibloc)
            do ii = 7,20
              imap(ii,nseg(ibloc),ibloc) = 0
            enddo
         end if
      enddo
c
c     determine nseg1(ibloc) that doesn't include redunant 1-1
c     segments
c
      do ibloc = 1,nbloc
         nseg1(ibloc) = nseg(ibloc)
         do is = 1,nseg(ibloc)
            if (iredundant(is,ibloc).eq.1) then
               nseg1(ibloc) = nseg1(ibloc) - 1
            end if
         end do
      end do
c
c     write out map file
c
      write(20,'("  nbloc")')
      write(20,'(i5)') nbloc
      write(20,'(" ")')
      do 600 ibloc = 1,nbloc
c
      write(20,'("  nseg  ivisb1 ivisb2 ivisb3",
     &           " iturbb itrb1 itrb2 jtrb1 jtrb2 ktrb1 ktrb2")')
      write(20,'(2i7,10i6)')  nseg1(ibloc),(ivisb(n,ibloc),n=1,3)
     &        , iturbb(ibloc),itrb1(ibloc),itrb2(ibloc),jtrb1(ibloc)
     &        ,               jtrb2(ibloc),ktrb1(ibloc),ktrb2(ibloc) 
      write(20,'("begin non-standard map data")')
      write(20,'(12i6)') (ivisb(n,ibloc),n=4,15)
      write(20,'(12i6)') (ivisb(n,ibloc),n=16,27)
      write(20,'(12i6)') (ivisb(n,ibloc),n=28,34)
      write(20,'("end non-standard map data")')
      write(20,'("  nbt  nst   typ  nft 1tmn 1tmx 2tmn 2tmx  ",
     &                       "nbs  nfs 1smn 1smx 2smn 2smx")')
      write(20,'(16x,"ntrb 1beg 1end 2beg 2end fmx1 fmx2 nvis twotrf")')
      write(20,'(" ")')
c
c
c----------  write out imap array ----------------------------------------
c
      iiseg = 0
      do 600 iseg  = 1,nseg(ibloc)
c
      if (iredundant(iseg,ibloc) .eq. 0) then
         iiseg = iiseg + 1
         write(20,'(2i5,i6,11i5)') ibloc,iiseg,
     .              (imap(n,iseg,ibloc),n=1,12)
         write(20,'(15x,8i5,f10.4)') (imap(n,iseg,ibloc),n=13,20)
     .                          ,twotref(iseg,ibloc)
         write(20,'(19x,"0.0   0.0  0.0  1.0  0.0")')
         write(20,'(15x,"begin non-standard map data")')
         write(20,'(15x,5f10.4)') (xmap(n,iseg,ibloc),n=1,5)
         write(20,'(15x,5f10.4)') (xmap(n,iseg,ibloc),n=6,10)
         write(20,'(15x,5f10.4)') (xmap(n,iseg,ibloc),n=11,15)
         write(20,'(15x,5f10.4)') (xmap(n,iseg,ibloc),n=16,20)
         write(20,'(15x,5f10.4)') (xmap(n,iseg,ibloc),n=21,25) 
         write(20,'(15x,5f10.4)') (xmap(n,iseg,ibloc),n=26,26)
         write(20,'(15x,8i5)') (imap(n,iseg,ibloc),n=21,23)
         write(20,'(15x,"end non-standard map data")')
         write(20,'(" ")')
      end if
c
  600 continue
c
c     rewind unsplit input file so it can be read again in order
c     to set keyword data for the split input file.
c
      rewind(10)
c
      return
      end
