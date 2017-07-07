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
#   ifdef FASTIO
      subroutine reass(maxgr,irdrea)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Reassembles the CFL3D output, if "fast I/O" option for
c     writ_buf used
c     Written by Bob Bergeron, NASA Ames, 2005
c***********************************************************************
c
cRJB...this code reassembles the cfl3d output 
cRJB...from files (fort.8000s) and the cfl3d.out
cRJB...written by the MLP version with to the following ppflags
cRJB...could also be used for the MPI version
cRJB...01 mgbl-memory
cRJB...02 mgblk-bc_info
cRJB...03 mgblk-bc connections and conditions
cRJB...04 mgblk-cfl number p/o & flux calc & time step distribution
cRJB...05 mgblk-ivmx-blomx printout
cRJB...06 mgblk-cfl number p/o
cRJB...07 mgblk
cRJB...08 mgblk
cRJB...09 mgblk
cRJB...10 mgblk
cRJB...11 mgblk-interpolating corrections
cRJB...13 qinter
cRJB...14 qinter
cRJB...15 resetg
cRJB...16 setup-grid/block dimensions
cRJB...17 setup-metrics
cRJB...18 setup-cell volumes
cRJB...19 setup-blank array
cRJB...20 setup-initial conditions
cRJB...21 setup
cRJB...22 setup
cRJB...23 setup-restart
cRJB...24 setup
cRJB...25 setup
cRJB...26 setup-minimum distance calc
cRJB...27 setup-turbulent initial conditions
cRJB...28 setup
cRJB...29 updatedg
cRJB...30 updatedg
cRJB...31 updateg
cRJB...32 dynptch
cRJB...33 dynptch
c...
cRJB...the code reads and writes the (short) cfl3d.out until 
cRJB...it finds a "ppflag" string and values written by the master onto cfl3d.out
cRJB...then it reads the fort.8000s written by the workers and upon appropriate matching
cRJB...of ppflag values, the code writes data for the value of ppflag(termed igo)
cRJB...from the fort.8000s onto a new cfl3d.out
cRJB...the fort.8000s contain a 'ppflag-value-nous' record where
cRJB...nous, the number of output strings, allows easy reading/writing of worker records
cRJB...generally,
cRJB...  reass writes the master's cfl3d data onto the new output, and
cRJB...process writes the worker's cfl3d.data onto the new output.
cRJB...this approach is straightforward until we encounter the special cases...
cRJB...code comments should document such cases
cRJB...allowing writ_buf to write ppflags of zero length records assists in the understanding of the IO
cRJB...I try to make the code clear because while I have run some problems,
cRJB...I have not run all problems, and I have not considered all special cases
cRJB...so the user may have to modify the reass code.
cRJB...or request that I do so
c
cRJB...successful completion of reass produces a cfl3d.out which is essentially
cRJB...the same as the corresponding MPI cfl3d output.
      parameter(nproc=2000) ! current hardwired maximum
      character*132 outstr
      character*8  num
      character*3 cigo
      character*80 grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .             output2,printout,pplunge,ovrlap,patch,restrt,
     .             subres,subtur,grdmov,alphahist,errfile,preout,
     .             aeinp,aeout,sdhist,avgg,avgq
      integer nm(nproc)
      integer irdrea(maxgr)
      common/rjbdbgi/lunfio0
      common /rjb/idbug,ierror,ngridd,ngridt,ncard(nproc),nlunitc,
     .  lunitc(nproc),nr,nww,nrec(nproc),igo23,igo27,nblockp,nstopr
      common /filenam/ grid,plt3dg,plt3dq,output,residual,turbres,blomx,
     .                 output2,printout,pplunge,ovrlap,patch,restrt,
     .                 subres,subtur,grdmov,alphahist,errfile,preout,
     .                 aeinp,aeout,sdhist,avgg,avgq
c
      ierror=0
      igo23=0
      igo27=0
      igop=-1
      nblockp=0
      nlunitc=0
      do m=1,maxgr
       ncard(m)=irdrea(m)+lunfio0
      enddo
      do n=1,nproc
       nm(n)=0
       nrec(n)=0
       lunitc(n)=0
      enddo

      ngridt=maxgr
c
      match=0
      do m=1,ngridt
       do n=1,ngridt
        if(nm(n).le.0)then
         if(n.ne.m)then
          if(ncard(n).eq.ncard(m))then
           match=match+1
           nm(n)=nm(n)+1
           nm(m)=nm(m)+1
          endif
         endif
        endif
       enddo
      enddo
      ngridd=ngridt-match
c
      open(lunfio0,file='cfl3d.out.reass',form='formatted',iostat=ios)
      open(10,file=output,form='formatted',iostat=ios)
c     write(lunfio0,'(" reassembling for ",i6," distinct grids ")')ngridd
      if(ngridd.le.0)stop 'ngridd.le.0'
      nr=0
      nw=0
    2 continue
      read(10,'(a)',err=998,end=999)outstr
      nr=nr+1
      j=index(outstr,"ppflag")
      if(j.gt.0 .and. nlunitc.lt.ngridd)then

cRJB..obtain the ppflag flag which will
cRJB..dictate how to read/write/insert the worker output 
        cigo(1:1)=outstr(9:9)
        cigo(2:2)=outstr(10:10)
        cigo(3:3)=outstr(11:11)
        read(cigo,*) igo
cRJB...the MLP shortened cfl3d.out contains duplicate igos and 
cRJB...at least one duplicate igo sequence: 4 8 9
cRJB...written by the master in writ_buf
cRJB...we process the worker files at the end of these sequences
cRJB... prolongation produces a series of 11's so igo=igop can allow process
        if(igo.eq.igop)then
         if(igo.eq.11)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
         elseif(igo.eq.23)then
          go to 2
         else ! igo.eq.11
cRJB...we have seen ppflag 26 written 3 times so we really have to call process
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         endif ! igo.eq.11
        else ! igo.eq.igop
         if(igop.eq.-1)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.1)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.2)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.3)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.4)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
cRJB...021005 for ivmx=2 ppflag 5 appears on 10-lets try first with 0 output on 14
         elseif(igop.eq.5)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.6)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.11)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.15)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.16)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.17)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.18)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.19)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.20)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.21)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.22)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.23)then
cRJB...when we see 23, we have already written first block info 
cRJB...            or  we have already written  last block info 
            backspace (10)
            nr=nr-1
            if(igo23.le.0)then
             call process(igo)
             if(ierror.gt.0)return
            igop=igo
            endif
         elseif(igop.eq.24)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.25)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.26)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.27)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.28)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.29)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.30)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.31)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         elseif(igop.eq.32)then
          backspace (10)
          nr=nr-1
          call process(igo)
          if(ierror.gt.0)return
          igop=igo
         else
          backspace (10)
          nr=nr-1
          igop=igo
          go to 2
         endif ! igop.eq.-1
         if(igo.eq.23)then
            igo23=igo23+1
         endif
        endif ! igo.eq.igop
        igop=igo
        go to 2
      else ! j.gt.0 .and...
cRJB...remove any reass trailing edge blanks
cRJB...ppflag   6 was written by 5-block b38b
       if(j.le.0)then
        write(lunfio0,'(a)')trim(outstr)
        nw=nw+1
       endif
       go to 2
      endif ! j.gt.0
  998 continue
      write(lunfio0,'(" reass: error at nr =",i6)')nr
      write(lunfio0,'(" reass: complete at nw =",i6)')nw
      write(lunfio0,'(" reass-error ",i6)')
      return
  999 continue
      write(lunfio0,'(" reass: complete at nr =",i6)')nr
      write(lunfio0,'(" reass: complete at nw =",i6)')nw
      write(lunfio0,'(" reass-complete ",i6)')
      return
      end
C
      subroutine process(igo)
      parameter(nproc=2000) ! current hardwired maximum
      character*132 outstr
      character*132 string
      character*132 string1
      character*132 string2
      character*3 cigop
      character*3 cnblock
      character*6 cnous
      character*3 cigop45
      character*6 cnous45
      character*3 cigop95
      character*6 cnous95
      common/rjbdbgi/lunfio0
      common /rjb/idbug,ierror,ngridd,ngridt,ncard(nproc),nlunitc,
     .  lunitc(nproc),nr,nww,nrec(nproc),igo23,igo27,nblockp,nstopr
c
      igop=0
      igop95=0
      l27done=0
      if(igo.eq.1)then
cRJB...MEMORY-lets read ONE line
cRJB...but 1 lun may contain memory for several grids
cRJB...the number of distinct units tells us
cRJB...the order is sequential
        do lun=1,ngridd
         lunit=lunfio0+lun
         lunitm=lunit-lunfio0
         ios=0
         read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
         nrec(lunitm)=nrec(lunitm)+1
cRJB...should be ppflag
         j1=index(outstr,"ppflag")
         if(j1.gt.0)then
cRJB..confirm the ppflag flag which will
          cigop(1:1)=outstr(9:9)
          cigop(2:2)=outstr(10:10)
          cigop(3:3)=outstr(11:11)
          read(cigop,*) igop
          cnous(1:1)=outstr(12:12)
          cnous(2:2)=outstr(13:13)
          cnous(3:3)=outstr(14:14)
          cnous(4:4)=outstr(15:15)
          cnous(5:5)=outstr(16:16)
          cnous(6:6)=outstr(17:17)
          read(cnous,*) nous
          if(igo.eq.igop)then
           read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
           nrec(lunitm)=nrec(lunitm)+1
           write(lunfio0,'(a)')trim(outstr)
           nw=nw+1
cRJB...we have found a ppflag on lunitm so move fort.10
           read(10,'(a)',err=998,end=997)outstr
           nr=nr+1
          else
cRJB...this error is serious so we'll let it stand
           write(lunfio0,'(" MATCH ERROR process: value igo has ",
     .      "found igop on lunit ",9i6)')igo,igop,lunit,nrec(lunitm)
           return
          endif
         endif
        enddo
        return
      else
cRJB...probably should be igo=2
cRJB...GRID and others-use cardinal order
  100   continue
        if(igo.eq.23.and.igo23.le.0)then
         lun1=2
         lun2=ngridt
        else
         lun1=1
         lun2=ngridt
        endif
cRJB...we write a ppflag  number to commence the worker text and
cRJB...we write a ppflag  0      to complete the worker text,
cRJB...so main loop will need two passes to process a given ppflag
        do n=lun1,lun2
         lunit=ncard(n)
         lunitm=lunit-lunfio0
         igop=0
         jres=0
         ngood=0
         nrec0=nrec(lunitm)
         if(lunitc(lunitm).le.0)then
         read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
         nrec(lunitm)=nrec(lunitm)+1
         nrec0=-1
cRJB...process does considerable backspacing, but experience indicates that 
c nrec(lunitm) should contain the ppflag
         j2=index(outstr,"ppflag")
         if(j2.gt.0)then
cRJB...we have found a ppflag on lunitm and so move 
cRJB...to the next ppflag on fort.10 and we expect the string to
cRJB...be a ppflag-but suppose not(mgt reading restart)
cRJB...lets test 
          read(10,'(a)',err=998,end=997)string
          j2=index(string,"ppflag")
          if(j2.gt.0)then
           nr=nr+1
          else
cRJB...might be goodstuff
           ngood=1
           nr=nr+1
           write(lunfio0,'(a)')trim(string)
          endif
cRJB..obtain number of output strings(nous) for this ppflag
          cigop(1:1)=outstr(9:9)
          cigop(2:2)=outstr(10:10)
          cigop(3:3)=outstr(11:11)
          read(cigop,*) igop
          cnous(1:1)=outstr(12:12)
          cnous(2:2)=outstr(13:13)
          cnous(3:3)=outstr(14:14)
          cnous(4:4)=outstr(15:15)
          cnous(5:5)=outstr(16:16)
          cnous(6:6)=outstr(17:17)
          read(cnous,*) nous
cRJB...if the igo (value of ppflag) from reass matches the igo from process...
cRJB...lets just read nous lines and make decisions after this loop
          if(igo.eq.igop)then
          if(nous.gt.0)then
           do no=1,nous
            read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
            nrec(lunitm)=nrec(lunitm)+1
cRJB...ppflag of 27 prints only sequential output
            if(igo.eq.27)then
cRJB...nblockp maintains sequential p/o for ppflag of 27
            jset=index(outstr,"setting")
            jtur=index(outstr,"turbulent")
            if(jset.gt.0 .and. jtur.gt.0)then
             cnblock(1:1)=outstr(46:46)
             cnblock(2:2)=outstr(47:47)
             cnblock(3:3)=outstr(48:48)
             read(cnblock,*) nblock
             if(nblock.eq.nblockp+1)then
               write(lunfio0,'(a)')trim(outstr)
               nw=nw+1
               nblockp=nblock
             else
cRJB...backspace the "ppflag  27" and the "setting" and we are done with this lunit access
               l27done=1
               backspace(lunit)
               backspace(lunit)
               nrec(lunitm)=nrec(lunitm)-2
               backspace(10)
               nr=nr-1
             endif
            endif
            else
             write(lunfio0,'(a)')trim(outstr)
             nw=nw+1
            endif ! igo.eq.27
           enddo
          endif ! if(nous.gt.0)
c
cRJB...we do not go to the next lunit because some igo and igop combos requires several 
c      sequential reads/writes off the same lunit to duplicate the output
cRJB...we are now going to test for the special cases and if the igo and igop combo is 
c      NOT one of these special combos, then we should backspace the current lunit as 
c      soon as we know and read the next lunit to preserve the loop counter
c
cRJB...COMMENCE THE TEST FOR THE SPECIAL CASES
          if(l27done.le.0)then
          read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
          nrec(lunitm)=nrec(lunitm)+1
          j3=index(outstr,"ppflag")
cRJB...when we finish reading nous records from lunit, we expect to find a ppflag
          if(j3.gt.0)then
cRJB...we have found a ppflag on lunitm and so move
cRJB...to the next ppflag on fort.10 and we expect the string to be a ppflag
cRJB...if its not(the case of the mgt reading restart) we should write the line to the output
cRJB...lets test
           read(10,'(a)',err=998,end=997)string
           j3=index(string,"ppflag")
           if(j3.gt.0)then
            nr=nr+1
           else
cRJB...the fort.10 string is output-but we should not backspace fort.10 if we write out
cRJB...except on bb1-092205
cRJB...one diff allow the correct bb1 output-
cRJB...skip(PUNT) the write to fort.10 and backspace both fort.10 and lunit to avoid 
cRJB...duplicating the fort.10 in reass
cRJB...(ngood is already 0 to allow backspace)
cRJB...mgt reading restart
             if(igo.eq.23 .and. igop.eq.23)then
              ngood=1
              nr=nr+1
              write(lunfio0,'(a)')trim(string)
cRJB...bb1
             else
              nr=nr+1
             endif
           endif
           cigop(1:1)=outstr(9:9)
           cigop(2:2)=outstr(10:10)
           cigop(3:3)=outstr(11:11)
           read(cigop,*) igop
           cnous(1:1)=outstr(12:12)
           cnous(2:2)=outstr(13:13)
           cnous(3:3)=outstr(14:14)
           cnous(4:4)=outstr(15:15)
           cnous(5:5)=outstr(16:16)
           cnous(6:6)=outstr(17:17)
           read(cnous,*) nous
           if(nous.gt.0)then
cRJB...decisions,decisions
cRJB...testing below occurs after reading the first ppflag for this lunit
cRJB...do we go to the next lunit or do we continue with this lunit?
cRJB...j3>0 => we have just found a ppflag, so we should backspace and complete reading on lunit  
cRJB...some special ppflags and cases make us pause before we do so...
cRJB...igop.eq.8(&nous > 0) -> "adding residual correction" printout follows igop.eq.4
cRJB...igop.eq.9 -> "AF" printout follows igop.eq.4
cRJB...igop.eq.11-> "interpolating correction" 3-line printout follows igop.eq.4
cRJB...igop.eq.27-> "turbulent initial conditions" for all blocks written

              if(igop.eq.8)then
                do no=1,nous
                 read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                 nrec(lunitm)=nrec(lunitm)+1
                 write(lunfio0,'(a)')trim(outstr)
                 nw=nw+1
                enddo
cRJB...and now read the ppflag=9
                read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                nrec(lunitm)=nrec(lunitm)+1
                j39=index(outstr,"ppflag")
cRJB...when we finish reading nous records from lunit, we expect to find a ppflag
                if(j39.gt.0)then
cRJB...we have found a ppflag on lunitm so we will move fort.10
                 cigop(1:1)=outstr(9:9)
                 cigop(2:2)=outstr(10:10)
                 cigop(3:3)=outstr(11:11)
                 read(cigop,*) igop
                 cnous(1:1)=outstr(12:12)
                 cnous(2:2)=outstr(13:13)
                 cnous(3:3)=outstr(14:14)
                 cnous(4:4)=outstr(15:15)
                 cnous(5:5)=outstr(16:16)
                 cnous(6:6)=outstr(17:17)
                 read(cnous,*) nous
                 read(10,'(a)',err=998,end=997)string
                 nr=nr+1
                 if(nous.gt.0)then
                  do no=1,nous
                   read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                   nrec(lunitm)=nrec(lunitm)+1
                   write(lunfio0,'(a)')trim(outstr)
                   nw=nw+1
                  enddo
                 else
                  write(lunfio0,'(" nous39mismatch-stop ",9i6)')igo,igop
     .             ,lunit,nrec(lunitm),10,nr
                  write(6,'(" process:nous39mismatch-stop ")')
                  ierror=1
                  return
                 endif
                else
cRJB...a ppflag mismatch backspaces if we have output data to write
                 write(lunfio0,'(" ppflag39mismatch-stop ",9i6)')igo,
     .            igop,lunit,nrec(lunitm),10,nr
                 write(6,'(" process:ppflag39mismatch-stop ")')
                 ierror=1
                 return
                endif

               elseif(igop.eq.9)then
               do no=1,nous
                read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                nrec(lunitm)=nrec(lunitm)+1
                write(lunfio0,'(a)')trim(outstr)
                nw=nw+1
               enddo

              elseif(igop.eq.11)then
cRJB...2 11's in a row on the same file (from prolongation loop)
cRJB...require separate passes thru the 8000s
cRJB...so we punt
                backspace(lunit)
                nrec(lunitm)=nrec(lunitm)-1
                backspace(10)
                nr=nr-1

c             elseif(igop.eq.27)then
              elseif(igop.eq.27 .and. igo.eq.27)then
cRJB...all turbulent initial condtions for this lunit written now
cRJB...we have read a ppflag of 27 and nous > 0
cRJB...if this unit begins the loop(n=1),
cRJB...read until we exhaust the (ppflag of 27) sequential blocks for this lunit
   27           continue
                if(n.eq.1 .or. igo27.gt.0)then
                if(n.eq.1)igo27=1
cRJB...nous should be 1
                 do no=1,nous
                  read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                  nrec(lunitm)=nrec(lunitm)+1
                  cnblock(1:1)=outstr(46:46)
                  cnblock(2:2)=outstr(47:47)
                  cnblock(3:3)=outstr(48:48)
                  read(cnblock,*) nblock
cRJB...decisions
c                 if(n.eq.1 .or. nblock.eq.nblockp+1)then
                  if(nblock.eq.nblockp+1)then
                   write(lunfio0,'(a)')trim(outstr)
                   nw=nw+1
                   nblockp=nblock
                  else
cRJB...backspace the "ppflag  27" and the "setting" and we are done with this lunit access
                   l27done=1
                   backspace(lunit)
                   backspace(lunit)
                   nrec(lunitm)=nrec(lunitm)-2
                   backspace(10)
                   nr=nr-1
                  endif
                 enddo
                 if(l27done.le.0)then
                 read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                 nrec(lunitm)=nrec(lunitm)+1
                 j27=index(outstr,"ppflag")
                 if(j27.gt.0)then
cRJB...we have found a ppflag on lunitm so we will move fort.10
                  cigop(1:1)=outstr(9:9)
                  cigop(2:2)=outstr(10:10)
                  cigop(3:3)=outstr(11:11)
                  read(cigop,*) igop
                  cnous(1:1)=outstr(12:12)
                  cnous(2:2)=outstr(13:13)
                  cnous(3:3)=outstr(14:14)
                  cnous(4:4)=outstr(15:15)
                  cnous(5:5)=outstr(16:16)
                  cnous(6:6)=outstr(17:17)
                  read(cnous,*) nous
                  read(10,'(a)',err=998,end=997)string
                  nr=nr+1
                  if(igop.eq.27)then
                   goto 27
                  else
                   backspace(lunit)
                   nrec(lunitm)=nrec(lunitm)-1
                   backspace(10)
                   nr=nr-1
                  endif
                else
cRJB...a ppflag mismatch backspaces if we have output data to write
                 write(lunfio0,'(" ppflag27mismatch-stop ",9i6)')igo,
     .            igop,lunit,nrec(lunitm),10,nr
                 write(6,'(" process:ppflag27mismatch-stop ")')
                 ierror=1
                 return
                endif ! j27.gt.0
                endif ! l27done.le.0
cRJB...COMPLETE THE TEST FOR THE SPECIAL CASES
                l27done=0
                else ! n.eq.1
                 backspace(lunit)
                 nrec(lunitm)=nrec(lunitm)-1
                 backspace(10)
                 nr=nr-1
                endif ! n.eq.1

              else  !if(igop.eq.8)
                if(igo.eq.23)then
                backspace(lunit)
                nrec(lunitm)=nrec(lunitm)-1
                else
                backspace(lunit)
                nrec(lunitm)=nrec(lunitm)-1
                backspace(10)
                nr=nr-1
                endif
              endif ! if(igop.eq.8)

            else ! if(nous.gt.0)then
cRJB...nous le 0 on lunit-we should backspace on 10 and read the next lunit
cRJB...special case involves igo=4 and igop=8 w/ nous=0
cRJB...          and we may write the residual and then write ppflag=9
cRJB...           or we may write directly the ppflag=9

cRJB...special case involves igo=4 and igop=5 w/ nous=0
cRJB...   read lunitm to see if we have an ppflag 8 (to be followed by ppflag 9)
cRJB...   if so, then proceed as in igo=4 and igop=8
             if(igo.eq.4 .and. igop.eq.5)then
              read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
              nrec(lunitm)=nrec(lunitm)+1
              j45=index(outstr,"ppflag")
              igop45=0
              nous45=0
              if(j45.gt.0)then
cRJB...looking for ppflag 8
cRJB...we have found a ppflag on lunitm so we will move fort.10
                cigop45(1:1)=outstr(9:9)
                cigop45(2:2)=outstr(10:10)
                cigop45(3:3)=outstr(11:11)
                read(cigop45,*) igop45
                cnous45(1:1)=outstr(12:12)
                cnous45(2:2)=outstr(13:13)
                cnous45(3:3)=outstr(14:14)
                cnous45(4:4)=outstr(15:15)
                cnous45(5:5)=outstr(16:16)
                cnous45(6:6)=outstr(17:17)
                read(cnous45,*) nous45
                read(10,'(a)',err=998,end=997)string
                nr=nr+1
                if(igop45.eq.8 .and. nous45.le.0)then
cRJB...reset igop and proceed as though the sequence was ppflag 4 and ppflag 8
                 igop=igop45
                else
                 do no=1,nous45
                  read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                  nrec(lunitm)=nrec(lunitm)+1
                  write(lunfio0,'(a)')trim(outstr)
                  nw=nw+1
                 enddo
cRJB...a ppflag of 9 may follow
                 read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                 nrec(lunitm)=nrec(lunitm)+1
                 j95=index(outstr,"ppflag")
                 igop95=0
                 nous95=0
                 if(j95.gt.0)then
cRJB...we have found a ppflag on lunitm so we will move fort.10
                  cigop95(1:1)=outstr(9:9)
                  cigop95(2:2)=outstr(10:10)
                  cigop95(3:3)=outstr(11:11)
                  read(cigop95,*) igop95
                  cnous95(1:1)=outstr(12:12)
                  cnous95(2:2)=outstr(13:13)
                  cnous95(3:3)=outstr(14:14)
                  cnous95(4:4)=outstr(15:15)
                  cnous95(5:5)=outstr(16:16)
                  cnous95(6:6)=outstr(17:17)
                  read(cnous95,*) nous95
                  read(10,'(a)',err=998,end=997)string
                  nr=nr+1
                  if(igop95.eq.9 .and. nous95.le.0)then
cRJB...reset igop and proceed through the ppflag 9 sequence
                   igop=igop95
                  else
                   do no=1,nous95
                    read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                    nrec(lunitm)=nrec(lunitm)+1
                    write(lunfio0,'(a)')trim(outstr)
                    nw=nw+1
                   enddo
                  endif ! if(igop95.eq.8 .and. nous95.le.0)
                 endif ! if(j95.gt.0)then
                endif ! if(igop45.eq.8 .and. nous45.le.0)
               endif !if(j45.gt.0)then
              endif !if(igo.eq.4 .and. igop.eq.5)

             if(igop .eq. 8 .and. igo.eq.4)then
cRJB...the 2nd line from fort.10 will tell the tale
              read(10,'(a)',err=998,end=997)outstr
              nr=nr+1
              read(10,'(a)',err=998,end=997)outstr
              nr=nr+1
              jlev=index(outstr,"level ")
              jres=index(outstr,"restricting")
              if(jlev.gt.0)then
               backspace(10)
               nr=nr-1
               backspace(10)
               nr=nr-1
cRJB...1st line of residual
               read(10,'(a)',err=998,end=997)outstr
               nr=nr+1
               write(lunfio0,'(a)')trim(outstr)
               nw=nw+1
cRJB...2nd line of residual
               read(10,'(a)',err=998,end=997)outstr
               nr=nr+1
               write(lunfio0,'(a)')trim(outstr)
               nw=nw+1
cRJB...3rd line of residual
               read(10,'(a)',err=998,end=997)outstr
               nr=nr+1
               write(lunfio0,'(a)')trim(outstr)
               nw=nw+1
              else
               backspace(10)
               nr=nr-1
               backspace(10)
               nr=nr-1
              endif

               read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
               nrec(lunitm)=nrec(lunitm)+1
               j489=index(outstr,"ppflag")
cRJB... we expect to find a ppflag
              if(j489.gt.0)then
cRJB...we have found a ppflag on lunitm so we will move fort.10
               cigop(1:1)=outstr(9:9)
               cigop(2:2)=outstr(10:10)
               cigop(3:3)=outstr(11:11)
               read(cigop,*) igop
               cnous(1:1)=outstr(12:12)
               cnous(2:2)=outstr(13:13)
               cnous(3:3)=outstr(14:14)
               cnous(4:4)=outstr(15:15)
               cnous(5:5)=outstr(16:16)
               cnous(6:6)=outstr(17:17)
               read(cnous,*) nous
cRJB...moving fort.10-if it's a ppflag
               read(10,'(a)',err=998,end=997)string
               nr=nr+1
               j489p=index(string,"ppflag")
               if(j489p.gt.0)then
               else
                backspace(10)
                nr=nr-1
               endif
cRJB...an unusual case: cfl3d prints max residual and max vorticity BEFORE
cRJB...the last block "updating with 3-D 3-factor scheme" wrapup
cRJB...we check fort.10 for the printout 
               if(n.eq.ngridt)then
                 read(10,'(a)',err=998,end=997)string1
                 nr=nr+1
                 read(10,'(a)',err=998,end=997)string2
                 nr=nr+1
                 jmax=index(string2,"max")
                 if(jmax.gt.0)then
cRJB...we have found the max res p/o-lets write until we find a ppflag
                  write(lunfio0,'(a)')trim(string1)
                  nw=nw+1
                  write(lunfio0,'(a)')trim(string2)
                  nw=nw+1
  489             continue
                  read(10,'(a)',err=998,end=997)outstr
                  nr=nr+1
                  jmax=index(outstr,"ppflag")
                  if(jmax.gt.0)then
cRJB...we move past the final ppflag
                  else
                   write(lunfio0,'(a)')trim(outstr)
                   nw=nw+1
                   go to 489
                  endif
                 else
                  backspace(10)
                  nr=nr-1
                  backspace(10)
                  nr=nr-1
                 endif
               endif ! n.eq.ngridt
               if(nous.gt.0)then
                do no=1,nous
                 read(lunit,'(a)',err=998,end=997,iostat=ios)outstr
                 nrec(lunitm)=nrec(lunitm)+1
                 write(lunfio0,'(a)')trim(outstr)
                 nw=nw+1
                enddo
               else
cRJB...surely nous=0 in this sequence is in error
                write(lunfio0,'(" nous489 error-stop ",9i6)')igo,igop,
     .           nous,lunit,nrec(lunitm),10,nr
                write(6,'(" process:nous489error-stop ")')
                ierror=1
                return
               endif
              else
cRJB...a ppflag mismatch backspaces if we have output data to write
               write(lunfio0,'(" ppflag489mismatch-stop ",9i6)')igo,
     .          igop,lunit,nrec(lunitm),10,nr
               write(6,'(" process:ppflag489mismatch-stop ")')
               ierror=1
               return
              endif ! j489.gt.0
             else
cRJB...test
cRJB...blomx avoid backspacing if j95 active
              if(igop95.gt.0 .or. ngood.gt.0)then
              else
               backspace(lunit)
               nrec(lunitm)=nrec(lunitm)-1
               backspace(10)
               nr=nr-1
              endif
             endif ! igop.eq.8 .and. igo.eq.4
            endif ! if(nous.gt.0)then

           else
cRJB...a ppflag mismatch backspaces if we have output data to write
cRJB...otherwise cfl3d entered writ_buf w/ no data to write
            write(lunfio0,'(" ppflag3mismatch-stop ",9i6)')igo,igop,
     .       lunit,nrec(lunitm),10,nr
            write(6,'(" process:ppflag3mismatch-stop ")')
            ierror=1
            return
           endif ! j3.le.0
          endif ! l27done.le.0
cRJB...COMPLETE THE TEST FOR THE SPECIAL CASES
          l27done=0
          else ! igo .eq.igop
cRJB...a ppflag mismatch
cRJB...entered with a ppflag value of igo from fort.10 and process read 
cRJB...a different value of igop on the lunit
cRJB...process should keep reading the lunits to find the proper igop
cRJB...a ppflag mismatch should backspace on both fort.10 and lunit
cRJB...otherwise cfl3d entered writ_buf w/ no data to write
              if(nous.gt.0)then
cRJB...are these 3 OK
               backspace(lunit)
               nrec(lunitm)=nrec(lunitm)-1
               read(lunit,'(a)')outstr
               backspace(lunit)
              else
               backspace(10)
               nr=nr-1
               backspace(lunit)
               nrec(lunitm)=nrec(lunitm)-1
              endif
          endif ! igo .eq.igop

         endif ! j2.gt.0

  997    continue
         if(nrec0.eq.nrec(lunitm))then
          lunitc(lunitm)=1
         endif
         endif ! lunitc(lunitm).le.0

        enddo
        nlunitc=0
        do lun=1,ngridt
         lunit=lunfio0+lun
         lunitm=lunit-lunfio0
         nlunitc=nlunitc+lunitc(lunitm)
        enddo


cRJB...decisions, decisions
cRJB...we have just completed reading all the lunits for ppflag igo
cRJB...we should go back to reass...
cRJB...some special ppflags and cases make us pause before we do so...
cRJB...a sequence of 4,8,9 may involve printing bc
        if(igo.eq.4 .and. igop .eq. 9)then
cRJB...ppflag, then a blank line, then a 5star
         read(10,'(a)',err=998,end=999)outstr
         nr=nr+1
         read(10,'(a)',err=998,end=999)outstr
         nr=nr+1
         jx=index(outstr(1:5),'*****')
         if(jx.gt.0)then
          backspace(10)
          nr=nr-1
          backspace(10)
          nr=nr-1
          return
cRJB...boundary conditions can follow type 9
         else
          backspace(10)
          nr=nr-1
          backspace(10)
          nr=nr-1
          igo=2
          goto 100 
         endif
cRJB...igop =2 may involve printing bc or  may involve an 5star ENDING header printout
cRJB...it will be a blank line and a line w/ 5star
cRJB...read the next line
cRJB...if the 5star test here fails print the bc
        elseif(igo.eq.4 .and. igop .eq. 2)then
         if(jres.gt.0)then
cRJB...5star test
           read(10,'(a)',err=998,end=999)outstr
           nr=nr+1
           read(10,'(a)',err=998,end=999)outstr
           nr=nr+1
           jx=index(outstr(1:5),'*****')
           if(jx.gt.0)then
cRJB...5star test passes-return
             backspace(10)
             nr=nr-1
             backspace(10)
             nr=nr-1
             return
           else
cRJB...5star test fails-print bc
            backspace(10)
            nr=nr-1
            backspace(10)
            nr=nr-1
            igo=2
            goto 100 
           endif
         else
          igo=2
          goto 100 
         endif
        else
         return
        endif ! if(igo.eq.4 .and. igop .eq. 9
      endif ! igo eq 1
  998 continue
      write(lunfio0,'(" process:error at nr =",9i6)')igo,igop,lunit,
     .   nrec(lunitm),10,nr
      write(6,'(" process:error ")')
      ierror=1
      return
  999 continue
      write(lunfio0,'(" process:  end at nr =",9i6)')igo,igop,lunit,
     .   nrec(lunitm),10,nr
      write(6,'(" process:end ")')
      return
      end
#else
      subroutine dummyreass
c**********************************************************************
c     Purpose: Provide a dummy routine to compile if FASTIO option is
c     not installed
c**********************************************************************
      return
      end
#   endif
