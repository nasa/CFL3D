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
      subroutine histout(ihstry,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .                   cmxw,cmyw,cmzw,
     .                   n_clcd,clcd,nblocks_clcd,blocks_clcd,
     .                   chdw,swetw,fmdotw,cfttotw,
     .                   cftmomw,cftpw,cftvw,rmstr,nneg,
     .                   ncycmax,aehist,aesrfdat,nmds,maxaes,timekeep,
     .                   nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Output convergence history for mean-flow equations and
c               turbulence equations. For aeroelastic cases, output
c               generalized displacement, velocity, and force.
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*80 filename, num
      character*241 avgq2,avgq2pert,clcds,clcdp,output_dir
      integer blocks_clcd
      dimension rms(ncycmax),clw(ncycmax),
     .          cdw(ncycmax),cdpw(ncycmax),cdvw(ncycmax),
     .          cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          chdw(ncycmax),swetw(ncycmax),
     .          fmdotw(ncycmax),cfttotw(ncycmax),
     .          cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .          rmstr(ncycmax,nummem),
     .          nneg(ncycmax,nummem),timekeep(ncycmax),
     .          aehist(ncycmax,3,nmds,maxaes), aesrfdat(5,maxaes),
     .          clcd(2,n_clcd,ncycmax), blocks_clcd(2,nblocks_clcd)
c
      common /filenam2/ avgq2,avgq2pert,clcds,clcdp,output_dir
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /igrdtyp/ ip3dgrd,ialph
      common /conversion/ radtodeg
      common /maxiv/ ivmx
      common /elastic/ ndefrm,naesrf
      common /reyue/ reue,tinf,ivisc(3)
      common /avgdata/ xnumavg,iteravg,xnumavg2,ipertavg,iclcd,isubit_r
c
      alphw = radtodeg*alpha
c
c     output convergence history for mean-flow equations
c
      write(12,2) (real(title(i)),i=1,20)
    2 format(2h# ,20a4)
c
      write(12,'(''# Mach='',e12.4,'', alpha='',e12.4,
     . '', ReUe='',e12.4)') real(xmach),real(alphw),real(reue)
      if (real(dt) .gt. 0) then
        write(12,'(''# Final res='',e12.4,
     .   ''   Final time (for unsteady)='',f10.4)') 
     .    real(rms(ntt)),real(time)
      else
        write(12,'(''# Final res='',e12.4)') real(rms(ntt))
      end if
c
c     mean-flow convergence history file options:
c     ihstry = 0....standard: a) ialph = 0 residual,cl,cd,cy,cmy
c                             b) ialph = 1 residual,cl,cd,cz,cmz
c            = 1....control surface: residual,mass flow,pressure force,
c                   viscous force,thrust (momentum) force 
c                   (forces are sums in x+y+z directions)
c            = 2....enhanced standard: ALL force/moment coefficients
c                   are output (thus there is no need to distinguish
c                   between ialpha=0 and ialpha=1)
c 
      if (ihstry.eq.0) then
         if (ialph .eq.0) then
            write(12,'(''# Final cl,cd,cy,cmy='',4e12.4)')
     .      real(clw(ntt)),real(cdw(ntt)),real(cyw(ntt)),
     .      real(cmyw(ntt))
            if (iclcd .ne. 1 .and. iclcd .ne. 2) then
               write(12,'(''#'',i5,'' it    log(res)        cl'',
     .            ''            cd            cy            cmy'')') ntt
            else
               write(12,'(''variables = it    log(res)        cl'',
     .            ''            cd            cz            cmz'')')
            end if
            do 3000 n=1,ntt
            if(real(rms(n)) .le. 0.) rms(n)=1.
c           write(12,'(3x,i6,5e14.5)') n,log10(real(rms(n))),
            write(12,'(3x,i6,6e15.7)') n,log10(real(rms(n))),
     .      real(clw(n)),real(cdw(n)),real(cyw(n)),real(cmyw(n))
 3000       continue
         else
            write(12,'(''# Final cl,cd,cz,cmz='',4e12.4)')
     .      real(clw(ntt)),real(cdw(ntt)),real(czw(ntt)),
     .      real(cmzw(ntt))
            if (iclcd .ne. 1 .and. iclcd .ne. 2) then
               write(12,'(''#'',i5,'' it    log(res)        cl'',
     .            ''            cd            cy            cmy'')') ntt
            else
               write(12,'(''variables = it    log(res)        cl'',
     .            ''            cd            cz            cmz'')')
            end if
            do 3001 n=1,ntt
            if(real(rms(n)) .le. 0.) rms(n)=1.
c           write(12,'(3x,i6,5e14.5)') n,log10(real(rms(n))),
            write(12,'(3x,i6,6e15.7)') n,log10(real(rms(n))),
     .      real(clw(n)),real(cdw(n)),real(czw(n)),real(cmzw(n))
 3001       continue
         end if

            if (iclcd .eq. 1 .or. iclcd .eq. 2) then
               do nn = 1, n_clcd
                  write(num,'(I2)') nn
                  num = adjustl(num)
                  filename = clcdp(1:len_trim(clcdp)) //
     $                 num(1:len_trim(num)) // ".dat"
                  open(unit=151,file=filename,form="formatted",
     .                 status="unknown") 
                  write(151,'(''variables=it,Cl,Cd,Clt,Cdt,Cdp,Cdv'')')
                  do n=1,ntt       
                     if (real(clcd(1,nn,n)) < 1e20) then
                        write(151,'(3x,i6,6e15.7)') n, 
     .                       real(clcd(1,nn,n)), real(clcd(2,nn,n)),
     .                       real(clw(n)), real(cdw(n)),
     .                       real(cdpw(n)), real(cdvw(n)) 
                     end if
                  enddo            
                  close(151)
               end do
            end if
      else if (ihstry.eq.1) then
         write(12,'(''# Final mass_flow,cftp,cftv,cftmom='',4e12.4)')
     .   real(fmdotw(ntt)),real(cftpw(ntt)),real(cftvw(ntt)),
     .   real(cftmomw(ntt))
         if (iclcd .ne. 1 .and. iclcd .ne. 2) then
            write(12,'(''#'',i5,'' it    log(res)      mass_flow'',
     .           ''       cftp           cftv        cftmom '')') ntt
         else
            write(12,'(i6,'' it    log(res)      mass_flow'',
     .           ''       cftp           cftv        cftmom '')') ntt
         end if
         do 3100 n=1,ntt
         if(real(rms(n)) .le. 0.) rms(n)=1.
c        write(12,'(3x,i6,5e14.5)') n,log10(real(rms(n))),
         write(12,'(3x,i6,6e15.7)') n,log10(real(rms(n))),
     .   real(fmdotw(n)),real(cftpw(n)),real(cftvw(n)),real(cftmomw(n))
 3100    continue
      else
         write(12,'(''# Final cl,cd       ='',2e13.5)')
     .   real(clw(ntt)),real(cdw(ntt))
         write(12,'(''# Final cx,cy,cz    ='',3e13.5)')
     .   real(cxw(ntt)),real(cyw(ntt)),real(czw(ntt))
         write(12,'(''# Final cmx,cmy,cmz ='',3e13.5)')
     .   real(cmxw(ntt)),real(cmyw(ntt)),real(cmzw(ntt))
         write(12,*)
         write(12,'(''#    it     log(res)           cl'',
     .   ''           cd           cx           cy'',
     .   ''           cz          cmx          cmy'',
     .   ''          cmz'')')
         write(12,*)
         do n=1,ntt
            if(real(rms(n)) .le. 0.) rms(n)=1.
            write(12,'(i6,9e13.5)') n,log10(real(rms(n))),
     .      real(clw(n)),real(cdw(n)),real(cxw(n)),real(cyw(n)),
     .      real(czw(n)),real(cmxw(n)),real(cmyw(n)),real(cmzw(n))
         end do
      end if
c
c     output convergence history for field equation turb. model
c
      if (ivmx.gt.2) then
c
         write(13,2) (real(title(i)),i=1,20)
         write(13,'(''# Mach='',e12.4,'', alpha='',e12.4,
     .    '', ReUe='',e12.4)') real(xmach),real(alphw),real(reue)
c
         do l=1,nummem
           write(13,'(''# Final turres'',i2,''='',e12.4)') 
     +       l,real(rmstr(ntt,l))
         enddo
         if (real(dt) .gt. 0) then
         write(13,'(''# Final time (for unsteady)='',f10.4)')real(time)
         end if
         if (iclcd .ne. 1 .and. iclcd .ne. 2) then
            if (nummem .eq. 2) then
            write(13,'(''#'',i6,'' it  log(turres1)  log(turres2)'',
     .           ''  nneg1  nneg2'')') ntt
            else if (nummem .eq. 3) then
            write(13,'(''#'',i6,'' it  log(turres1)  log(turres2)'',
     .           ''  log(turres3)  nneg1  nneg2  nneg3'')') ntt
            else if (nummem .eq. 4) then
            write(13,'(''#'',i6,'' it  log(turres1)  log(turres2)'',
     .           ''  log(turres3)  log(turres4)  nneg1'',
     .           ''  nneg2  nneg3  nneg4'')') ntt
            else
            write(13,'(''#'',i6,'' it  log(turres1)  log(turres2)'',
     .           ''  log(turres3)  log(turres4)  log(turres5)'',
     .           ''  log(turres6)  log(turres7)  nneg1  nneg2'',
     .           ''  nneg3  nneg4  nneg5  nneg6  nneg7'')') ntt
            end if
         else
            if (nummem .eq. 2) then
            write(13,'(''variables= it  log(turres1)  log(turres2)'',
     .           ''  nneg1  nneg2'')')
            else if (nummem .eq. 3) then
            write(13,'(''variables= it  log(turres1)  log(turres2)'',
     .           ''  log(turres3)  nneg1  nneg2  nneg3'')')
            else if (nummem .eq. 4) then
            write(13,'(''variables= it  log(turres1)  log(turres2)'',
     .           ''  log(turres3)  log(turres4)  nneg1'',
     .           ''  nneg2  nneg3  nneg4'')')
            else
            write(13,'(''variables= it  log(turres1)  log(turres2)'',
     .           ''  log(turres3)  log(turres4)  log(turres5)'',
     .           ''  log(turres6)  log(turres7)  nneg1  nneg2'',
     .           ''  nneg3  nneg4  nneg5  nneg6  nneg7'')')
            end if
         end if
         do 3200 n=1,ntt
           do l=1,nummem
             if(real(rmstr(n,l)) .le. 0.) rmstr(n,l)=1.
           enddo
         if (nummem .eq. 2) then
         write(13,'(3x,i6,2e14.5,1x,i6,1x,i6)') n,
     .   log10(real(rmstr(n,1))),log10(real(rmstr(n,2))),nneg(n,1),
     .   nneg(n,2)
         else if (nummem .eq. 3) then
         write(13,'(3x,i6,3e14.5,1x,i6,1x,i6,1x,i6)') n,
     .   log10(real(rmstr(n,1))),log10(real(rmstr(n,2))),
     .   log10(real(rmstr(n,3))),nneg(n,1),
     .   nneg(n,2),nneg(n,3)
         else if (nummem .eq. 4) then
         write(13,'(3x,i6,4e14.5,1x,i6,1x,i6,1x,i6,1x,i6)') n,
     .   log10(real(rmstr(n,1))),log10(real(rmstr(n,2))),
     .   log10(real(rmstr(n,3))),log10(real(rmstr(n,4))),
     .   nneg(n,1),nneg(n,2),nneg(n,3),nneg(n,4)
         else
         write(13,'(3x,i6,7e14.5,1x,i6,1x,i6,1x,i6,
     .   1x,i6,1x,i6,1x,i6,1x,i6)') n,
     .   log10(real(rmstr(n,1))),log10(real(rmstr(n,2))),
     .   log10(real(rmstr(n,3))),log10(real(rmstr(n,4))),
     .   log10(real(rmstr(n,5))),log10(real(rmstr(n,6))),
     .   log10(real(rmstr(n,7))),
     .   nneg(n,1),nneg(n,2),nneg(n,3),nneg(n,4),
     .   nneg(n,5),nneg(n,6),nneg(n,7)
         end if
 3200    continue
c
      end if 
c
c     output history of aeroelastic data 
c
      if (iunst.gt.1 .and. naesrf.gt.0) then
         write(34,2) (real(title(i)),i=1,20)
         write(34,'('' Mach='',e12.4,'', alpha='',e12.4,
     .    '', ReUe='',e12.4)') real(xmach),real(alphw),real(reue)
         write(34,'('' Number of aeroelastic surfaces ='',i3)') naesrf
         do iaes=1,naesrf
             write(34,'('' Data for aeroelastic surface '',i3)') iaes
             nmodes = aesrfdat(5,iaes)
             do nm=1,nmodes
                write(34,'(''   mode number'',i4)') nm
                write(34,'(''       it      time        xs(2*n-1)'',
     .             ''      xs(2*n)     gforcn(2*n)'')')
                do n=1,ntt
c                  write(34,'(3x,i6,5e14.5)') n,real(timekeep(n)),
                   write(34,'(3x,i6,6e15.7)') n,real(timekeep(n)),
     .             real(aehist(n,1,nm,iaes)),real(aehist(n,2,nm,iaes)),
     .             real(aehist(n,3,nm,iaes))
                end do
             end do
         end do
      end if
c
      return
      end

