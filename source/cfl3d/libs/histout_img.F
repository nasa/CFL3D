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
      subroutine histout_img(ihstry,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .                   cmxw,cmyw,cmzw,chdw,swetw,fmdotw,cfttotw,
     .                   cftmomw,cftpw,cftvw,rmstr,nneg,
     .                   ncycmax,aehist,aesrfdat,nmds,maxaes,timekeep,
     .                   nummem)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Output convergence history for derivatives of the
c               mean-flow equations.
c               The derivative is given by Im(z)/h, where h is the
c               magnitude of the perturbation in the complex plane,
c               and z is a solution quantity (e.g. Cl)
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension rms(ncycmax),clw(ncycmax),
     .          cdw(ncycmax),cdpw(ncycmax),cdvw(ncycmax),
     .          cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          chdw(ncycmax),swetw(ncycmax),
     .          fmdotw(ncycmax),cfttotw(ncycmax),
     .          cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .          rmstr(ncycmax,nummem),
     .          nneg(ncycmax,nummem),timekeep(ncycmax),
     .          aehist(ncycmax,3,nmds,maxaes), aesrfdat(5,maxaes)
c
      common /complx/ xmach_img,alpha_img,beta_img,reue_img,tinf_img,
     .                geom_img,surf_img,xrotrate_img,yrotrate_img,
     .                zrotrate_img
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax
      common /igrdtyp/ ip3dgrd,ialph
      common /conversion/ radtodeg
      common /maxiv/ ivmx
      common /elastic/ ndefrm,naesrf
      common /reyue/ reue,tinf,ivisc(3)
c
      alphw = radtodeg*alpha
c
      write(35,*)
      write(35,'(''***********************************************'',
     .  ''********************************'')')
c
c     determine which variable the derivative is with respect to
c
      iunit = 35
      call gradinfo(delh,iunit)
c
      write(35,'(''                        with step size = '',
     . e12.5)') real(delh)
      write(35,'(''***********************************************'',
     .  ''********************************'')')
      write(35,*)
c
c     output derivative convergence history for mean-flow equations
c
      write(35,2) (real(title(i)),i=1,20)
    2 format(20a4)
c
      write(35,'('' Mach='',e12.4,'', alpha='',e12.4,
     . '', ReUe='',e12.4)') real(xmach),real(alphw),real(reue)
      if (real(dt) .gt. 0) then
        write(35,'('' Final Im(res)='',e12.4,
     .   ''   Final time (for unsteady)='',f10.4)')
     .    imag(rms(ntt)),real(time)
      else
        write(35,'('' Final Im(res)='',e12.4)') imag(rms(ntt))
      end if
c
c     mean-flow derivative convergence history file options:
c     ihstry = 0....standard: a) ialph = 0 Im(residual),d(cl)/d(),
c                                          d(cd)/d(),d(cy)/d(),d(cmy)/d()
c                             b) ialph = 1 Im(residual),d(cl)/d(),
c                                          d(cd)/d(),d(cz)/d(),d(cmz)/d()
c            = 1....control surface: Im(residual),d(mass flow)/d(),
c                   d(pressure force)/d(),d(viscous force)/d(),
c                   d(viscous force)/d(),d(thrust (momentum) force)/d()
c                   (forces are sums in x+y+z directions)
c            = 2....enhanced standard: ALL force/moment derivative
c                   coefficients are output (thus there is no need
c                   to distinguish between ialpha=0 and ialpha=1)
c
c           where d() denotes the derivative with respect to whatever 
c           variable has been given an imaginary component on input 
c 
      if (ihstry.eq.0) then
         if (ialph .eq.0) then
            write(35,'('' Final d/d() of cl,cd,cy,cmy='',4e12.4)')
     .      imag(clw(ntt))/real(delh),imag(cdw(ntt))/real(delh),
     .      imag(cyw(ntt))/real(delh),imag(cmyw(ntt))/real(delh)
            write(35,'(i6,'' it  log|Im(res)|    d(cl)/d()'',
     .      ''     d(cd)/d()     d(cy)/d()     d(cmy)/d()'')') ntt
            do 3000 n=1,ntt
            if(imag(rms(n)) .eq. 0.) rms(n)=cmplx(real(rms(n)),1.)
            write(35,'(3x,i6,5e14.5)') n,log10(abs(imag(rms(n)))),
     .      imag(clw(n))/real(delh),imag(cdw(n))/real(delh),
     .      imag(cyw(n))/real(delh),imag(cmyw(n))/real(delh)
 3000       continue
         else
            write(35,'('' Final d/d() of cl,cd,cz,cmz='',4e12.4)')
     .      imag(clw(ntt))/real(delh),imag(cdw(ntt))/real(delh),
     .      imag(czw(ntt))/real(delh),imag(cmzw(ntt))/real(delh)
            write(35,'(i6,'' it  log|Im(res)|    d(cl)/d()'',
     .      ''     d(cd)/d()     d(cz)/d()     d(cmz)/d()'')') ntt
            do 3001 n=1,ntt
            if(imag(rms(n)) .eq. 0.) rms(n)=cmplx(real(rms(n)),1.)
            write(35,'(3x,i6,5e14.5)') n,log10(abs(imag(rms(n)))),
     .      imag(clw(n))/real(delh),imag(cdw(n))/real(delh),
     .      imag(czw(n))/real(delh),imag(cmzw(n))/real(delh)
 3001       continue
         end if
      else if (ihstry.eq.1) then
         write(35,'('' Final d/d() of mdot,cftp,cftv,'',
     .   ''cftmom='',4e12.4)')
     .   imag(fmdotw(ntt))/real(delh),imag(cftpw(ntt))/real(delh),
     .   imag(cftvw(ntt))/real(delh),imag(cftmomw(ntt))/real(delh)
         write(35,'(i6,'' it  log|Im(res)|   d(mdot)/d()'',
     .   ''    d(cftp/d()   d(cftv)/d() d(cftmom)/d()'')') ntt
         do 3100 n=1,ntt
         if(imag(rms(n)) .eq. 0.) rms(n)=cmplx(real(rms(n)),1.)
         write(35,'(3x,i6,5e14.5)') n,log10(abs(imag(rms(n)))),
     .   imag(fmdotw(n))/real(delh),imag(cftpw(n))/real(delh),
     .   imag(cftvw(n))/real(delh),imag(cftmomw(n))/real(delh)
 3100    continue
      else
         write(35,'('' Final d/d() of cl,cd       ='',2e13.5)')
     .   imag(clw(ntt))/real(delh),imag(cdw(ntt))/real(delh)
         write(35,'('' Final d/d() of cx,cy,cz    ='',3e13.5)')
     .   imag(cxw(ntt))/real(delh),imag(cyw(ntt))/real(delh),
     .   imag(czw(ntt))/real(delh)
         write(35,'('' Final d/d() of cmx,cmy,cmz ='',3e13.5)')
     .   imag(cmxw(ntt))/real(delh),imag(cmyw(ntt))/real(delh),
     .   imag(cmzw(ntt))/real(delh)
         write(35,*)
         write(35,'(''    it log|Im(res)|    d(cl)/d()'',
     .   ''    d(cd)/d()    d(cx)/d()    d(cy)/d()    d(cz)/d()'',
     .   ''   d(cmx)/d()   d(cmy)/d()   d(cmz)/d()'')')
         write(35,*)
         do n=1,ntt
            if(imag(rms(n)) .eq. 0.) rms(n)=cmplx(real(rms(n)),1.)
            write(35,'(i6,9e13.5)') n,log10(abs(imag(rms(n)))),
     .      imag(clw(n))/real(delh),imag(cdw(n))/real(delh),
     .      imag(cxw(n))/real(delh),imag(cyw(n))/real(delh),
     .      imag(czw(n))/real(delh),imag(cmxw(n))/real(delh),
     .      imag(cmyw(n))/real(delh),imag(cmzw(n))/real(delh)
         end do
      end if
c
      return
      end
