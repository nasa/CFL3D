      subroutine histout(ihstry,rms,clw,cdw,cdpw,cdvw,cxw,cyw,czw,
     .                   cmxw,cmyw,cmzw,chdw,swetw,fmdotw,cfttotw,
     .                   cftmomw,cftpw,cftvw,rmstr1,rmstr2,nneg1,nneg2,
     .                   ncycmax,aehist,aesrfdat,nmds,maxaes,timekeep)
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
      dimension rms(ncycmax),clw(ncycmax),
     .          cdw(ncycmax),cdpw(ncycmax),cdvw(ncycmax),
     .          cxw(ncycmax),cyw(ncycmax),czw(ncycmax),
     .          cmxw(ncycmax),cmyw(ncycmax),cmzw(ncycmax),
     .          chdw(ncycmax),swetw(ncycmax),
     .          fmdotw(ncycmax),cfttotw(ncycmax),
     .          cftmomw(ncycmax),cftpw(ncycmax),cftvw(ncycmax),
     .          rmstr1(ncycmax),rmstr2(ncycmax),
     .          nneg1(ncycmax),nneg2(ncycmax),timekeep(ncycmax),
     .          aehist(ncycmax,3,nmds,maxaes), aesrfdat(5,maxaes)
c
      common /info/ title(20),rkap(3),xmach,alpha,beta,dt,fmax,nit,ntt,
     .        idiag(3),nitfo,iflagts,iflim(3),nres,levelb(5),mgflag,
     .        iconsf,mseq,ncyc1(5),levelt(5),nitfo1(5),ngam,nsm(5),iipv
      common /unst/ time,cfltau,ntstep,ita,iunst
      common /igrdtyp/ ip3dgrd,ialph
      common /conversion/ radtodeg
      common /maxiv/ ivmx
      common /elastic/ ndefrm,naesrf
      common /reyue/ reue,tinf,ivisc(3)
c
      alphw = radtodeg*alpha
c
c     output convergence history for mean-flow equations
c
      write(12,2) (real(title(i)),i=1,20)
    2 format(20a4)
c
      write(12,'('' Mach='',e12.4,'', alpha='',e12.4,
     . '', ReUe='',e12.4)') real(xmach),real(alphw),real(reue)
      if (real(dt) .gt. 0) then
        write(12,'('' Final res='',e12.4,
     .   ''   Final time (for unsteady)='',f10.4)') 
     .    real(rms(ntt)),real(time)
      else
        write(12,'('' Final res='',e12.4)') real(rms(ntt))
      end if
c
c     mean-flow convergence history file options:
c     ihstry = 0....standard: a) ialph = 0 residual,cl,cd,cy,cmy
c                             b) ialph = 1 residual,cl,cd,cz,cmz
c            = 1....control surface: residual,mass flow,pressure force,
c                   viscous force,thrust (momentum) force 
c                   (forces are sums in x+y+z directions)
c 
      if (ihstry.eq.0) then
         if (ialph .eq.0) then
            write(12,'('' Final cl,cd,cy,cmy='',4e12.4)')
     .      real(clw(ntt)),real(cdw(ntt)),real(cyw(ntt)),
     .      real(cmyw(ntt))
            write(12,'(i6,'' it    log(res)        cl'',
     .      ''            cd            cy            cmy'')') ntt
            do 3000 n=1,ntt
            if(real(rms(n)) .le. 0.) rms(n)=1.
            write(12,'(3x,i6,5e14.5)') n,log10(real(rms(n))),
     .      real(clw(n)),real(cdw(n)),real(cyw(n)),real(cmyw(n))
 3000       continue
         else
            write(12,'('' Final cl,cd,cz,cmz='',4e12.4)')
     .      real(clw(ntt)),real(cdw(ntt)),real(czw(ntt)),
     .      real(cmzw(ntt))
            write(12,'(i6,'' it    log(res)        cl'',
     .      ''            cd            cz            cmz'')') ntt
            do 3001 n=1,ntt
            if(real(rms(n)) .le. 0.) rms(n)=1.
            write(12,'(3x,i6,5e14.5)') n,log10(real(rms(n))),
     .      real(clw(n)),real(cdw(n)),real(czw(n)),real(cmzw(n))
 3001       continue
         end if
      else
         write(12,'('' Final mass_flow,cftp,cftv,cftmom='',4e12.4)')
     .   real(fmdotw(ntt)),real(cftpw(ntt)),real(cftvw(ntt)),
     .   real(cftmomw(ntt))
         write(12,'(i6,'' it    log(res)      mass_flow'',
     .   ''       cftp           cftv        cftmom '')') ntt
         do 3100 n=1,ntt
         if(real(rms(n)) .le. 0.) rms(n)=1.
         write(12,'(3x,i6,5e14.5)') n,log10(real(rms(n))),
     .   real(fmdotw(n)),real(cftpw(n)),real(cftvw(n)),real(cftmomw(n))
 3100    continue

      end if
c
c     output convergence history for field equation turb. model
c
      if (ivmx.gt.2) then
c
         write(13,2) (real(title(i)),i=1,20)
         write(13,'('' Mach='',e12.4,'', alpha='',e12.4,
     .    '', ReUe='',e12.4)') real(xmach),real(alphw),real(reue)
c
         write(13,'('' Final turres1='',e12.4)') real(rmstr1(ntt))
         write(13,'('' Final turres2='',e12.4)') real(rmstr2(ntt))
         if (real(dt) .gt. 0) then
           write(13,'('' Final time (for unsteady)='',f10.4)')real(time)
         end if
         write(13,'(i6,'' it  log(turres1)  log(turres2)'',
     .         ''  nneg1  nneg2'')') ntt
         do 3200 n=1,ntt
         if(real(rmstr1(n)) .le. 0.) rmstr1(n)=1.
         if(real(rmstr2(n)) .le. 0.) rmstr2(n)=1.
         write(13,'(3x,i6,2e14.5,1x,i6,1x,i6)') n,
     .   log10(real(rmstr1(n))),log10(real(rmstr2(n))),nneg1(n),
     .   nneg2(n)
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
                   write(34,'(3x,i6,5e14.5)') n,real(timekeep(n)),
     .             real(aehist(n,1,nm,iaes)),real(aehist(n,2,nm,iaes)),
     .             real(aehist(n,3,nm,iaes))
                end do
             end do
         end do
      end if
c
      return
      end
