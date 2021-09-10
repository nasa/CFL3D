c
c********************************************************************
c
c     Reads in plot3d grid and q files for a wing
c     surface as well as an experimental data file. 
c     Calculates (-)Cp vs x/c distributions for the
c     spanwise stations dictated in the exp. data file.
c     Outputs the results to a (formatted) tecplot file
c     for plotting. 
c
c     The input plot3d files are expected to be named plot3dg.bin
c     and plot3dq.bin, and the experimental data fille is expected
c     to be named exp.dat
c
c     The experimental data file is assumed to have the following
c     form:
c
c        dummy text line           (title)
c        nsta                      (number of spanwise data stations)
c        dummy text line           (start of new span station)
c        nus, nls, eta             (# upper pts, # lower pts, % span)
c        dummy text                (start of upper surface data)
c        i, xupper, cpupper        (nus of these lines)
c        dummy text line           (start of lower surface data)
c        i, xlower, cplower        (nls of these lines)
c        dummy text line           (start of new span station)
c             .
c             .
c        (repeated for all remaing stations)
c             .
c             .
c
c     case-dependent hardwired indicies...it is assumed j=1 is
c     the wing surface:
c
c     ite1....lower surface te point (chordwise direction - i)
c     ite2....upper surface te point (chordwise direction - i)
c     igrid...zone in the plot3d file that contains the wing surface
c     jroot...wing root point (spanwise direction - j)
c
c********************************************************************
c
      parameter(idim=225,jdim=97,kdim=1,ndim=10)

      dimension x(idim,jdim,kdim,3),q(idim,jdim,kdim,5),xle(jdim),
     .          yle(jdim),zle(jdim),c(jdim),cp(idim,jdim)
      dimension xcmp(idim,ndim),ycmp(idim,ndim),zcmp(idim,ndim),
     .          cpcmp(idim,ndim)
      dimension xexp(idim,ndim),cpexp(idim,ndim),nus(ndim),nls(ndim),
     .          etaexp(ndim)
      dimension id1(ndim),jd1(ndim),kd1(ndim)
c
      character*80 pgrdfn,pqfn,expdat
c
c     write(6,*)'name of plot3d grid file: '
c     read(5,'(a80)') pgrdfn
c     write(6,*)'name of plot3d q file: '
c     read(5,'(a80)') pqfn
c     write(6,*)'name of experimental data file'
c     read(5,'(a80)') expdat
c
      pgrdfn = 'plot3dg.bin_surf'
      pqfn   = 'plot3dq.bin_surf'
      expdat = 'exp.dat'
c
      open(2,file=pgrdfn,status='old',form='unformatted')
      open(3,file=pqfn,status='old',form='unformatted')
      open(4,file=expdat,status='old',form='formatted')
c
c     read and verify header info. in plot3d files
c
      read(2) ngrid
      read(3)
      read(2)(id1(n),jd1(n),kd1(n),n=1,ngrid)
      read(3)
c
      do n=1,ngrid
         write(6,'(''zone '',i3,'' i x j x k = '',3i4)') n,id1(n),
     .   jd1(n),kd1(n)
      end do
c
c     case-dependent hardwired indicies...it is assumed j=1 is
c     the wing surface
c
      igrid  = 3
c     igrid2 = 6
      igrid2 = 3
      ite1  = 1
      ite2  = id1(igrid)
      jroot = 33
c
c     echo parameters
c
      write(6,'(''Cp data to be extracted from zones '',i3,
     .'' and '',i3,'' of '',i3)') igrid,igrid2,ngrid
      write(6,'(''ite1 = '',i3,''  ite2 = ''i3)') ite1,ite2
      write(6,'(''jroot = ''i3)') jroot
c
      idim1 = id1(igrid)
      jdim1 = jd1(igrid)
      kdim1 = kd1(igrid)
      if (idim1 .gt. idim .or. jdim1 .gt. jdim .or.
     .    kdim1 .gt. kdim) then
         write(6,*)'need to increase dimensions'
         stop
      end if
      if (kdim1 .ne. 1) then
         write(6,*)'plot3d file must contain k=const ',
     .   'surface'
         stop
      end if
c
c     read in grid data, concatenate inboard (zone 3) and
c     outboard wing panels into one zone
c
      do n=1,ngrid
         id = id1(n)
         jd = jd1(n)
         kd = kd1(n)
         if (n.eq.igrid .or. n.eq.igrid2) then
            if (n.eq.igrid) then
               jj = 0
            else
               jj = jd1(igrid)
            end if
            read(2)(((x(i,j+jj,k,1),i=1,id),j=1,jd),k=1,kd),
     .             (((x(i,j+jj,k,2),i=1,id),j=1,jd),k=1,kd),
     .             (((x(i,j+jj,k,3),i=1,id),j=1,jd),k=1,kd)
c    .             ,(((idum,i=1,id),j=1,jd),k=1,kd)
         else
            read(2) (((dum,i=1,id),j=1,jd),k=1,kd),
     .              (((dum,i=1,id),j=1,jd),k=1,kd),
     .              (((dum,i=1,id),j=1,jd),k=1,kd) 
c    .             ,(((idum,i=1,id),j=1,jd),k=1,kd)
         end if
      end do
c
      jdim1 = jd1(igrid) + jd1(igrid2)
c
c     read in q data
c
      do n=1,ngrid
         read(3)smach,alp,reue,tm
         id = id1(n)
         jd = jd1(n)
         kd = kd1(n)
         if (n.eq.igrid .or. n.eq.igrid2) then
            if (n.eq.igrid) then
               jj = 0
            else
               jj = jd1(igrid)
            end if
            read(3)((((q(i,j+jj,k,nx),i=1,id),j=1,jd),k=1,kd),nx=1,5)
         else
            read(3)((((dum,i=1,id),j=1,jd),k=1,kd),nx=1,5)
         end if
      end do
c
c     flip y
c
      do i=1,idim1
         do j=1,jdim1
            x(i,j,1,2) = -x(i,j,1,2)
         end do
      end do
c
c     calculate necessary nondimensional quantities
c
      gamma=1.4
      gamm1=gamma-1
      imax  = ite1
      do i=1,idim1
         do j=1,jdim1
            p       = gamm1*q(i,j,1,5)-gamm1/2*(q(i,j,1,2)**2+
     .                q(i,j,1,3)**2+q(i,j,1,4)**2)/q(i,j,1,1)
            cp(i,j) = (p-1/gamma)/(.5*smach**2)
         end do
      end do
      yroot = x(ite1,jroot,1,2)
      do j=1,jroot
         imax = ite1
         xm = x(ite1,j,1,1)
         do i=ite1,ite2
            x(i,j,1,2) = x(i,j,1,2)-yroot
            if (x(i,j,1,1) .le. xm) then
               xm   = x(i,j,1,1)
               imax = i
            end if
         end do
         xle(j) = x(imax,j,1,1)
         yle(j) = x(imax,j,1,2)
         zle(j) = x(imax,j,1,3)
         c(j)   = x(ite1,j,1,1)-xle(j)
      end do
      span = yle(1)
      if (abs(span) .eq. 0.) then
         write(6,*)'stopping...span = 0'
         stop
      end if
      do j=1,jroot
         yle(j) = yle(j)/span
      end do
      do j=1,jroot
         xte = x(ite1,j,1,1)
         do i=ite1,ite2
            x(i,j,1,1) = 1.+ (x(i,j,1,1)-xte)/c(j)
            x(i,j,1,2) = x(i,j,1,2)/span
         end do
      end do
c
c     read exp. data file to detemine station locations
c
      read(4,*)
      read(4,*) nsta
      do nn=1,nsta
         read(4,*)
         read(4,*) nus(nn),nls(nn),etaexp(nn)
         read(4,*)
         do n=1,nus(nn)
            read(4,*) ii,xexp(n,nn),cpexp(n,nn)
         end do
         read(4,*)
         do n=1,nls(nn)
            read(4,*) ii,xexp(n+nus(nn),nn),cpexp(n+nus(nn),nn)
         end do
      end do 
c
c     output x/c and -cp for plotting
c
      do nn=1,nsta
         eta = etaexp(nn)
         write(77,*)'zone'
         write(77,'(''# Computed Data, eta = '',f6.4)') eta
         do i=ite1,ite2
            jsta = jroot
            do j=jroot-1,1,-1
               if (x(i,j+1,1,2) .le. eta .and.
     .            x(i,j,1,2) .gt. eta) then
                  jsta = j
               end if
            end do
            deta =        (eta - x(i,jsta,1,2))/
     .                    (x(i,jsta+1,1,2) - x(i,jsta,1,2))
            xcmp(i,nn)  = x(i,jsta,1,1) + 
     .                    deta*(x(i,jsta+1,1,1)-x(i,jsta,1,1))
            ycmp(i,nn)  = x(i,jsta,1,2) + 
     .                    deta*(x(i,jsta+1,1,2)-x(i,jsta,1,2))
            zcmp(i,nn)  = x(i,jsta,1,3) + 
     .                    deta*(x(i,jsta+1,1,3)-x(i,jsta,1,3))
            cpcmp(i,nn) = cp(i,jsta) + 
     .                    deta*(cp(i,jsta+1)-cp(i,jsta))
            cpcmp(i,nn) = -cpcmp(i,nn)
            write (77,*)xcmp(i,nn),cpcmp(i,nn)
         end do
         write(77,*)'zone'
         write(77,'(''# Experimental Data, eta = '',f6.4)') eta
         do i=1,nus(nn)+nls(nn)
            cpexp(i,nn) = -cpexp(i,nn)
            write (77,*)xexp(i,nn),cpexp(i,nn)
         end do
      end do
c
      close(2)
      close(3)
      close(4)
c
      write(6,'(''Done...tecplot file written to unit 77 (fort.77)'')')
c
      stop
      end
