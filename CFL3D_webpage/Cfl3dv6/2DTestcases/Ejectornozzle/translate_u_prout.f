      program translate_u_prout
c
c   read cfl3d.prout file and write u-vel files for plotting
c   This code is HARDWIRED for this particular case!
c
      dimension u(121),y(121)
c
      open (21,file='cfl3d.prout',form='formatted',status='old')
c
      open (31,file='ucfd.x3.dat',form='formatted',status='unknown')
      open (32,file='ucfd.x5.dat',form='formatted',status='unknown')
      open (33,file='ucfd.x7.dat',form='formatted',status='unknown')
      open (34,file='ucfd.x10.5.dat',form='formatted',status='unknown')
c
      write(31,'(''variables="y/H","axial velocity (ft/s)"'')')
      write(32,'(''variables="y/H","axial velocity (ft/s)"'')')
      write(33,'(''variables="y/H","axial velocity (ft/s)"'')')
      write(34,'(''variables="y/H","axial velocity (ft/s)"'')')
c
      do n=1,10
        read(21,*)
      enddo
      do n=1,121
        read(21,*) i,j,k,x,z,y(n),u(n)
c       re-scale u to be in ft/s:
        u(n)=u(n)*1117.*.22
      enddo
      do n=121,2,-1
        write(31,'(2e16.6)') -y(n)/y(121),u(n)
      enddo
      do n=1,121
        write(31,'(2e16.6)') y(n)/y(121),u(n)
      enddo
c
      do n=1,10
        read(21,*)
      enddo
      do n=1,121
        read(21,*) i,j,k,x,z,y(n),u(n)
c       re-scale u to be in ft/s:
        u(n)=u(n)*1117.*.22
      enddo
      do n=121,2,-1
        write(32,'(2e16.6)') -y(n)/y(121),u(n)
      enddo
      do n=1,121
        write(32,'(2e16.6)') y(n)/y(121),u(n)
      enddo
c
      do n=1,10
        read(21,*)
      enddo
      do n=1,121
        read(21,*) i,j,k,x,z,y(n),u(n)
c       re-scale u to be in ft/s:
        u(n)=u(n)*1117.*.22
      enddo
      do n=121,2,-1
        write(33,'(2e16.6)') -y(n)/y(121),u(n)
      enddo
      do n=1,121
        write(33,'(2e16.6)') y(n)/y(121),u(n)
      enddo
c
      do n=1,10
        read(21,*)
      enddo
      do n=1,121
        read(21,*) i,j,k,x,z,y(n),u(n)
c       re-scale u to be in ft/s:
        u(n)=u(n)*1117.*.22
      enddo
      do n=121,2,-1
        write(34,'(2e16.6)') -y(n)/y(121),u(n)
      enddo
      do n=1,121
        write(34,'(2e16.6)') y(n)/y(121),u(n)
      enddo
c
      write(6,'('' output to ucfd.x3.dat, ucfd.x5.dat, '',
     + ''ucfd.x7.dat, ucfd.x10.5.dat'')')
c
      stop
      end
