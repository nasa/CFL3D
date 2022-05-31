      program postproc
c
c   extracts cp and from cfl3d.prout file, specifically for the
c   backward facing step test case. the cfl3d.prout file must be
c   be generated with NPRINT=2 and:
c
c PRINT OUT:
c  GRID IPTYPE ISTART   IEND   IINC JSTART   JEND   JINC KSTART   KEND   KINC
c     1      0      1      1      1      0      0      0      1      1      1
c     2      0      1      1      1      0      0      0      1      1      1
c
c
      open(21,file='cfl3d.prout',form='formatted',status='old')
c
c     first read and write cp data
c
      open(22,file='cp.dat',form='formatted',status='unknown')
      write(22,'(''# X,Cp'')')
      nskip = 10
      do n=1,nskip
         read(21,*)
      end do
      nread = 25
      do n=1,nread
        read(21,*) idum,idum,idum,x,y,dum,dum,dum,dum,dum,dum,dum,cp
        write(22,'(2e15.5)') x,cp
      enddo
      nskip = 31
      do n=1,nskip
         read(21,*)
      end do
      nread = 129
      do n=1,nread
        read(21,*) idum,idum,idum,x,y,dum,dum,dum,dum,dum,dum,dum,cp
        write(22,'(2e15.5)') x,cp
      enddo
      write(6,'('' Output to cp.dat'')')
c
c     now read and write cf data
c
      rewind(21)
      open(23,file='cf.dat',form='formatted',status='unknown')
      write(23,'(''# X,Cf'')')
      nskip = 37
      do n=1,nskip
         read(21,*)
      end do
      nread = 23
      do n=1,nread
        read(21,*) idum,idum,idum,x,y,dum,dum,dum,dum,cf
        write(23,'(2e15.5)') x,cf
      enddo
      nskip = 137
      do n=1,nskip
         read(21,*)
      end do
      nread = 127
      do n=1,nread
        read(21,*) idum,idum,idum,x,y,dum,dum,dum,dum,cf
        write(23,'(2e15.5)') x,cf
      enddo
      write(6,'('' Output to cf.dat'')')
      stop
      end
