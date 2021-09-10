      program cfget_somellor
c
c   extracts cf from cfl3d.prout file
c   This program is HARDWIRED for the specific So Mellor
c   testcases somellor_sa.inp and somellor_sarc.inp
c
      open(2,file='cfl3d.prout',form='formatted',status='old')
      open(3,file='cf.dat',form='formatted',status='unknown')
c  read front matter that we do not need
      do n=1,269
        read(2,*)
      enddo
      sswall=0.
c  read cf
      write(3,'(''variables="s","cf"'')')
      do n=1,255
        read(2,*) idum,idum,idum,x,y,z,dum,dum,dum,cf
        if (n .eq. 1) then
          dist=x
        else
          dist=sqrt((x-xold)**2+(z-zold)**2)
        end if
        sswall=sswall+dist
        if (cf .lt. 0.) cf=-cf
        write(3,'(2e15.5)') sswall,cf
        xold=x
        zold=z
      enddo
      write(6,'('' tecplot file output to cf.dat'')')
      stop
      end
