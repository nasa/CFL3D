      program cpget_somellor
c
c   extracts cp from cfl3d.prout file
c   This program is HARDWIRED for the specific So Mellor
c   testcases somellor_sa.inp and somellor_sarc.inp
c
      open(2,file='cfl3d.prout',form='formatted',status='old')
      open(3,file='cp.dat',form='formatted',status='unknown')
c  read front matter that we do not need
      do n=1,10
        read(2,*)
      enddo
      sswall=0.
c  read cp
      write(3,'(''variables="s","cp"'')')
      do n=1,257
        read(2,*) idum,idum,idum,x,y,z,dum,dum,dum,dum,dum,dum,cp
        if (n .eq. 1) then
          dist=x
          cpref=cp
        else
          dist=sqrt((x-xold)**2+(z-zold)**2)
        end if
        sswall=sswall+dist
        write(3,'(2e15.5)') sswall,cp-cpref
        xold=x
        zold=z
      enddo
      write(6,'('' Output (s vs cp) to cp.dat'')')
      stop
      end
