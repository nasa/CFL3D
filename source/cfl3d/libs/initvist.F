      subroutine initvist(nbl,jdim,kdim,idim,vist3d,tursav,smin,cmuv)
c
c     $Id$
c
c***********************************************************************
c     Purpose:  Initialize the turbulent initial conditions on a mesh
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      dimension tursav(jdim,kdim,idim,2),vist3d(jdim,kdim,idim),
     .          smin(jdim-1,kdim-1,idim-1),cmuv(jdim-1,kdim-1,idim-1)
c
      common /ivals/ p0,rho0,c0,u0,v0,w0,et0,h0,pt0,rhot0,qiv(5),
     .        tur10,tur20
      common /reyue/ reue,tinf,ivisc(3)
c
c   scaling factor for distance-based I.C.s:
      refac=reue/5.e6
c
c     write(15,904)nbl
  904 format(39h turbulent initial conditions for block,i3)
c
      if (ivisc(1).eq.8  .or. ivisc(2).eq.8  .or. ivisc(3).eq.8  .or.
     +    ivisc(1).eq.9  .or. ivisc(2).eq.9  .or. ivisc(3).eq.9  .or.
     +    ivisc(1).eq.13 .or. ivisc(2).eq.13 .or. ivisc(3).eq.13 .or.
     +    ivisc(1).eq.14 .or. ivisc(2).eq.14 .or. ivisc(3).eq.14) then
      do 4020 i=1,idim-1
        do 4020 k=1,kdim-1
          do 4020 j=1,jdim-1
            cmuv(j,k,i)=-.09
 4020 continue
      end if
c   For field eqn turbulence models:
      if (ivisc(1).ge.4 .or. ivisc(2).ge.4 .or. ivisc(3).ge.4) then
        if (ivisc(1).eq.4 .or. ivisc(2).eq.4 .or. ivisc(3).eq.4) then
          do 2000 i=1,idim-1
            do 2000 k=1,kdim-1
              do 2000 j=1,jdim-1
                tursav(j,k,i,1)=tur10
                vist3d(j,k,i)=0.009
 2000     continue
        else if (ivisc(1).eq.5 .or. ivisc(2).eq.5 .or.
     .           ivisc(3).eq.5) then
          do 2001 i=1,idim-1
            do 2001 k=1,kdim-1
              do 2001 j=1,jdim-1
                tursav(j,k,i,1)=tur10
                vist3d(j,k,i)=0.009
 2001     continue
        else if (ivisc(1).eq.11 .or. ivisc(2).eq.11 .or.
     .           ivisc(3).eq.11 .or.
     .           ivisc(1).eq.10 .or. ivisc(2).eq.10 .or.
     .           ivisc(3).eq.10 .or.
     .           ivisc(1).eq. 9 .or. ivisc(2).eq. 9 .or.
     .           ivisc(3).eq. 9 .or. ivisc(1).eq.13 .or.
     .           ivisc(2).eq.13 .or. ivisc(3).eq.13) then
          do 2002 i=1,idim-1
            do 2002 k=1,kdim-1
              do 2002 j=1,jdim-1
                zk1=tur20
                exponent=-471.*(smin(j,k,i)*refac)+.47
                if(real(exponent) .lt. -40.) then
                  zk2=0.
                else
                  zk2=10.**(exponent)
                end if
                exponent=-37.5*(smin(j,k,i)*refac)-3.7
                if(real(exponent) .lt. -40.) then
                  zk3=0.
                else
                  zk3=10.**(exponent)
                end if
                zk4=6.7*(smin(j,k,i)*refac)
                tursav(j,k,i,2)=ccmin(zk2,zk3)
                tursav(j,k,i,2)=ccmax(tursav(j,k,i,2),zk1)
                tursav(j,k,i,2)=ccmin(tursav(j,k,i,2),zk4)
c
                ep1=tur10
                exponent=-555.*(smin(j,k,i)*refac)-6.
                if(real(exponent) .lt. -20.) then
                  ep2=0.
                else
                  ep2=10.**(exponent)
                end if
                exponent=-280.*(smin(j,k,i)*refac)-9.2
                if(real(exponent) .lt. -20.) then
                  ep3=0.
                else
                  ep3=10.**(exponent)
                end if
                exponent=13333.*(smin(j,k,i)*refac)-9.8
                if(real(exponent) .gt. 20.) then
                  ep4=1.e20
                else
                  ep4=10.**(exponent)
                end if
                tursav(j,k,i,1)=ccmin(ep2,ep3)
                tursav(j,k,i,1)=ccmax(tursav(j,k,i,1),ep1)
                tursav(j,k,i,1)=ccmin(tursav(j,k,i,1),ep4)
                vist3d(j,k,i)=.09*rho0*tursav(j,k,i,2)**2/
     .                        tursav(j,k,i,1)

 2002     continue
        else
          const1=45.8
          const2=1.68
          smax=const2/(2.*const1)
          tmax=-const1*smax**2 + const2*smax
          cmu=1.
          if (ivisc(1).eq. 8 .or. ivisc(2).eq. 8 .or.
     .        ivisc(3).eq. 8 .or.
     .        ivisc(1).eq.12 .or. ivisc(2).eq.12 .or.
     .        ivisc(3).eq.12 .or.
     .        ivisc(1).eq.14 .or. ivisc(2).eq.14 .or.
     .        ivisc(3).eq.14) cmu=.09
          v3dset=cmu*rho0*tur20/tur10
          do 2009 i=1,idim-1
            do 2009 k=1,kdim-1
              do 2009 j=1,jdim-1
                zk1=tur20
                zk2=-const1*(smin(j,k,i)*refac)**2 +
     +               const2*smin(j,k,i)*refac
                tursav(j,k,i,2)=ccmax(zk1,zk2)
                v3d=tursav(j,k,i,2)*100./tmax
                v3d=ccmax(v3d,v3dset)
                om1=-12444.*(smin(j,k,i)*refac) + .54
                om2=cmu*tursav(j,k,i,2)/v3d
                tursav(j,k,i,1)=ccmax(om1,om2)
                vist3d(j,k,i)=cmu*rho0*tursav(j,k,i,2)/
     .                        tursav(j,k,i,1)
 2009     continue
        end if
c
      else
c
      do 4021 i=1,idim-1
        do 4021 k=1,kdim-1
          do 4021 j=1,jdim-1
            vist3d(j,k,i)=0.
 4021 continue
c
      end if
c
      return
      end
