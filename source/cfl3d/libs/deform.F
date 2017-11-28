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
      subroutine deform_el(islavept,nslave,nmaster,u,xst,yst,zst,
     .                     nt,myhost,mycomm,myid,nnodes,
     .                     mblk2nd,nblelst,maxbl,iseqr)
c
c***********************************************************************
c     Purpose: Compute mesh deformation via Elastic Equations
c***********************************************************************
c

      integer stats
c
      dimension islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave)
      dimension u(3*nslave)
      dimension mblk2nd(maxbl)
      dimension nblelst(maxbl,2)

      allocatable :: sa(:)
      allocatable :: stiffl(:,:)
      allocatable :: ija(:)
      allocatable :: b(:)
      allocatable :: ei(:)
      allocatable :: ripm(:,:)
      allocatable :: ej(:)
      allocatable :: ek(:)
      allocatable :: gij(:)
      allocatable :: gjk(:)
      allocatable :: gik(:)
      allocatable :: xix(:)
      allocatable :: xiy(:)
      allocatable :: xiz(:)
      allocatable :: etax(:)
      allocatable :: etay(:)
      allocatable :: etaz(:)
      allocatable :: zetax(:)
      allocatable :: zetay(:)
      allocatable :: zetaz(:)
      allocatable :: ooj(:)
      allocatable :: volij(:)
      allocatable :: volik(:)
      allocatable :: tinv(:,:)
      allocatable :: cjm1(:,:)
      allocatable :: cjp1(:,:)
      allocatable :: cim1(:,:)
      allocatable :: cip1(:,:)
      allocatable :: ckm1(:,:)
      allocatable :: ckp1(:,:)
      allocatable :: c(:,:)
      allocatable :: t(:,:)
      allocatable :: ul(:,:)
      allocatable :: n1(:)

c
      common /twod/ i2d
      common /zero/ iexp
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax 
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit

c
      memuse = 0
c
      allocate( sa(245*nslave+2), stat=stats )
      call umalloc(245*nslave+2,0,'sa',memuse,stats)
      allocate( stiffl(24,24), stat=stats )
      call umalloc(576,0,'stiffl',memuse,stats)
      allocate( ija(245*nslave+2), stat=stats )
      call umalloc(245*nslave+2,1,'ija',memuse,stats)
      allocate( b(3*nslave), stat=stats )
      call umalloc(3*nslave,0,'b',memuse,stats)
      allocate( ripm(7,nslave), stat=stats )
      call umalloc(7*nslave,0,'ripm',memuse,stats)
      allocate( ei(nslave), stat=stats )
      call umalloc(nslave,0,'ei',memuse,stats)
      allocate( ej(nslave), stat=stats )
      call umalloc(nslave,0,'ej',memuse,stats)
      allocate( ek(nslave), stat=stats )
      call umalloc(nslave,0,'ek',memuse,stats)
      allocate( gjk(nslave), stat=stats )
      call umalloc(nslave,0,'gjk',memuse,stats)
      allocate( gij(nslave), stat=stats )
      call umalloc(nslave,0,'gij',memuse,stats)
      allocate( gik(nslave), stat=stats )
      call umalloc(nslave,0,'gik',memuse,stats)
      allocate( xix(nslave), stat=stats )
      call umalloc(nslave,0,'xix',memuse,stats)
      allocate( xiy(nslave), stat=stats )
      call umalloc(nslave,0,'xiy',memuse,stats)
      allocate( xiz(nslave), stat=stats )
      call umalloc(nslave,0,'xiz',memuse,stats)
      allocate( etax(nslave), stat=stats )
      call umalloc(nslave,0,'etax',memuse,stats)
      allocate( etay(nslave), stat=stats )
      call umalloc(nslave,0,'etay',memuse,stats)
      allocate( etaz(nslave), stat=stats )
      call umalloc(nslave,0,'etaz',memuse,stats)
      allocate( zetax(nslave), stat=stats )
      call umalloc(nslave,0,'zetax',memuse,stats)
      allocate( zetay(nslave), stat=stats )
      call umalloc(nslave,0,'zetay',memuse,stats)
      allocate( zetaz(nslave), stat=stats )
      call umalloc(nslave,0,'zetaz',memuse,stats)
      allocate( ooj(nslave), stat=stats )
      call umalloc(nslave,0,'ooj',memuse,stats)
c     allocate( volij(nslave), stat=stats )
c     call umalloc(nslave,0,'volij',memuse,stats)
      allocate( volij(3*nslave), stat=stats )
      call umalloc(3*nslave,0,'volij',memuse,stats)
      allocate( volik(nslave), stat=stats )
      call umalloc(nslave,0,'volik',memuse,stats)
      allocate( tinv(6,6), stat=stats )
      call umalloc(36,0,'tinv',memuse,stats)
      allocate( cjm1(6,6), stat=stats )
      call umalloc(36,0,'cjm1',memuse,stats)
      allocate( cjp1(6,6), stat=stats )
      call umalloc(36,0,'cjp1',memuse,stats)
      allocate( cim1(6,6), stat=stats )
      call umalloc(36,0,'cim1',memuse,stats)
      allocate( cip1(6,6), stat=stats )
      call umalloc(36,0,'cip1',memuse,stats)
      allocate( ckm1(6,6), stat=stats )
      call umalloc(36,0,'ckm1',memuse,stats)
      allocate( ckp1(6,6), stat=stats )
      call umalloc(36,0,'ckp1',memuse,stats)
      allocate( c(6,6), stat=stats )
      call umalloc(36,0,'c',memuse,stats)
      allocate( t(6,6), stat=stats )
      call umalloc(36,0,'t',memuse,stats)
      allocate( ul(6,6), stat=stats )
      call umalloc(36,0,'ul',memuse,stats)
      allocate( n1(20), stat=stats )
      call umalloc(20,1,'n1',memuse,stats)

      eps  = 100.*10.**(-iexp)
c 
c     n       = i,j,k
c     islavept(2) = j - 1
c     islavept(3) = j + 1
c     islavept(4) = k - 1
c     islavept(5) = k + 1
c     islavept(6) = i - 1
c     islavept(7) = i + 1
c     islavept(10) = number of ija and sa array values to reserve for this
c                    point
c     islavept(11) = max number of coincident points at block interfaces
c     islavept(12) = coincident pt 2 (pt n is the first)   
c     islavept(13) = coincident pt 3
c     islavept(14) = coincident pt 4
c     islavept(15) = coincident pt 5
c     islavept(16) = coincident pt 6
c     islavept(17) = coincident pt 7
c     islavept(18) = coincident pt 8
c     islavept(19) = coincident pt 9
c     islavept(20) = coincident pt 10
c     islavept(21) = number of the nearest surface pt
c
      ipt = 3*nslave+2 
      ija(1) = ipt
      do n = 1,nslave 
         ija(3*(n-1)+2)= ija(3*(n-1)+1) + islavept(n,10,iseqr)
         ija(3*(n-1)+3)= ija(3*(n-1)+2) + islavept(n,10,iseqr) 
         ija(3*n    +1)= ija(3*(n-1)+3) + islavept(n,10,iseqr) 
      enddo
      do n = 3*nslave+2,245*nslave+2
        ija(n) = 0
      enddo 
c
c
       if(nt.eq.1.and.meshdef.eq.1) then
          do n = 1,nslave
           iimax = islavept(n,11,iseqr)
           if(iimax.gt.1) then
                do ii2 = 2,iimax
                  n2 = islavept(n,12+ii2-2,iseqr)
                  test  = (xst(n2)-xst(n))*(xst(n2)-xst(n))
     .                   +(yst(n2)-yst(n))*(yst(n2)-yst(n))
     .                   +(zst(n2)-zst(n))*(zst(n2)-zst(n))
                  test = sqrt(test)
                  if(test.gt.100.*eps) then
                   write(11,29102) n,n2,test
29102              format(' WARNING: Macro-elements ',i6,' and ',i6,
     .           ' at 1-1 blocking interface have geometric mismatch'
     .              ,' = ',e16.8)
                  end if
                enddo
           end if
          enddo
       end if

c
      if(isktyp.gt.0) then 
        gini = 10.d0 
        eini = 10.d0 

        mxiter1 = 1      
c
        call elrhs(b,u,islavept,nslave,nmaster,iseqr) 

        do itr1 = 1,mxiter1 

          do n = 1,245*nslave+2
            sa(n) = 0.d0 
          enddo


          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21022) nslave 
21022        format(' Calculating macroelement moduli ',
     .               'for ',i8,' element nodes')
           end if
          end if
          call hookefe(ei,ej,ek,gij,gjk,gik,xst,yst,zst,volij,
     .                 volik,eps,eini,gini,arg1,arg2,islavept,
     .                 nslave,nmaster,nt,nnodes,myhost,myid,
     .                 mycomm,mblk2nd,maxbl,iseqr)

          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21023) 
21023        format(' Calculating macroelement metrics ')
           end if
          end if

          call elmetricsfe(xix,xiy,xiz,etax,etay,etaz,zetax,zetay,
     .                     zetaz,ooj,xst,yst,zst,eps,islavept,nslave,
     .                     nmaster,nnodes,myhost,myid,mycomm,
     .                     mblk2nd,maxbl,iseqr)
 
          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21024) 
21024        format(' Calculating and assembling macroelement stiffness'
     .                ,'  matrix ')
           end if
          end if

          call elglobfe(sa,xst,yst,zst,xix,xiy,xiz,etax,etay,
     .                  etaz,zetax,zetay,zetaz,ei,ej,ek,gij,gjk,
     .                  gik,ooj,eps,stiffl,islavept,ija,nslave,
     .                  nmaster,nnodes,myhost,myid,mycomm,mblk2nd,
     .                  nblelst,maxbl,iseqr)

          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21025) 
21025        format(' Solving the macroelement system ')
           end if
          end if
         
          ncount = 245*nslave+2
          ndim   = 3*nslave
c
          call gaussseidel(ncount,ndim,b,u,sa,ija,1e-6,9*nslave,err
     .                     ,volij,myid,myhost) 
c         call dprec(sa,ija,b,ndim,ncount)

c         call linbcg(ncount,ndim,b,u,sa,ija,1,1e-6,nslave,iter,err
c    .                ,time,eps,myid) 
          
          call coincdef(volij,u,islavept,nslave,nmaster,iseqr,n1)
        enddo
      else

          if(meshdef.eq.1) then
           if(myid.eq.myhost) then
             write(11,21027) 
21027        format(' Solving nodal displacements using exponential'
     .              ,' decay  ')
           end if
          end if
         
          call expdecay(xst,yst,zst,u,volij,ei,ripm,eps,
     .                  islavept,nslave,nmaster,ndim,nt,
     .                  nnodes,myhost,myid,mycomm,mblk2nd,
     .                  maxbl,iseqr)
          call coincdef(volij,u,islavept,nslave,nmaster,iseqr,n1)
      end if

        if(i2d.ne.0) then 
          do n = 1,nslave
            if(islavept(n,6,iseqr).eq.n) then
              nip1     = islavept(n,7,iseqr)
              u(3*(n-1)+1)    = .5*(u(3*(n-1)+1) + u(3*(nip1-1)+1))
              u(3*(nip1-1)+1) = u(3*(n-1)+1)
              u(3*(nip1-1)+2) = 0.
              u(3*(n-1)+2)    = 0. 
              u(3*(n-1)+3)    = .5*(u(3*(n-1)+3) + u(3*(nip1-1)+3))
              u(3*(nip1-1)+3) = u(3*(n-1)+3) 
            end if
          enddo
        end if
c
c     release memory
c
      deallocate(sa)
      deallocate(stiffl)
      deallocate(ija)
      deallocate(b)
      deallocate(ei)
      deallocate(ripm)
      deallocate(ej)
      deallocate(ek)
      deallocate(gij)
      deallocate(gjk)
      deallocate(gik)
      deallocate(xix)
      deallocate(xiy)
      deallocate(xiz)
      deallocate(etax)
      deallocate(etay)
      deallocate(etaz)
      deallocate(zetax)
      deallocate(zetay)
      deallocate(zetaz)
      deallocate(ooj)
      deallocate(volij)
      deallocate(volik)
      deallocate(tinv)
      deallocate(cjm1)
      deallocate(cjp1)
      deallocate(ckm1)
      deallocate(ckp1)
      deallocate(cim1)
      deallocate(cip1)
      deallocate(c)
      deallocate(t)
      deallocate(ul)
      deallocate(n1)


      return
      end
c
      subroutine coincdef(volij,u,islavept,nslave,nmaster,iseqr,n1)
      dimension islavept(nslave,nmaster,5)
      dimension u(3*nslave),volij(3*nslave)
      dimension n1(20) 
c
          volij = 0.
          do n = 1,nslave
           iimax = islavept(n,11,iseqr)
           if(iimax.gt.1) then
              n1(1) = n
              if(islavept(n,8,iseqr).ne.0) then
                ii3 = 0
                do ii2 = 2,iimax
                  n1(ii2) = islavept(n,12+ii2-2,iseqr)
                  if(islavept(n1(ii2),8,iseqr).eq.0) ii3 = ii2
                enddo
                if(ii3.eq.0) then
                  rooiim = 1./real(iimax)
                  do ii2 = 1,iimax
                    do j = 1,3
                      volij(3*(n-1)+j) = volij(3*(n-1)+j) + 
     .                                   u(3*(n1(ii2)-1)+j)*rooiim
                    enddo
                  enddo
                else
                  do j = 1,3
                    volij(3*(n-1)+j) = u(3*(n1(ii3)-1)+j)
                  enddo
                end if
              else
                 do j = 1,3
                   volij(3*(n-1)+j) = u(3*(n-1)+j)
                 enddo
              end if
           end if
          enddo
          do n = 1,nslave
           iimax = islavept(n,11,iseqr)
           if(iimax.gt.1) then
              n1(1) = n
              do ii2 = 2,iimax
                n1(ii2) = islavept(n,12+ii2-2,iseqr)
              enddo
              do ii2 = 1,iimax
                do j = 1,3
                  u(3*(n1(ii2)-1)+j) = volij(3*(n-1)+j)
                enddo
              enddo
           end if
          enddo

      return
      end
c
      SUBROUTINE gaussseidel(nc,n,b,x,sa,ija,tol,itmax,err, 
     .                       xt,myid,myhost) 
      parameter(EPS=1.d-12)  
      INTEGER iter,itmax,itol,n 
      INTEGER j 
      dimension b(n),x(n)
      dimension xt(n)
      dimension ija(nc)
      dimension sa(nc)
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt

      it1 = 0
      itmx1 =200
      err   = 1.e+10
      do it = 1,itmax
        sum = 0.
        resid = 0.
        do i = 1,n
          ut = b(i)/sa(i)  
          do j = ija(i),ija(i+1)-1
            ut = ut - sa(j)*x(ija(j))/sa(i)
          enddo
          xt(i) = alpha1*ut + (1.-alpha1)*x(i) 
          sum = sum + (xt(i)-x(i))*(xt(i)-x(i))
          resid1= sqrt((xt(i)-x(i))*(xt(i)-x(i)))
          if(resid1.gt.resid) then
           resid = resid1
           nresmx= i
          end if
          x(i) = xt(i)
        enddo
        errold = err
        err = sqrt(sum)
        nresmx = 1+(nresmx-1)/3
        if(meshdef.eq.1.and.myid.eq.myhost) 
     .               write(1000+myid,31029) it,err,nresmx,resid   
31029   format(i8,e16.8,i8,e16.8)
        if(err.lt.tol) goto 2000 
        if(errold.lt.err) then
          it1 = it1 + 1
          if(it1.gt.itmx1) then
             write(1000,21022) 
             stop
          end if
21022     format(' Stopping, Gauss-Seidel scheme not converging: '
     .        ,/,' Reduce BETA1 or change Macro-element definition')
        else
          it1 = 0
        end if
      enddo
2000  continue 
      if(myid.eq.myhost) write(1000+myid,31029) it,err   
 
      return
      end
      SUBROUTINE linbcg(nc,n,b,x,sa,ija,itol,tol,itmax,iter,err
     .                  ,time,eps,myid) 
      INTEGER iter,itmax,itol,n 
      INTEGER j 
      dimension b(n),x(n)
      dimension ija(nc)
      dimension sa(nc)
      dimension p(n),pp(n),r(n),rr(n),z(n),zz(n)
CU    USES atimes,asolve,snrm 
      iter=0 

      test = 0.
      do 5 j = 1,n
        if(test.lt.abs(b(j))) test = abs(b(j))
5     continue
      if(test.le.eps) return
      call atimes(n,nc,x,r,sa,ija,0,myid) 
      do 11 j=1,n 
        r(j)=b(j)-r(j) 
        rr(j)=r(j) 
11    continue 
      call atimes(n,nc,r,rr,sa,ija,0,myid) 
      znrm=1.d0 
      if(itol.eq.1) then 
        bnrm=snrm(n,b,itol) 
      else if (itol.eq.2) then 
        call asolve(n,nc,b,z,sa,0)
        bnrm=snrm(n,z,itol) 
      else if (itol.eq.3.or.itol.eq.4) then 
        call asolve(n,nc,b,z,sa,0)
        bnrm=snrm(n,z,itol) 
        call asolve(n,nc,r,z,sa,0)
        znrm=snrm(n,z,itol) 
      else 
        write(1000+myid,31922)
31922   format('illegal itol in linbcg') 
        stop
      endif 
      if(bnrm.lt.EPS) then
        err = 0.
        znrm= 0.
        goto 200 
      end if
      call asolve(n,nc,r,z,sa,0)
100   if (iter.le.itmax) then 
        iter=iter+1 
        call asolve(n,nc,rr,zz,sa,1)
        bknum=0.d0 
        do 12 j=1,n 
          bknum=bknum+z(j)*rr(j) 
12      continue 
        if(iter.eq.1) then 
          do 13 j=1,n 
            p(j)=z(j) 
            pp(j)=zz(j) 
13        continue 
        else 
          bk=bknum/bkden 
          do 14 j=1,n 
            p(j)=bk*p(j)+z(j) 
            pp(j)=bk*pp(j)+zz(j) 
14        continue 
        endif 
        bkden=bknum 
        call atimes(n,nc,p,z,sa,ija,0,myid) 
        akden=0.d0 
        do 15 j=1,n 
          akden=akden+z(j)*pp(j) 
15      continue 
        ak=bknum/akden 
        call atimes(n,nc,pp,zz,sa,ija,1,myid) 
c
c
        do 16 j=1,n 
          x(j)=x(j)+ak*p(j) 
          r(j)=r(j)-ak*z(j) 
          rr(j)=rr(j)-ak*zz(j) 
16      continue 
c
        
        call asolve(n,nc,r,z,sa,0)
        if(itol.eq.1.) then
          if(bnrm.gt.EPS) then
            err=snrm(n,r,itol)/bnrm 
          else
            err = 0.
          end if
             write(1000+myid,5040) iter,log10(err),znrm  
        else if(itol.eq.2) then 
          if(bnrm.gt.EPS) then
            err=snrm(n,z,itol)/bnrm 
          else
            err = 0.
          end if
             write(1000+myid,5040) iter,log10(err),znrm  
        else if(itol.eq.3.or.itol.eq.4)then 
          zm1nrm=znrm 
          znrm=snrm(n,z,itol)  
c
          if(abs(zm1nrm-znrm).gt.EPS*znrm) then 
            dxnrm=abs(ak)*snrm(n,p,itol) 
            if(abs(zm1nrm-znrm).gt.EPS) then
              err=znrm/abs(zm1nrm-znrm)*dxnrm 
            else
              err = 0.
            end if
          else
            if(bnrm.gt.EPS) then
              err=znrm/bnrm 
            else
              err = 0.
              znrm= 0. 
            end if
c            write(1000+myid,5042) iter,log10(err),znrm  
            goto 100 
          endif 
          xnrm=snrm(n,x,itol) 
c
          if(err.le.0.5d0*xnrm) then 
              if(xnrm.gt.EPS) then
                err=err/xnrm 
              else
                err = 0.
              end if
             write(1000+myid,5043) iter,log10(err),znrm  
          else 
            if(bnrm.gt.EPS) then
              err=znrm/bnrm 
            else
              err = 0.
              znrm= 0. 
            end if
c            write(1000+myid,5044) iter,log10(err),znrm  
            goto 100 
          endif 
        endif 
      if(err.gt.tol) goto 100 
      endif 
200   continue
5040    format('1',i8,2x,2(1x,e12.5)) 
5041    format('2',i8,2x,6(1x,e12.5)) 
5042    format('3',i8,2x,2(1x,e12.5)) 
5043    format('4',i8,2x,2(1x,e12.5)) 
5044    format('5',i8,2x,2(1x,e12.5)) 

      return 
      END 


C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3. 



      SUBROUTINE atimes(n,nc,x,r,sa,ija,itrnsp,myid) 
      INTEGER n,itrnsp
      dimension ija(nc)
      dimension x(n),r(n),sa(nc) 
CU    USES dsprsax,dsprstx 
      if (itrnsp.eq.0) then 
        call dsprsax(sa,ija,x,r,n,nc,myid) 
      else 
        call dsprstx(sa,ija,x,r,n,nc,myid) 
      endif 
      return 
      END 
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      SUBROUTINE dsprstx(sa,ija,x,b,n,nc,myid)
      INTEGER n,ija(nc)
      dimension b(n),sa(nc),x(n)
      INTEGER i,j,k
      if (ija(1).ne.n+2) then
        write(1000+myid,3111) ija(1),n+2,nc 
3111    format('mismatched vector and matrix in sprstx',3i8) 
        stop
      end if 
      do 11 i=1,n
        b(i)=sa(i)*x(i)
11    continue
      do 13 i=1,n
        do 12 k=ija(i),ija(i+1)-1
          j=ija(k)
          b(j)=b(j)+sa(k)*x(i)
12      continue
13    continue
c
      return
      END
c
c  Bartels preconditioner
c
      SUBROUTINE dprec(sa,ija,b,n,nc)
      INTEGER ija(nc)
      dimension sa(nc),b(n) 
      INTEGER i,k
      if (ija(1).ne.n+2) then
        write(1000+myid,3111) ija(1),n+2,nc 
3111    format('mismatched vector and matrix in sprsax',3i8) 
        stop
      end if 
      do 12 i=1,n
        b(i) = b(i)/sa(i)
        do 11 k=ija(i),ija(i+1)-1
          sa(k)=sa(k)/sa(i)
11      continue
        sa(i)=1.
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      SUBROUTINE dsprsax(sa,ija,x,b,n,nc,myid)
      INTEGER n,ija(nc)
      dimension b(n),sa(nc),x(n)
      INTEGER i,k
      if (ija(1).ne.n+2) then
        write(1000+myid,3111) ija(1),n+2,nc 
3111    format('mismatched vector and matrix in sprsax',3i8) 
        stop
      end if 
      do 12 i=1,n
        b(i)=sa(i)*x(i)
        do 11 k=ija(i),ija(i+1)-1
          b(i)=b(i)+sa(k)*x(ija(k))
11      continue
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      SUBROUTINE asolve(n,nc,b,x,sa,itrnsp) 
      INTEGER n,itrnsp,i 
      dimension x(n),b(n),sa(nc)
      do 11 i=1,n 
        x(i)=b(i)/sa(i) 
11    continue 

      return 
      END 
C  (C) Copr. 1986-92 Numerical Recipes Software ">u,3.

      FUNCTION snrm(n,sx,itol) 
      INTEGER n,itol,i,isamax 
      dimension sx(n)
      if (itol.le.3)then 
        snrm=0. 
        do 11 i=1,n 
          snrm=snrm+sx(i)**2 
11      continue 
        snrm=sqrt(snrm) 
      else 
        isamax=1 
        do 12 i=1,n 
          if(abs(sx(i)).gt.abs(sx(isamax))) isamax=i 
12      continue 
        snrm=abs(sx(isamax)) 
      endif 
      return 
      END 
c
************************************************************************
*
*         SUBROUTINE FOR THE INVERSION OF COMIN USING GAUSS-JORDAN
*         WITH PIVOTING.
*
************************************************************************
      SUBROUTINE INVDET(UL,COMIN,N,DTNRM,DETM)
      DIMENSION  UL(N,N),COMIN(N,N)
      DIMENSION C(N,N),J(550)
      DO 10 K = 1,n
        DO 12 I = 1,n
          C(I,K) = COMIN(I,K)
12      CONTINUE
10    CONTINUE
      PD = 1.D0
      DO 124 L = 1,N
        DD = 0.D0
        DO 123 K = 1,N
123     DD = DD + C(L,K)*C(L,K)
        DD = SQRT(DD)
124   PD = PD*DD
      DETM = 1.D0
      DO 125 L = 1,N
125   J(L+20) = L
      DO 144 L = 1,N
        CC = 0.D0
        M  = L
        DO 135 K = L,N
          IF((ABS(CC)-ABS(C(L,K))).GE.0.D0) GOTO 135
126       M = K
          CC = C(L,K)
135     CONTINUE
127     IF(L.EQ.M) GOTO 138
128     K = J(M+20)
        J(M+20) = J(L+20)
        J(L+20) = K
        DO 137 K = 1,N
          S = C(K,L)
          C(K,L) = C(K,M)
137     C(K,M) = S
138     C(L,L) = 1.D0
        DETM = DETM*CC
        DO 139 M = 1,N
139     C(L,M) = C(L,M)/CC
        DO 142 M = 1,N
          IF(L.EQ.M) GOTO 142
129       CC = C(M,L)
          IF(CC.EQ.0.D0) GOTO 142
130       C(M,L) = 0.D0
          DO 141 K = 1,N
141       C(M,K) = C(M,K) - CC*C(L,K)
142     CONTINUE
144   CONTINUE
      DO 143 L = 1,N
        IF(J(L+20).EQ.L) GOTO 143
131     M = L
132     M = M + 1
        IF(J(M+20).EQ.L) GOTO 133
136     IF(N.GT.M) GOTO 132
133     J(M+20) = J(L+20)
        DO 163 K = 1,N
          CC = C(L,K)
          C(L,K) = C(M,K)
163     C(M,K) = CC
        J(L+20) = L
143   CONTINUE
      DETM = ABS(DETM)
      DTNRM = DETM/PD
      DO 210 k = 1,N
        DO 212 I = 1,N
          UL(I,K) = C(I,K)
212     CONTINUE
210   CONTINUE
      RETURN
      END
      subroutine expdecay(xst,yst,zst,u,ut,r,ripm,eps,
     .                    islavept,nslave,nmaster,ndim,nt,
     .                    nnodes,myhost,myid,mycomm,mblk2nd,
     .                    maxbl,iseqr)
      dimension islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave)
      dimension ut(ndim),u(ndim)
      dimension r(nslave),ripm(7,nslave)
      dimension mblk2nd(maxbl)
      dimension arg(7),coef(8),rsn(7)
      dimension nipm(7),nsf(7),n3(7)
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt
      common /sgminmax/ rmin,rmax

      relax = 1.0
c
c    Do this at every time step because the closest surface
c    point to a control point will potentially change as a
c    surface moves.
c
        do n1 = 1,nslave
          r(n1) = 1.e6 
          if(islavept(n1,8,iseqr).ne.0) then
           do n2 = 1,nslave
            if(islavept(n2,8,iseqr).eq.0) then
              test  = (xst(n2)-xst(n1))*(xst(n2)-xst(n1))
     .               +(yst(n2)-yst(n1))*(yst(n2)-yst(n1))
     .               +(zst(n2)-zst(n1))*(zst(n2)-zst(n1))
              test  = sqrt(test)
              if(r(n1).gt.test) then
                r(n1)             = test
                islavept(n1,21,iseqr) = n2
              end if
            end if
           enddo
          else
           r(n1) = 0.
           islavept(n1,21,iseqr) = n1
          end if
        enddo
        if(nt.eq.1) then
          rmax = 0. 
          rmin = 1.e+30  
          do n = 1,nslave
            rave =  r(n)
            if(rave.gt.rmax) rmax = rave   
            if(rave.lt.rmin.and.rave.gt.eps) rmin = rave  
          enddo
        end if
        eps1 = rmin/10000.
        do n = 1,nslave
          nsf(1)  = islavept(n,8,iseqr)
          if(nsf(1).ne.0) then
            nipm(1)= n
            nipm(2)= islavept(n,2,iseqr)
            nipm(3)= islavept(n,3,iseqr)
            nipm(4)= islavept(n,4,iseqr)
            nipm(5)= islavept(n,5,iseqr)
            nipm(6)= islavept(n,6,iseqr)
            nipm(7)= islavept(n,7,iseqr)
            do ii = 1,7
              nsf(ii)  = islavept(nipm(ii),8,iseqr)
              if(n.ne.nipm(ii)) then 
                ripm(ii,n)=(xst(n)-xst(nipm(ii)))*(xst(n)-xst(nipm(ii)))
     .                    +(yst(n)-yst(nipm(ii)))*(yst(n)-yst(nipm(ii)))
     .                    +(zst(n)-zst(nipm(ii)))*(zst(n)-zst(nipm(ii)))
     .                     +eps1 
                ripm(ii,n) = sqrt(ripm(ii,n))
              else
                ripm(ii,n) = 1.e+30
              end if  
            enddo
            coef   = 0.
            coef(8)= 1.
            coef(1)= 1.
            do ii = 2,7
              if(n.ne.nipm(ii)) then
                coef(ii)= exp(-100.*ripm(ii,n)/rmax)
                coef(8)= coef(8) + coef(ii)
              end if
            enddo
            n3(1) = islavept(n,21,iseqr)
            rsn(1)= r(n)
            do ii = 2,7
              n3(ii) = islavept(nipm(ii),21,iseqr)
              rsn(ii)= r(nipm(ii))
            enddo
            n2    = n3(1)
            arg(1) = - beta2*(rsn(1)/rmax  - alpha2) 
            if(arg(1).gt.0.) arg(1) = 0.
            coef(1)= coef(1)/coef(8) 
            do ii = 2,7
              if(n.ne.nipm(ii)) then 
                n2    = n3(ii)
                arg(ii) = - beta2*(rsn(ii)/rmax  - alpha2) 
                if(arg(ii).gt.0.) arg(ii) = 0.
                coef(ii)= coef(ii)/coef(8) 
              else
                coef(ii) = 0.
                arg(ii)  = 0.
              end if
            enddo
            do j = 1,3
              n2 = n3(1)
              u(3*(n-1)+j) = u(3*(n2-1)+j)*exp(arg(1))
            enddo
          end if
        enddo
       do it = 1,nsprgit   
        do n = 1,nslave
          nipm(2)= islavept(n,2,iseqr)
          nipm(3)= islavept(n,3,iseqr)
          nipm(4)= islavept(n,4,iseqr)
          nipm(5)= islavept(n,5,iseqr)
          nipm(6)= islavept(n,6,iseqr)
          nipm(7)= islavept(n,7,iseqr)
          nsf(1) = islavept(n,8,iseqr)
c         do ii = 2,7
c           nsf(ii) = islavept(nipm(ii),8,iseqr)
c         enddo
          if(nsf(1).ne.0) then
            coef1 = 0. 
            do ii = 2,7
              coef1= coef1 + 1./ripm(ii,n)
            enddo
            do j = 1,3
              ut(3*(n-1)+j) = 0.
              do ii = 2,7
               ut(3*(n-1)+j)= ut(3*(n-1)+j)
     .             + u(3*(nipm(ii)-1)+j)/ripm(ii,n)/coef1   
              enddo
            enddo
          end if
        enddo
        do n = 1,nslave
          nipm(2)= islavept(n,2,iseqr)
          nipm(3)= islavept(n,3,iseqr)
          nipm(4)= islavept(n,4,iseqr)
          nipm(5)= islavept(n,5,iseqr)
          nipm(6)= islavept(n,6,iseqr)
          nipm(7)= islavept(n,7,iseqr)
          nsf(1) = islavept(n,8,iseqr)
c         do ii = 2,7
c           nsf(ii) = islavept(nipm(ii),8,iseqr)
c         enddo
          if(nsf(1).ne.0) then
            u(3*(n-1)+1) = relax*ut(3*(n-1)+1)+(1.-relax)*u(3*(n-1)+1) 
            u(3*(n-1)+2) = relax*ut(3*(n-1)+2)+(1.-relax)*u(3*(n-1)+2) 
            u(3*(n-1)+3) = relax*ut(3*(n-1)+3)+(1.-relax)*u(3*(n-1)+3) 
          end if
        enddo
       enddo

      return
      end

      subroutine hookefe(ei,ej,ek,gij,gjk,gik,xst,yst,zst,volij,
     .                   volik,eps,eini,gini,arg1,arg2,islavept,
     .                   nslave,nmaster,nt,nnodes,myhost,myid,
     .                   mycomm,mblk2nd,maxbl,iseqr)
      dimension islavept(nslave,nmaster,5)
      dimension xst(nslave),yst(nslave),zst(nslave)
      dimension ei(nslave),ej(nslave),ek(nslave),gij(nslave),
     .          gjk(nslave),gik(nslave)
      dimension volik(nslave),volij(3*nslave)
      dimension r(nslave),r1(0:8) 
      dimension ni1(0:8),ni0(0:8),test1(0:8) 
      dimension mblk2nd(maxbl)
      common /sgminmax/ rmin,rmax
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt

c
c    Do this at every time step because the closest surface
c    point to a control point will potentially change as a
c    surface moves.
c
        rmax = 0. 
        rmin = 1.e+30  
        do n1 = 1,nslave
          r(n1) = 1.e6 
          islavept(n1,21,iseqr) = 0
          ni0(0)    = n1
          ni0(1)    = islavept(n1    ,3,iseqr)
          ni0(2)    = islavept(n1    ,5,iseqr)
          ni0(3)    = islavept(ni0(1),5,iseqr)
          ni0(4)    = islavept(n1    ,7,iseqr)
          ni0(5)    = islavept(ni0(2),7,iseqr)
          ni0(6)    = islavept(ni0(1),7,iseqr)
          ni0(7)    = islavept(ni0(6),5,iseqr)
          if(ni0(1).ne.n1.and.ni0(2).ne.n1.and.ni0(4).ne.n1) then
           xcent = 0.
           ycent = 0.
           zcent = 0.
           do ii = 0,7
             xcent = xcent + xst(ni0(ii))/8. 
             ycent = ycent + yst(ni0(ii))/8. 
             zcent = zcent + zst(ni0(ii))/8. 
           enddo
           do n2 = 1,nslave
            if(islavept(n2,8,iseqr).eq.0) then
              ni1(0)    = n2
              ni1(1)    = islavept(n2    ,3,iseqr)
              ni1(2)    = islavept(n2    ,5,iseqr)
              ni1(3)    = islavept(ni1(1),5,iseqr)
              ni1(4)    = islavept(n2    ,7,iseqr)
              ni1(5)    = islavept(ni1(2),7,iseqr)
              ni1(6)    = islavept(ni1(1),7,iseqr)
              ni1(7)    = islavept(ni1(6),5,iseqr)
              test  = (xst(n2)-xcent)*(xst(n2)-xcent)
     .               +(yst(n2)-ycent)*(yst(n2)-ycent)
     .               +(zst(n2)-zcent)*(zst(n2)-zcent)
              test    = sqrt(test)
              if(r(n1).gt.test) then
                r(n1)                 = test   
                islavept(n1,21,iseqr) = n2
              end if
              if(islavept(ni1(1),8,iseqr).eq.0.and.
     .           islavept(ni1(2),8,iseqr).eq.0) then 
                 xsct = (xst(n2)+xst(ni1(1))+xst(ni1(2))
     .                    +xst(ni1(3)))/4.
                 ysct = (yst(n2)+yst(ni1(1))+yst(ni1(2))
     .                    +yst(ni1(3)))/4.
                 zsct = (zst(n2)+zst(ni1(1))+zst(ni1(2))
     .                    +zst(ni1(3)))/4.
                 test  = (xsct-xcent)*(xsct-xcent)
     .                  +(ysct-ycent)*(ysct-ycent)
     .                  +(zsct-zcent)*(zsct-zcent)
                 test    = sqrt(test)
                 if(r(n1).gt.test) then
                    r(n1)                 = test   
                    islavept(n1,21,iseqr) = n2
                 end if
              end if
              if(islavept(ni1(2),8,iseqr).eq.0.and.
     .           islavept(ni1(4),8,iseqr).eq.0) then 
                 xsct = (xst(n2)+xst(ni1(2))+xst(ni1(4))
     .                    +xst(ni1(5)))/4.
                 ysct = (yst(n2)+yst(ni1(2))+yst(ni1(4))
     .                    +yst(ni1(5)))/4.
                 zsct = (zst(n2)+zst(ni1(2))+zst(ni1(4))
     .                    +zst(ni1(5)))/4.
                 test  = (xsct-xcent)*(xsct-xcent)
     .                  +(ysct-ycent)*(ysct-ycent)
     .                  +(zsct-zcent)*(zsct-zcent)
                 test    = sqrt(test)
                 if(r(n1).gt.test) then
                    r(n1)                 = test   
                    islavept(n1,21,iseqr) = n2
                 end if
              end if
              if(islavept(ni1(1),8,iseqr).eq.0.and.
     .           islavept(ni1(4),8,iseqr).eq.0) then 
                 xsct = (xst(n2)+xst(ni1(1))+xst(ni1(4))
     .                    +xst(ni1(6)))/4.
                 ysct = (yst(n2)+yst(ni1(1))+yst(ni1(4))
     .                    +yst(ni1(6)))/4.
                 zsct = (zst(n2)+zst(ni1(1))+zst(ni1(4))
     .                    +zst(ni1(6)))/4.
                 test  = (xsct-xcent)*(xsct-xcent)
     .                  +(ysct-ycent)*(ysct-ycent)
     .                  +(zsct-zcent)*(zsct-zcent)
                 test    = sqrt(test)
                 if(r(n1).gt.test) then
                    r(n1)                 = test   
                    islavept(n1,21,iseqr) = n2
                 end if
              end if
            end if
           enddo
          end if 
          if(islavept(n1,21,iseqr).eq.0) then
            r(n1) = 0.
            islavept(n1,21,iseqr) = n1
          end if
        enddo
        rmax = 0.
        rmin = 1.e+30
        do n1 = 1,nslave
          n2    = islavept(n1,21,iseqr)
          if(r(n1).gt.rmax) rmax = r(n1)
          if(r(n1).lt.rmin.and.r(n1).gt.eps) rmin = r(n1)
        enddo
      e0 = 500000.
      r0 = (rmax*alog(1.-eini/e0)+beta1*rmin)/
     .     (     alog(1.-eini/e0)+beta1     ) 
  
      do n = 1,nslave
       nbl = islavept(n,9,iseqr)
        ni0(1)    = islavept(n     ,3,iseqr)
        ni0(2)    = islavept(n     ,5,iseqr)
        ni0(4)    = islavept(n     ,7,iseqr)
        if(ni0(1).ne.n.and.ni0(2).ne.n.and.ni0(4).ne.n) then
          rave = r(n)
          denom    = (1.-exp(-beta1*(rave-r0)/(rmax-r0))) 
          if(denom.lt.eps) denom = eps
          fact  = 1./denom
          ej(n) = eini*fact
          ek(n) = eini*fact
          ei(n) = eini*fact
          gij(n)= gini*fact
          gik(n)= gini*fact
          gjk(n)= gini*fact
        end if
      enddo
      return
      end

      subroutine elmetricsfe(xix,xiy,xiz,etax,etay,etaz,zetax,zetay,
     .                       zetaz,ooj,xs,ys,zs,eps,islavept,nslave,
     .                       nmaster,nnodes,myhost,myid,mycomm,
     .                       mblk2nd,maxbl,iseqr)

      dimension islavept(nslave,nmaster,5)
      dimension xs(nslave),ys(nslave),zs(nslave)
      dimension xix(nslave),xiy(nslave),xiz(nslave),
     .          etax(nslave),etay(nslave),etaz(nslave),
     .          zetax(nslave),zetay(nslave),zetaz(nslave),
     .          ooj(nslave)
      dimension mblk2nd(maxbl)

c
      do n = 1,nslave
       nbl = islavept(n,9,iseqr)
c      if (myid .eq. mblk2nd(nbl)) then 
        njp    = islavept(n,3,iseqr)
        nkp    = islavept(n,5,iseqr)
        njpkp  = islavept(nkp,3,iseqr)
        nip    = islavept(n,7,iseqr)
        nipkp  = islavept(nkp,7,iseqr)
        njpip  = islavept(njp,7,iseqr)
        njpipkp= islavept(njpip,5,iseqr)
        if(njp.ne.n.and.nip.ne.n.and.nkp.ne.n) then
          xxi    =.25*(xs(njp)  -xs(n)  +xs(njpkp)  -xs(nkp)
     .                +xs(njpip)-xs(nip)+xs(njpipkp)-xs(nipkp))
          xeta   =.25*(xs(nkp)  -xs(n)  +xs(njpkp)  -xs(njp)
     .                +xs(nipkp)-xs(nip)+xs(njpipkp)-xs(njpip))
          xzeta  =.25*(xs(nip)  -xs(n)  +xs(nipkp)  -xs(nkp)
     .                +xs(njpip)-xs(njp)+xs(njpipkp)-xs(njpkp))
          yxi    =.25*(ys(njp)  -ys(n)  +ys(njpkp)  -ys(nkp)
     .                +ys(njpip)-ys(nip)+ys(njpipkp)-ys(nipkp))
          yeta   =.25*(ys(nkp)  -ys(n)  +ys(njpkp)  -ys(njp)
     .                +ys(nipkp)-ys(nip)+ys(njpipkp)-ys(njpip))
          yzeta  =.25*(ys(nip)  -ys(n)  +ys(nipkp)  -ys(nkp)
     .                +ys(njpip)-ys(njp)+ys(njpipkp)-ys(njpkp))
          zxi    =.25*(zs(njp)  -zs(n)  +zs(njpkp)  -zs(nkp)
     .                +zs(njpip)-zs(nip)+zs(njpipkp)-zs(nipkp))
          zeta   =.25*(zs(nkp)  -zs(n)  +zs(njpkp)  -zs(njp)
     .                +zs(nipkp)-zs(nip)+zs(njpipkp)-zs(njpip))
          zzeta  =.25*(zs(nip)  -zs(n)  +zs(nipkp)  -zs(nkp)
     .                +zs(njpip)-zs(njp)+zs(njpipkp)-zs(njpkp))
          ooj(n)  = 1./(  xxi*(yeta*zzeta-yzeta*zeta)
     .                  -xeta*(yxi *zzeta-yzeta*zxi )
     .                 +xzeta*(yxi *zeta -yeta *zxi )) 
          xix(n)  =(yeta*zzeta-yzeta*zeta)*ooj(n)
          xiy(n)  =(xzeta*zeta-xeta*zzeta)*ooj(n)
          xiz(n)  =(xeta*yzeta-xzeta*yeta)*ooj(n)
          etax(n) =(yzeta*zxi -yxi *zzeta)*ooj(n)
          etay(n) =(xxi*zzeta -xzeta* zxi)*ooj(n)
          etaz(n) =(xzeta*yxi -xxi *yzeta)*ooj(n)
          zetax(n)=(yxi *zeta -yeta*  zxi)*ooj(n)
          zetay(n)=(xeta*zxi  -xxi * zeta)*ooj(n)
          zetaz(n)=(xxi *yeta -xeta*  yxi)*ooj(n)
        else
          ooj(n)  = 0.
          xix(n)  = 0.
          xiy(n)  = 0.
          xiz(n)  = 0.
          etax(n) = 0.
          etay(n) = 0.
          etaz(n) = 0.
          zetax(n)= 0.
          zetay(n)= 0.
          zetaz(n)= 0.
        end if
c      end if
      enddo
  
      return
      end
c
      subroutine elrhs(b,u,islavept,nslave,nmaster,iseqr) 
      dimension islavept(nslave,nmaster,5)
      dimension b(3*nslave)
      dimension u(3*nslave)
c
c
      do n = 1,nslave
         ii4  = 1
         ni3  = 0 
         if(islavept(n,8,iseqr).eq.0) ii4 = 0
         if(ii4.ne.0) then
           iimax = islavept(n,11,iseqr)
           do ii2 = 2,iimax
             ni3 = islavept(n,12+ii2-2,iseqr)
             if(islavept(ni3,8,iseqr).eq.0) then
               ii4 = 0
               goto 1500
             end if
           enddo
         end if
1500     continue

       if(ii4.eq.0.and.ni3.eq.0) then
         b(3*(n-1)+1) =  u(3*(n-1)+1)
         b(3*(n-1)+2) =  u(3*(n-1)+2)
         b(3*(n-1)+3) =  u(3*(n-1)+3)
       else if(ii4.eq.0.and.ni3.ne.0) then
         b(3*(n-1)+1) =  u(3*(ni3-1)+1)
         b(3*(n-1)+2) =  u(3*(ni3-1)+2)
         b(3*(n-1)+3) =  u(3*(ni3-1)+3)
       else
         b(3*(n-1)+1) =  0.
         b(3*(n-1)+2) =  0.
         b(3*(n-1)+3) =  0.
       end if
      enddo

      return
      end
      subroutine deform(nbl,idim,jdim,kdim,x,y,z,xnm2,ynm2,znm2,
     .                  xnm1,ynm1,znm1,deltj,deltk,delti,u,
     .                  vel,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                  maxbl,dt,nou,bou,nbuf,ibufdim,myid,
     .                  idefrm,nbci0,nbcidim,nbcj0,nbcjdim,nbck0,
     .                  nbckdim,ibcinfo,jbcinfo,kbcinfo,maxseg,wk,nsurf,
     .                  irst,iflag,islavept,nslave,iskip,jskip,
     .                  kskip,nsegdfrm,idfrmseg,iaesurf,maxsegdg,
     .                  nmaster,iseq,iskmax,jskmax,kskmax,nt)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute mesh deformation via Transfinite Interpolation
c              using arc-length blending functions, compute grid speeds
c              at grid points, store of the current grid as the one at
c              the previous time step, and update the grid shape.
c
c     irst    = 0 standard update; do all steps described above
c               1 restart - compute the mesh deformation and update
c                 only the grid shape; do not update the grid speeds.
c
c     isktyp  = 0 use only transfinite interpolation 
c             > 0 do isktyp smoothing steps using the spring analogy 
c                 and mesh shape preservation near the k=1 and k=kdim
c                 surfaces
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      integer stats
c
      dimension delti(jdim,kdim,3,2)
      dimension deltj(kdim,idim,3,2)
      dimension deltk(jdim,idim,3,2)
      dimension iaesurf(maxbl,maxsegdg)
      dimension ibcinfo(maxbl,maxseg,7,2)
      dimension icsf(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg)
      dimension idefrm(maxbl)
      dimension idfrmseg(maxbl,maxsegdg)
      dimension iskip(maxbl,500)
      dimension jbcinfo(maxbl,maxseg,7,2)
      dimension jcsf(maxbl,maxsegdg)
      dimension jcsi(maxbl,maxsegdg)
      dimension jskip(maxbl,500)
      dimension kbcinfo(maxbl,maxseg,7,2)
      dimension kcsf(maxbl,maxsegdg)
      dimension kcsi(maxbl,maxsegdg)
      dimension kskip(maxbl,500)
      dimension nbci0(maxbl)
      dimension nbcidim(maxbl)
      dimension nbcj0(maxbl)
      dimension nbcjdim(maxbl)
      dimension nbck0(maxbl)
      dimension nbckdim(maxbl)
      dimension nou(nbuf)
      dimension nsegdfrm(maxbl)
      dimension islavept(nslave,nmaster,5)
      dimension vel(jdim,kdim,idim,3)
      dimension wk(9*nsurf)
      dimension x(jdim*kdim*idim)
      dimension xnm1(jdim*kdim*idim)
      dimension xnm2(jdim*kdim*idim)
      dimension y(jdim*kdim*idim)
      dimension ynm1(jdim*kdim*idim)
      dimension ynm2(jdim*kdim*idim)
      dimension z(jdim*kdim*idim)
      dimension znm1(jdim*kdim*idim)
      dimension znm2(jdim*kdim*idim)
      dimension u(3*nslave)
      dimension iskmax(maxbl)
      dimension jskmax(maxbl)
      dimension kskmax(maxbl)
c
      allocatable :: arci(:,:,:)
      allocatable :: arcj(:,:,:)
      allocatable :: arck(:,:,:)
      allocatable :: dx(:)
      allocatable :: dx1(:,:,:)
      allocatable :: dx2(:,:,:)
      allocatable :: dx3(:,:,:)
      allocatable :: dx4(:,:,:)
      allocatable :: dx5(:,:,:)
      allocatable :: dx6(:,:,:)
      allocatable :: dx7(:,:,:)
      allocatable :: dy(:)
      allocatable :: dy1(:,:,:)
      allocatable :: dy2(:,:,:)
      allocatable :: dy3(:,:,:)
      allocatable :: dy4(:,:,:)
      allocatable :: dy5(:,:,:)
      allocatable :: dy6(:,:,:)
      allocatable :: dy7(:,:,:)
      allocatable :: dz(:)
      allocatable :: dz1(:,:,:)
      allocatable :: dz2(:,:,:)
      allocatable :: dz3(:,:,:)
      allocatable :: dz4(:,:,:)
      allocatable :: dz5(:,:,:)
      allocatable :: dz6(:,:,:)
      allocatable :: dz7(:,:,:)
      allocatable :: ibl(:)

      common /twod/ i2d
      common /deformz/ beta1,beta2,alpha1,alpha2,isktyp,negvol,meshdef,
     .                 nsprgit,ndgrd,ndwrt
      common /unst/ time,cfltau,ntstep,ita,iunst,cfltau0,cfltauMax 
      common /moov/movie,nframes,icall1,lhdr,icoarsemovie,i2dmovie 
c
c     allocate memory
c
      memuse = 0
c
      allocate( arci(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'arci',memuse,stats)
      allocate( arcj(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'arcj',memuse,stats)
      allocate( arck(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'arck',memuse,stats)
      allocate( dx(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx',memuse,stats)
      allocate( dx1(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx1',memuse,stats)
      allocate( dx2(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx2',memuse,stats)
      allocate( dx3(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx3',memuse,stats)
      allocate( dx4(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx4',memuse,stats)
      allocate( dx5(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx5',memuse,stats)
      allocate( dx6(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx6',memuse,stats)
      allocate( dx7(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx7',memuse,stats)
      allocate( dy(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy',memuse,stats)
      allocate( dy1(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy1',memuse,stats)
      allocate( dy2(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy2',memuse,stats)
      allocate( dy3(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy3',memuse,stats)
      allocate( dy4(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy4',memuse,stats)
      allocate( dy5(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy5',memuse,stats)
      allocate( dy6(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy6',memuse,stats)
      allocate( dy7(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy7',memuse,stats)
      allocate( dz(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz',memuse,stats)
      allocate( dz1(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz1',memuse,stats)
      allocate( dz2(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz2',memuse,stats)
      allocate( dz3(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz3',memuse,stats)
      allocate( dz4(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz4',memuse,stats)
      allocate( dz5(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz5',memuse,stats)
      allocate( dz6(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz6',memuse,stats)
      allocate( dz7(jdim,kdim,idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz7',memuse,stats)
      allocate( ibl(nsurf), stat=stats )
      call umalloc(nsurf,1,'ibl',memuse,stats)
c
c     determine arc lengths for use in TFI blending functions
c
      call arclen(idim,jdim,kdim,arci,arcj,arck,xnm1,ynm1,
     .            znm1,nbl,nou,bou,nbuf,ibufdim,myid)
c
c     initialize deltas
c
      do j=1,idim*jdim*kdim
        dx(j) = 0.
        dy(j) = 0.
        dz(j) = 0.
      end do
      do n = 1,nslave
        nbl1= islavept(n,9,iseq)
        ll  = islavept(n,1,iseq)
        if(nbl1.eq.nbl) then
          dx(ll+1) = u(3*(n-1)+1)
          dy(ll+1) = u(3*(n-1)+2)
          dz(ll+1) = u(3*(n-1)+3)
        end if
      enddo
      if(meshdef.eq.1) then
        inc = ntstep
        if(movie.ne.0) then
          rinc1 = real(nt)/real(movie)
          inc1  = nt/movie
          if(real(inc1).eq.rinc1) inc = nt
        end if
        if(nt.eq.inc) then
          write(4000+myid,81621) jskmax(nbl),kskmax(nbl),iskmax(nbl),nbl
     .                        ,nt
81621 format('ZONE I = ',i8,'  J = ',i8,' K = ',i8,
     .               ' T="Block ',i5,' Time Step = ',i6,'"')
31920   format(3(1x,e16.8),1x,e16.8,2(1x,e16.8),2i8)
          do n = 1,nslave
            nbl1= islavept(n,9,iseq)
            ll  = islavept(n,1,iseq)
            n2  = islavept(n,21,iseq)
            if(nbl1.eq.nbl) then
              write(4000+myid,31920) x(ll+1)+u(3*(n-1)+1)
     .                             ,y(ll+1)+u(3*(n-1)+2)
     .                             ,z(ll+1)+u(3*(n-1)+3)
     .                             ,u(3*(n-1)+1)
     .                             ,u(3*(n-1)+2)
     .                             ,u(3*(n-1)+3),n,n2
            end if
          enddo
        end if
      end if
c
c
      if(abs(isktyp).eq.1) then
        iskp = iskip(nbl,1)
        jskp = jskip(nbl,1)
        kskp = kskip(nbl,1)
      end if
c
c     i=1 to idim  subfaces
c
      if(abs(isktyp).eq.1) then
       do i = 1,idim,iskp     
         js = 1
         je = jdim
         ks = 1
         ke = kdim
         do j=js,je-jskp,jskp
            do k=ks,ke,kskp
               call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                      j,j+jskp,k,k,arci,arcj,arck,
     .                      nou,bou,nbuf,ibufdim,myid,nbl)
            end do
         end do
         do j=js,je,jskp
            do k=ks,ke-kskp,kskp
               call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                      j,j,k,k+kskp,arci,arcj,arck,
     .                      nou,bou,nbuf,ibufdim,myid,nbl)
            end do
         end do
       enddo
      else
       do i1 = 1,iskmax(nbl)     
         i = iskip(nbl,i1)
         do j1=1,jskmax(nbl)-1 
            j = jskip(nbl,j1)
            jskp = jskip(nbl,j1+1)
            do k1=1,kskmax(nbl) 
               k = kskip(nbl,k1)
               call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                      j,jskp,k,k,arci,arcj,arck,
     .                      nou,bou,nbuf,ibufdim,myid,nbl)
            end do
         end do
         do j1=1,jskmax(nbl)
            j = jskip(nbl,j1)
            do k1=1,kskmax(nbl)-1
               k = kskip(nbl,k1)
               kskp = kskip(nbl,k1+1)
               call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i,
     .                      j,j,k,kskp,arci,arcj,arck,
     .                      nou,bou,nbuf,ibufdim,myid,nbl)
            end do
         end do
       enddo
      end if
c
c            preserve deltas on edges of solid surfaces
c
             if (idefrm(nbl).lt.999) then
                  call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                         jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                         idim,maxbl,maxsegdg,nsegdfrm)
             end if
c
c     i=1 subfaces
c
      do nseg=1,nbci0(nbl)
         ii     = 1
         ibctyp = ibcinfo(nbl,nseg,1,ii)
         if (abs(ibctyp).ne.2004 .and.
     .       abs(ibctyp).ne.2014 .and.
     .       abs(ibctyp).ne.2024 .and.
     .       abs(ibctyp).ne.2034 .and.
     .       abs(ibctyp).ne.2016 .and.
     .       abs(ibctyp).ne.1005 .and.
     .       abs(ibctyp).ne.1006) then 
             js = ibcinfo(nbl,nseg,2,ii)
             je = ibcinfo(nbl,nseg,3,ii)
             ks = ibcinfo(nbl,nseg,4,ii)
             ke = ibcinfo(nbl,nseg,5,ii)
             if(abs(isktyp).eq.1) then
               do j=js,je-jskp,jskp
                  do k=ks,ke-kskp,kskp
                     call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                            dx2,dy2,dz2,arci,arcj,arck,1,1,
     .                            j,j+jskp,k,k+kskp,nou,bou,nbuf,
     .                            ibufdim,myid)
                  end do
               end do
             else
               do j1=1,jskmax(nbl)-1
                  j = jskip(nbl,j1)
                  jskp = jskip(nbl,j1+1)
                  do k1=1,kskmax(nbl)-1
                     k = kskip(nbl,k1)
                     kskp = kskip(nbl,k1+1)
                     if((j.ge.js.and.j.lt.je).and.(k.ge.ks.
     .                  and.k.lt.ke)) then
                       call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                              dx2,dy2,dz2,arci,arcj,arck,1,1,
     .                              j,jskp,k,kskp,nou,bou,nbuf,
     .                              ibufdim,myid)
                     end if
                  enddo
               enddo
             end if
         end if
      enddo
c
c     i=idim subfaces
c
      do nseg=1,nbcidim(nbl)
         ii     = 2
         ibctyp = ibcinfo(nbl,nseg,1,ii)
         if (abs(ibctyp).ne.2004 .and.
     .       abs(ibctyp).ne.2014 .and.
     .       abs(ibctyp).ne.2024 .and.
     .       abs(ibctyp).ne.2034 .and.
     .       abs(ibctyp).ne.2016 .and.
     .       abs(ibctyp).ne.1005 .and.
     .       abs(ibctyp).ne.1006) then
             js = ibcinfo(nbl,nseg,2,ii)
             je = ibcinfo(nbl,nseg,3,ii)
             ks = ibcinfo(nbl,nseg,4,ii)
             ke = ibcinfo(nbl,nseg,5,ii)
             if(abs(isktyp).eq.1) then
               do j=js,je-jskp,jskp
                  do k=ks,ke-kskp,kskp
                     call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                            dx2,dy2,dz2,arci,arcj,arck,idim,idim,
     .                            j,j+jskp,k,k+kskp,nou,bou,nbuf,
     .                            ibufdim,myid)
                  end do
               end do
             else
               do j1=1,jskmax(nbl)-1
                  j = jskip(nbl,j1)
                  jskp = jskip(nbl,j1+1)
                  do k1=1,kskmax(nbl)-1
                     k = kskip(nbl,k1)
                     kskp = kskip(nbl,k1+1)
                     if((j.ge.js.and.j.lt.je).and.(k.ge.ks.
     .                  and.k.lt.ke)) then
                       call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                              dx2,dy2,dz2,arci,arcj,arck,idim,     
     .                              idim,j,jskp,k,kskp,nou,bou,nbuf,
     .                              ibufdim,myid)
                     end if
                  enddo
               enddo
             end if
         end if
      enddo
c
c     Intermediate i subfaces
c
      if(idim.gt.2) then
        if(abs(isktyp).eq.1) then 
         do i = 1+iskp,idim-iskp       
          js = 1
          je = jdim
          ks = 1
          ke = kdim
          do j=js,je-jskp,jskp
            do k=ks,ke-kskp,kskp
               call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                      dx2,dy2,dz2,arci,arcj,arck,i,i,
     .                      j,j+jskp,k,k+kskp,nou,bou,nbuf,
     .                      ibufdim,myid)
            end do
          end do
         enddo
        else
         do i1 = 2,iskmax(nbl)-1     
          i = iskip(nbl,i1)
          do j1=1,jskmax(nbl)-1 
            j = jskip(nbl,j1)
            jskp = jskip(nbl,j1+1)
            do k1=1,kskmax(nbl)-1 
               k = kskip(nbl,k1)
               kskp = kskip(nbl,k1+1)
               call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                      dx2,dy2,dz2,arci,arcj,arck,i,i,
     .                      j,jskp,k,kskp,nou,bou,nbuf,
     .                      ibufdim,myid)
            enddo
          enddo
         enddo
        end if
      end if
c
      if (i2d .eq. 0) then
c
c     j=1 to jdim  subfaces
c
        if(abs(isktyp).eq.1) then
         do j = 1,jdim,jskp     
           is = 1
           ie = idim
           ks = 1
           ke = kdim
           do i=is,ie-iskp,iskp
             do k=ks,ke,kskp
               call tfiedge(idim,jdim,kdim,dx,dy,dz,i,i+iskp,
     .                      j,j,k,k,arci,arcj,arck,
     .                      nou,bou,nbuf,ibufdim,myid,nbl)
             end do
           end do
         enddo
        else
         do j1 = 1,jskmax(nbl)     
           j = jskip(nbl,j1)
           do i1=1,iskmax(nbl)-1 
             i = iskip(nbl,i1)
             iskp = iskip(nbl,i1+1)
             do k1=1,kskmax(nbl) 
               k = kskip(nbl,k1)
               call tfiedge(idim,jdim,kdim,dx,dy,dz,i,iskp,
     .                      j,j,k,k,arci,arcj,arck,
     .                      nou,bou,nbuf,ibufdim,myid,nbl)
             end do
           end do
         enddo
        end if
c
c            preserve deltas on edges of solid surfaces
c
        if (idefrm(nbl).lt.999) then
           call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                  jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                  idim,maxbl,maxsegdg,nsegdfrm)
        end if
c
c        j=1 subfaces
c
         do nseg=1,nbcj0(nbl)
           jj     = 1
           jbctyp = jbcinfo(nbl,nseg,1,jj)
           if (abs(jbctyp).ne.2004 .and.
     .         abs(jbctyp).ne.2014 .and.
     .         abs(jbctyp).ne.2024 .and.
     .         abs(jbctyp).ne.2034 .and.
     .         abs(jbctyp).ne.2016 .and.
     .         abs(jbctyp).ne.1005 .and.
     .         abs(jbctyp).ne.1006) then
               is = jbcinfo(nbl,nseg,2,jj)
               ie = jbcinfo(nbl,nseg,3,jj)
               ks = jbcinfo(nbl,nseg,4,jj)
               ke = jbcinfo(nbl,nseg,5,jj)
               if(abs(isktyp).eq.1) then
                do i=is,ie-iskp,iskp
                   do k=ks,ke-kskp,kskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             1,1,k,k+kskp,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                end do
               else
                do i1=1,iskmax(nbl)-1
                  i = iskip(nbl,i1)
                  iskp = iskip(nbl,i1+1)
                  do k1=1,kskmax(nbl)-1
                     k = kskip(nbl,k1)
                     kskp = kskip(nbl,k1+1)
                     if((i.ge.is.and.i.lt.ie).and.(k.ge.ks.
     .                  and.k.lt.ke)) then
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,iskp,
     .                             1,1,k,kskp,nou,bou,nbuf,
     .                             ibufdim,myid)
                     end if
                  enddo
                enddo
               end if
           end if
         end do
c
c        j=jdim subfaces
c
         do nseg=1,nbcjdim(nbl)
           jj     = 2
           jbctyp = jbcinfo(nbl,nseg,1,jj)
           if (abs(jbctyp).ne.2004 .and.
     .         abs(jbctyp).ne.2014 .and.
     .         abs(jbctyp).ne.2024 .and.
     .         abs(jbctyp).ne.2034 .and.
     .         abs(jbctyp).ne.2016 .and.
     .         abs(jbctyp).ne.1005 .and.
     .         abs(jbctyp).ne.1006) then
               is = jbcinfo(nbl,nseg,2,jj)
               ie = jbcinfo(nbl,nseg,3,jj)
               ks = jbcinfo(nbl,nseg,4,jj)
               ke = jbcinfo(nbl,nseg,5,jj)
               if(abs(isktyp).eq.1) then
                do i=is,ie-iskp,iskp
                   do k=ks,ke-kskp,kskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             jdim,jdim,k,k+kskp,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                end do
               else
                do i1=1,iskmax(nbl)-1
                  i = iskip(nbl,i1)
                  iskp = iskip(nbl,i1+1)
                  do k1=1,kskmax(nbl)-1
                     k = kskip(nbl,k1)
                     kskp = kskip(nbl,k1+1)
                     if((i.ge.is.and.i.lt.ie).and.(k.ge.ks.
     .                  and.k.lt.ke)) then
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,iskp,
     .                             jdim,jdim,k,kskp,nou,bou,nbuf,
     .                             ibufdim,myid)
                     end if
                  enddo
                enddo
               end if
            end if
         end do
c
c     Intermediate j subfaces
c
         if(abs(isktyp).eq.1) then
           do j = 1+jskp,jdim-jskp       
             is = 1
             ie = idim
             ks = 1
             ke = kdim
             do i=is,ie-iskp,iskp
                do k=ks,ke-kskp,kskp
                   call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                          dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                          j,j,k,k+kskp,nou,bou,nbuf,
     .                          ibufdim,myid)
                end do
             end do
           end do
         else
           do j1 = 2,jskmax(nbl)-1     
             j = jskip(nbl,j1)
             do i1=1,iskmax(nbl)-1
               i = iskip(nbl,i1)
               iskp = iskip(nbl,i1+1)
               do k1=1,kskmax(nbl)-1
                  k = kskip(nbl,k1)
                  kskp = kskip(nbl,k1+1)
                  call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                         dx2,dy2,dz2,arci,arcj,arck,i,iskp,
     .                         j,j,k,kskp,nou,bou,nbuf,
     .                         ibufdim,myid)
               enddo
             enddo
           enddo
         end if
c
c        k=1 subfaces
c
         do nseg=1,nbck0(nbl)
           kk     = 1
           kbctyp = kbcinfo(nbl,nseg,1,kk)
           if (abs(kbctyp).ne.2004 .and.
     .         abs(kbctyp).ne.2014 .and.
     .         abs(kbctyp).ne.2024 .and.
     .         abs(kbctyp).ne.2034 .and.
     .         abs(kbctyp).ne.2016 .and.
     .         abs(kbctyp).ne.1005 .and.
     .         abs(kbctyp).ne.1006) then
               is = kbcinfo(nbl,nseg,2,kk)
               ie = kbcinfo(nbl,nseg,3,kk)
               js = kbcinfo(nbl,nseg,4,kk)
               je = kbcinfo(nbl,nseg,5,kk)
               if(abs(isktyp).eq.1) then
                 do i=is,ie-iskp,iskp
                   do j=js,je-jskp,jskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             j,j+jskp,1,1,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                 end do
               else
                do i1=1,iskmax(nbl)-1
                  i = iskip(nbl,i1)
                  iskp = iskip(nbl,i1+1)
                  do j1=1,jskmax(nbl)-1
                     j = jskip(nbl,j1)
                     jskp = jskip(nbl,j1+1)
                     if((i.ge.is.and.i.lt.ie).and.(j.ge.js.
     .                  and.j.lt.je)) then
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,iskp,
     .                             j,jskp,1,1,nou,bou,nbuf,
     .                             ibufdim,myid)
                     end if
                  enddo
                enddo
               end if
           end if
         end do
c
c        k=kdim subfaces
c
         do nseg=1,nbckdim(nbl)
           kk     = 2
           kbctyp = kbcinfo(nbl,nseg,1,kk)
           if (abs(kbctyp).ne.2004 .and.
     .         abs(kbctyp).ne.2014 .and.
     .         abs(kbctyp).ne.2024 .and.
     .         abs(kbctyp).ne.2034 .and.
     .         abs(kbctyp).ne.2016 .and.
     .         abs(kbctyp).ne.1005 .and.
     .         abs(kbctyp).ne.1006) then
               is = kbcinfo(nbl,nseg,2,kk)
               ie = kbcinfo(nbl,nseg,3,kk)
               js = kbcinfo(nbl,nseg,4,kk)
               je = kbcinfo(nbl,nseg,5,kk)
               if(abs(isktyp).eq.1) then
                 do i=is,ie-iskp,iskp
                   do j=js,je-jskp,jskp
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                             j,j+jskp,kdim,kdim,nou,bou,nbuf,
     .                             ibufdim,myid)
                   end do
                 end do
               else
                do i1=1,iskmax(nbl)-1
                  i = iskip(nbl,i1)
                  iskp = iskip(nbl,i1+1)
                  do j1=1,jskmax(nbl)-1
                     j = jskip(nbl,j1)
                     jskp = jskip(nbl,j1+1)
                     if((i.ge.is.and.i.lt.ie).and.(j.ge.js.
     .                  and.j.lt.je)) then
                      call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                             dx2,dy2,dz2,arci,arcj,arck,i,iskp,
     .                             j,jskp,kdim,kdim,nou,bou,nbuf,
     .                             ibufdim,myid)
                     end if
                  enddo
                enddo
               end if
           end if
         end do
c
c     Intermediate k subfaces
c
         if(abs(isktyp).eq.1) then
           do k = 1+kskp,kdim-kskp       
             is = 1
             ie = idim
             js = 1
             je = jdim
             do i=is,ie-iskp,iskp
               do j=js,je-jskp,jskp
                  call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                         dx2,dy2,dz2,arci,arcj,arck,i,i+iskp,
     .                         j,j+jskp,k,k,nou,bou,nbuf,
     .                         ibufdim,myid)
               end do
             end do
           end do
         else
           do k1=2,kskmax(nbl)-1
             k = kskip(nbl,k1)
             do i1=1,iskmax(nbl)-1
               i = iskip(nbl,i1)
               iskp = iskip(nbl,i1+1)
               do j1=1,jskmax(nbl)-1
                  j = jskip(nbl,j1)
                  jskp = jskip(nbl,j1+1)
                   call tfiface(idim,jdim,kdim,dx,dy,dz,dx1,dy1,dz1,
     .                          dx2,dy2,dz2,arci,arcj,arck,i,iskp,
     .                          j,jskp,k,k,nou,bou,nbuf,
     .                          ibufdim,myid)
               enddo
             enddo
           enddo
         end if
c
c     TFI to get deltas throughout the volume
c
         call tfivol(idim,jdim,kdim,iskp,jskp,kskp,iskmax,jskmax,
     .               kskmax,iskip,jskip,kskip,isktyp,dx,dy,dz,
     .               dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3,dx4,dy4,
     .               dz4,dx5,dy5,dz5,dx6,dy6,dz6,dx7,dy7,dz7,
     .               arci,arcj,arck,nou,bou,nbuf,ibufdim,myid,
     .               maxbl,nbl)
c
      end if
c
      if (irst.eq.0) then
c
c        compute grid velocites before overwriting current grid
c
         if (real(time).le.real(dt) .or. abs(ita).eq.1) then
c
c           first order backward derivatives
c
            fact = 1./dt
            do i=1,idim
               do k=1,kdim
                  do j=1,jdim
                     jj = j + (k-1)*jdim + (i-1)*jdim*kdim
                     vel(j,k,i,1) = fact*((dx(jj)+x(jj)) - xnm1(jj))
                     vel(j,k,i,2) = fact*((dy(jj)+y(jj)) - ynm1(jj))
                     vel(j,k,i,3) = fact*((dz(jj)+z(jj)) - znm1(jj))
                  end do
               end do
            end do
c
         else
c
c           second order backward derivatives
c
            fact = 1./(2.*dt)
            do i=1,idim
               do k=1,kdim
                  do j=1,jdim
                     jj = j + (k-1)*jdim + (i-1)*jdim*kdim
                     vel(j,k,i,1) = fact*(3.*(dx(jj)+x(jj))
     .                            - 4*xnm1(jj)+xnm2(jj))
                     vel(j,k,i,2) = fact*(3.*(dy(jj)+y(jj))
     .                            - 4*ynm1(jj)+ynm2(jj))
                     vel(j,k,i,3) = fact*(3.*(dz(jj)+z(jj))
     .                            - 4*znm1(jj)+znm2(jj))
                  end do
               end do
            end do
c
         end if
c
      end if
c
c     for second order case, store current--->old before
c     updating grid
c
      if (abs(ita).gt.1) then
         do i=1,jdim*kdim*idim
           xnm2(i) = xnm1(i)
           ynm2(i) = ynm1(i)
           znm2(i) = znm1(i)
         end do
      end if
c
c     add deltas to the current mesh to get the new one
c
      do i=1,jdim*kdim*idim
        x(i) = x(i) + dx(i)
        y(i) = y(i) + dy(i)
        z(i) = z(i) + dz(i)
      end do
c
c     release memory
c
      deallocate(dx1)
      deallocate(dy1)
      deallocate(dz1)
      deallocate(dx2)
      deallocate(dy2)
      deallocate(dz2)
      deallocate(dx3)
      deallocate(dy3)
      deallocate(dz3)
      deallocate(dx4)
      deallocate(dy4)
      deallocate(dz4)
      deallocate(dx5)
      deallocate(dy5)
      deallocate(dz5)
      deallocate(dx6)
      deallocate(dy6)
      deallocate(dz6)
      deallocate(dx7)
      deallocate(dy7)
      deallocate(dz7)
      deallocate(dx)
      deallocate(dy)
      deallocate(dz)
      deallocate(arci)
      deallocate(arcj)
      deallocate(arck)
      deallocate(ibl)
c
      return
      end
      subroutine deform_surf(nbl,idim,jdim,kdim,deltj,deltk,delti,
     .                       lw,lw2,icsi,icsf,jcsi,jcsf,kcsi,kcsf,
     .                       maxbl,mseq,time,dt,ita,nou,bou,nbuf,
     .                       ibufdim,myid,idefrm,nbci0,nbcidim,nbcj0,
     .                       nbcjdim,nbck0,nbckdim,ibcinfo,jbcinfo,
     .                       kbcinfo,maxseg,wk,u,nsurf,irst,iflag,
     .                       islavept,nslave,nsegdfrm,idfrmseg,iaesurf,
     .                       maxsegdg,nmaster,iseq)
c
c     $Id$
c
c***********************************************************************
c     Purpose: Compute solid surface deformations for initializing the
c              mesh deformation scheme.
c
c***********************************************************************
c
#   ifdef CMPLX
      implicit complex(a-h,o-z)
#   endif
c
      character*120 bou(ibufdim,nbuf)
c
      integer stats
c
      dimension lw(65,maxbl),lw2(43,maxbl)
      dimension delti(jdim,kdim,3,2)
      dimension deltj(kdim,idim,3,2)
      dimension deltk(jdim,idim,3,2)
      dimension iaesurf(maxbl,maxsegdg)
      dimension ibcinfo(maxbl,maxseg,7,2)
      dimension icsf(maxbl,maxsegdg)
      dimension icsi(maxbl,maxsegdg)
      dimension idefrm(maxbl)
      dimension idfrmseg(maxbl,maxsegdg)
      dimension jbcinfo(maxbl,maxseg,7,2)
      dimension jcsf(maxbl,maxsegdg)
      dimension jcsi(maxbl,maxsegdg)
      dimension kbcinfo(maxbl,maxseg,7,2)
      dimension kcsf(maxbl,maxsegdg)
      dimension kcsi(maxbl,maxsegdg)
      dimension nbci0(maxbl)
      dimension nbcidim(maxbl)
      dimension nbcj0(maxbl)
      dimension nbcjdim(maxbl)
      dimension nbck0(maxbl)
      dimension nbckdim(maxbl)
      dimension nou(nbuf)
      dimension nsegdfrm(maxbl)
      dimension islavept(nslave,nmaster,5)
      dimension wk(9*nsurf)
      dimension u(3*nslave)
c
      allocatable :: dx(:)
      allocatable :: dy(:)
      allocatable :: dz(:)

      common /twod/ i2d
c
c     allocate memory
c
      memuse = 0
c
      allocate( dx(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dx',memuse,stats)
      allocate( dy(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dy',memuse,stats)
      allocate( dz(jdim*kdim*idim), stat=stats )
      call umalloc(jdim*kdim*idim,0,'dz',memuse,stats)
c
c     initialize deltas
c
      do i=1,jdim*kdim*idim
        dx(i) = 0.
        dy(i) = 0.
        dz(i) = 0.
      end do
c
c               preserve deltas on edges of solid surfaces
c
       if (idefrm(nbl).lt.999) then
          call bc_delt(nbl,dx,dy,dz,deltj,deltk,delti,jcsi,
     .                 jcsf,kcsi,kcsf,icsi,icsf,jdim,kdim,
     .                 idim,maxbl,maxsegdg,nsegdfrm)
c
        do n = 1,nslave
          nbl1= islavept(n,9,iseq)
          ll  = islavept(n,1,iseq)
          isrf= islavept(n,8,iseq)
          if(nbl1.eq.nbl.and.isrf.eq.0) then
            u(3*(n-1)+1) = dx(ll+1)
            u(3*(n-1)+2) = dy(ll+1)
            u(3*(n-1)+3) = dz(ll+1)
          end if
        enddo
      end if
c
c     release memory
c
      deallocate(dx)
      deallocate(dy)
      deallocate(dz)


      return
      end
