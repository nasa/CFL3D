      PROGRAM getuplus_flat
c
c  This is an OLD program
c  There are no guarantees as to its usability or correctness
c  Use at your own risk
c
      PARAMETER(II=451,JJ=199,KK=JJ-1)
      DIMENSION UVAL(15),YVAL1(15),VVAL(16),YVAL2(16),U1(JJ),Y1(JJ) 
      DIMENSION V1(JJ),U2(JJ),Y2(JJ),V2(JJ),XX(2),T1(JJ),T2(JJ) 
      DIMENSION X(II,KK),Y(II,KK),U(II,KK),V(II,KK),RHO(II,KK),P(II,KK) 
      DIMENSION TITLE(8),TVAL(16),YVAL3(16),PATTERN(2),TATTERN(4) 
      DIMENSION UBOY(40),YBOY(40),UGIRL(40),YGIRL(40),UVALM5(26), 
     +YVALM5(26),USPAL(52),YSPAL(52)
      character filnam*40
C 
C   PROGRAM TO READ FLAT PLATE DATA AND PLOT IT VS BLASIUS SOLUTION 
C   NOTE THAT U AND V VELOCITIES ARE NOT CONTRAVARIANT (ALONG AND NORMAL
C   TO THE BODY), BUT IN THE X AND Y DIRECTIONS 
C 
      DATA UVAL/0.,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.99,1.,-.90,.2/ 
      DATA YVAL1/0.,.06,.12,.182,.246,.307,.376,.45,.545,.687,.789,1.,
     +1.6,-.3,.2/ 
      DATA VVAL/0.,.02,.05,.1,.2,.3,.4,.5,.6,.7,.75,.8,.83,.85, 
     +-.90,.2/
      DATA YVAL2/0.,.1,.154,.221,.317,.394,.466,.541,.624,.727,.785,
     +.896,.992,1.1,-.3,.2/ 
      DATA TVAL/0.,0.,0.,0.,0.,0.,0., 
     +0.,0.,0.,0.,0.,0.,0.,.25,.1/
      DATA YVAL3/0.,0.,0.,0.,0.,0.,0.,
     +0.,0.,0.,0.,0.,0.,0.,-.3,.2/
      DATA UVALM5/0.,.0945,.189,.282,.374,.464,.549,.628,.7,.765, 
     +.82,.866,.904,.933,.955,.971,.982,.989,.993,.996,.998,.999, 
     +.999,1.0,-.90,.2/ 
      DATA YVALM5/0.,.0551,.110,.165,.220,.275,.329,.384,.438,.492, 
     +.546,.6,.653,.707,.760,.813,.866,.919,.972,1.02,1.08,1.13,
     +1.18,1.24,-.3,.2/ 
C
      filnam = 'plot3dg.bin'
      open(unit=2,file=filnam,form='unformatted',status='old') 
      filnam = 'plot3dq.bin'
      open(unit=3,file=filnam,form='unformatted',status='old')
      open(unit=17,file='flatplate.dat',form='formatted',
     +  status='unknown') 
      RHOE=1. 
      TINF=460. 
      C2B=198.6/TINF
      C2BP=C2B+1. 
      GAMMA=1.4 
      PATTERN(1)=.25
      PATTERN(2)=-.25 
      TATTERN(1)=.20
      TATTERN(2)=-.20 
      TATTERN(3)=.5 
      TATTERN(4)=-.20 
C 
c   read grid
      read(2) nbl
      if (nbl .ne. 1) then
        write(6,'('' currently can only read 1 zone!'')')
        write(6,'(''    there are '',i5)') nbl
        stop
      end if
      read(2) ni,nj
      if(ni .gt. ii .or. nj .gt. jj) then
        write(6,'('' error, must increase ii and jj to:'',2i5)') ni,nj
        stop
      end if
      read(2) ((x(i,j),i=1,ni),j=1,nj),
     +        ((y(i,j),i=1,ni),j=1,nj)
c   read q file (conserved)
      read(3) nbl
      if (nbl .ne. 1) then
        write(6,'('' currently can only read 1 zone!'')')
        write(6,'(''    there are '',i5)') nbl
        stop
      end if
      read(3) nix,njx
      if (nix .ne. ni .or. njx .ne. nj) then
        write(6,'('' error, indices dont match grid file'')')
        stop
      end if
      read(3) tmach,alp,reuew,time
      read(3)  ((rho(i,j),i=1,ni),j=1,nj),
     +         ((  u(i,j),i=1,ni),j=1,nj),
     +         ((  v(i,j),i=1,ni),j=1,nj),
     +         ((  p(i,j),i=1,ni),j=1,nj)
c   convert to primitive (p goes to p/pinf)
      do i=1,ni
      do j=1,nj
        u(i,j)=u(i,j)/rho(i,j)
        v(i,j)=v(i,j)/rho(i,j)
        p(i,j)=(gamma*(gamma-1.))*(p(i,j)-0.5*rho(i,j)*(u(i,j)**2+
     +         v(i,j)**2))
      enddo
      enddo

c 
      IQUIT  = 1
      IF(IQUIT .EQ. 0) GOTO 8888
C 
      REX = 6000000.
      XIND1 = .3
      XIND2 = .8
      IND1=1
      IND2=1
      XMR1=ABS(X(1,1)-XIND1)
      XMR2=ABS(X(1,1)-XIND2)
      DO 4022 I=2,NI
      IF(ABS(X(I,1)-XIND1) .LT. XMR1) THEN
        IND1=I
        XMR1=ABS(X(I,1)-XIND1)
      END IF
      IF(ABS(X(I,1)-XIND2) .LT. XMR2) THEN
        IND2=I
        XMR2=ABS(X(I,1)-XIND2)
      END IF
 4022 CONTINUE
      WRITE(6,'('' I-INDICES ARE'',2I6,'' FOR X ='',2E15.8)') IND1, 
     +IND2,X(IND1,1),X(IND2,1)
      CCORD = 1.
C 
      IUD = 1
 111  DELDET=1.
      DELDET = 0. 
      IF(IUD.EQ.0 .AND. DELDET.EQ.0.) GOTO 111
      IF(IUD.EQ.1 .AND. DELDET.EQ.0) THEN 
      ISPAL  = 1
      END IF
      IEXP = 0
C   THE 2 SETS OF EXP DATA MUST BE IN TAPE4 
C   FIRST LINE: NO OF PTS IN THE FIRST SET (ZERO = NONE)
C   NEXT: FIRST SET OF DATA, U AND Y....
C   REPEAT FOR SECOND SET OF DATA 
      IF(IEXP.EQ.1) THEN
      READ(4,*) NPT1
      IF(NPT1.GT.0) THEN
      DO 862 J=1,NPT1 
 862  READ(4,*) UBOY(J),YBOY(J) 
      END IF
      READ(4,*) NPT2
      IF(NPT2.GT.0) THEN
      DO 863 J=1,NPT2 
 863  READ(4,*) UGIRL(J),YGIRL(J) 
      END IF
      ELSE
      NPT1=0
      NPT2=0
      END IF
C 
      IF(CCORD.LE.0.) THEN
      WRITE(6,'('' DO YOU WANT BLASIUS ALSO (1=YES)?'')') 
      READ(5,*) IBLAS 
      ELSE
      IBLAS=0 
      END IF
C 
      NJJ=NJ
      IF(NJ .GT. (JJ-2)) NJJ=JJ-2 
C   RE WILL EQUAL UE*X/VNU
C 
      DO 300 I=1,2
      IF(I .EQ. 1) K=IND1 
      IF(I .EQ. 2) K=IND2 
      YNN=Y(K+1,1)-Y(K-1,1) 
      XNN=X(K+1,1)-X(K-1,1) 
      YNORM=SQRT(XNN*XNN+YNN*YNN) 
      XN=XNN/YNORM
      YN=YNN/YNORM
C 
      IF(IUD.EQ.1) THEN 
      UE=TMACH
      else
      HINF=(1.+.5*(GAMMA-1.)*TMACH*TMACH)/(GAMMA-1.)
      RHOEE=P(K,1)**(1./GAMMA)
      UETEMP=((GAMMA-1.)*RHOEE*HINF-P(K,1))/ 
     +       (.5*(GAMMA-1.)*RHOEE) 
      if (uetemp .le. 0.) then
        write(6,'('' Error finding UE.  Do you wish to scale with'',
     + '' UINF instead?'')')
        write(6,'('' (1=yes, 0=no,stop program):'')')
        read(5,*) ippol
        if(ippol .eq. 1) then
          UE=TMACH
        else
          stop
        end if
      else
        ue=sqrt(uetemp)
      end if
      ZZ=UE/TMACH 
      WRITE(6,'('' UE/UINF = '',E15.8)') ZZ 
      END IF
C 
      XX(I)=X(K,1)
      VNU=UE/REX
      RE=SQRT(UE*XX(I)/VNU) 
      IF(CCORD.LE.0.) THEN
      D99=5.*XX(I)/RE 
      ELSE
      D99=CCORD 
      END IF
      IF(I .EQ. 1) THEN 
         QSQR=P(K,1)/RHO(K,1) 
         FMU1=C2BP*(QSQR**1.5)/(C2B+QSQR) 
         DO 100 J=1,NJJ 
         U1(J)=(U(K,J)*XN+V(K,J)*YN)/UE 
         Y1(J)=(SQRT((Y(K,J)-Y(K,1))**2+(X(K,J)-X(K,1))**2))/D99
         V1(J)=(-U(K,J)*YN+V(K,J)*XN)*RE/UE 
         T1(J)=P(K,J)/RHO(K,J)
 100  CONTINUE
         CF1=FMU1*(U1(2)-U1(1))/(0.5*Y1(2)*RHOE*REX)
      ELSE
         QSQR=P(K,1)/RHO(K,1) 
         FMU2=C2BP*(QSQR**1.5)/(C2B+QSQR) 
         DO 200 J=1,NJJ 
         U2(J)=(U(K,J)*XN+V(K,J)*YN)/UE 
         Y2(J)=(SQRT((Y(K,J)-Y(K,1))**2+(X(K,J)-X(K,1))**2))/D99
         V2(J)=(-U(K,J)*YN+V(K,J)*XN)*RE/UE 
         T2(J)=P(K,J)/RHO(K,J)
 200  CONTINUE
         CF2=FMU2*(U2(2)-U2(1))/(0.5*Y2(2)*RHOE*REX)
      END IF
 300  CONTINUE
      IF(DELDET .EQ. 0.) THEN 
      NJJ=NJJ-1 
      DO 8671 J=1,NJJ 
c   ***rumsey
c     U1(J)=U1(J+1)*SQRT(2./CF1)
c     U2(J)=U2(J+1)*SQRT(2./CF2)
      U1(J)=U1(J+1)*SQRT(2.*rho(ind1,1)/CF1)
      U2(J)=U2(J+1)*SQRT(2.*rho(ind2,1)/CF2)
c   *********
      V1(J)=V1(J+1) 
      V2(J)=V2(J+1) 
      T1(J)=T1(J+1) 
      T2(J)=T2(J+1) 
c   ***rumsey
c     Y1(J)=ALOG10(Y1(J+1)*CCORD*(SQRT(CF1/2.))*REX)
c8671 Y2(J)=ALOG10(Y2(J+1)*CCORD*(SQRT(CF2/2.))*REX)
      Y1(J)=ALOG10(Y1(J+1)*CCORD*(SQRT(CF1*rho(ind1,1)/
     + (2.*fmu1**2)))*REX)
      Y2(J)=ALOG10(Y2(J+1)*CCORD*(SQRT(CF2*rho(ind2,1)/
     + (2.*fmu2**2)))*REX)
 8671 continue
c   *********
      END IF
C 
C   U-VELOCITY PROFILES 
C 
      NBL=13
      IF(IUD.EQ.1) THEN 
      IF(DELDET .GT. 0) THEN
c     CALL AXES(1.5,1.5,0.,12.,-.6,.2,1.,0.,'U/$(U INF$)',
c    + .28,-11,2) 
      ELSE
c     CALL AXES(1.5,1.5,90.,8.,-5.,5.,1.,0.,'U+',.28,2,2) 
      END IF
      ELSE
c     CALL AXES(1.5,1.5,0.,12.,-.6,.2,1.,0.,'U/$(U EDGE$)',
c    + .28,-12,2) 
      END IF
      IF(CCORD .LE. 0.) THEN
c     CALL AXES(1.5,1.5,90.,8.,0.,.2,1.,0.,'Y/DELTA',.28,7,2) 
      ELSE
      IF(DELDET .GT. 0) THEN
c     CALL AXES(1.5,1.5,90.,8.,0.,DELDET,1.,0.,'Y/CHORD',.28,7,3) 
      ELSE
c     CALL AXESLOG(1.5,1.5,0.,12.,-1.,2.,1,0,'Y+',.28,-2) 
      END IF
      END IF
      IF(DELDET .GT. 0.) THEN 
      U1(NJJ+1)=-.6-(1.5*.2)
      U1(NJJ+2)=.2
      U2(NJJ+1)=-.6-(1.5*.2)
      U2(NJJ+2)=.2
      ELSE
      U1(NJJ+1)=-5.-(1.5*5.)
      U1(NJJ+2)=5.
      U2(NJJ+1)=U1(NJJ+1) 
      U2(NJJ+2)=U1(NJJ+2) 
      END IF
      IF(CCORD.LE.0.) THEN
      Y1(NJJ+1)=-.3 
      Y1(NJJ+2)=.2
      Y2(NJJ+1)=-.3 
      Y2(NJJ+2)=.2
      ELSE
      IF(DELDET .GT. 0.) THEN 
      Y1(NJJ+1)=-DELDET*1.5 
      Y1(NJJ+2)=DELDET
      Y2(NJJ+1)=-DELDET*1.5 
      Y2(NJJ+2)=DELDET
      ELSE
      Y1(NJJ+1)=-1.-(1.5*.5)
      Y1(NJJ+2)=.5
      Y2(NJJ+1)=Y1(NJJ+1) 
      Y2(NJJ+2)=Y1(NJJ+2) 
      END IF
      END IF
      IF(NPT1.NE.0) THEN
      UBOY(NPT1+1)=U1(NJJ+1)
      UBOY(NPT1+2)=U1(NJJ+2)
      YBOY(NPT1+1)=Y1(NJJ+1)
      YBOY(NPT1+2)=Y1(NJJ+2)
      END IF
      IF(NPT2.NE.0) THEN
      UGIRL(NPT2+1)=U1(NJJ+1) 
      UGIRL(NPT2+2)=U1(NJJ+2) 
      YGIRL(NPT2+1)=Y1(NJJ+1) 
      YGIRL(NPT2+2)=Y1(NJJ+2) 
      END IF
      IF(IBLAS .EQ. 1) THEN 
c     CALL DASHPLT(UVAL,YVAL1,NBL,1,PATTERN,2)
c     CALL CHARACT(2.,6.5,.2,'DASH = BLASIUS',0.,14,0)
      END IF
      IF(NPT1.NE.0) THEN
      IF(DELDET .GT. 0.) THEN 
c     CALL LINPLT(UBOY,YBOY,NPT1,1,-1,1,1,1)
      ELSE
c     CALL LINPLT(YBOY,UBOY,NPT1,1,-1,1,1,1)
      END IF
      END IF
      IF(NPT2.NE.0) THEN
      IF(DELDET .GT. 0.) THEN 
c     CALL LINPLT(UGIRL,YGIRL,NPT2,1,-1,2,1,1)
      ELSE
c     CALL LINPLT(YGIRL,UGIRL,NPT2,1,-1,2,1,1)
      END IF
      END IF
      IF(DELDET.NE.0.) THEN 
c     CALL LINPLT(U1,Y1,NJJ,1,+1,2,1,1)
c     CALL LINPLT(U2,Y2,NJJ,1,+1,3,1,1)
      ELSE
c     CALL LINPLT(Y1,U1,NJJ,1,-1,2,1,1)
c     CALL LINPLT(Y2,U2,NJJ,1,-1,3,1,1)
      IF(ISPAL .EQ. 1) THEN 
      XPL=0.5 
      DO 6231 J=1,50
      XPL=XPL+0.5 
      USPAL(J)=XPL
      YSPAL(J)=XPL+EXP(-.4*5.5)*(EXP(0.4*XPL)-1.-0.4*XPL- 
     + (((0.4*XPL)**2)/2.)-(((0.4*XPL)**3)/6.)) 
      YSPAL(J)=ALOG10(YSPAL(J)) 
 6231 CONTINUE
      USPAL(51)=-5.-(1.5*5.)
      USPAL(52)=5.
      YSPAL(51)=-1.-(1.5*.5)
      YSPAL(52)=.5
c      CALL DASHPLT(YSPAL,USPAL,50,1,PATTERN,2)
c     call linplt(yspal,uspal,50,1,0,0,3,1)
c     CALL CHARACT(2.,6.5,.2,'LINE   = SPALDING',0.,17,0) 
      END IF
      END IF
c     CALL CHARACT(2.,6.1,.2,'SQUARE = X AT ',0.,14,0)
c     CALL NUMBER(4.5,6.1,.2,XX(1),0.,3)
c     CALL CHARACT(2.,5.7,.2,'DIAMOND= X AT ',0.,14,0)
c     CALL NUMBER(4.5,5.7,.2,XX(2),0.,3)
      if(deldet .ne. 0.) then
      none=1
c     WRITE(7,'(2I6,'' u vs y at x='',f12.5)') NJJ,none,xx(1) 
      write(17,'(''variables="u","y"'')')
      write(17,'(''zone,t="x='',f12.5,''"'')') xx(1)
      DO 500 J=1,NJJ
c     WRITE(7,'(2E15.6)') U1(J),Y1(J) 
      WRITE(17,'(2E15.6)') U1(J),Y1(J) 
 500  continue
c     WRITE(7,'(2I6,'' u vs y at x='',f12.5)') NJJ,none,xx(2) 
      write(17,'(''zone,t="x='',f12.5,''"'')') xx(2)
      DO 501 J=1,NJJ
c     WRITE(7,'(2E15.6)') U2(J),Y2(J) 
      WRITE(17,'(2E15.6)') U2(J),Y2(J) 
 501  continue
      IF(IBLAS .EQ. 1) THEN 
c       WRITE(7,'(I6,'' Balsius soln u vs y'')') NBL 
        write(17,'(''zone,t="Blasius"'')')
        DO 502 J=1,NBL
c       WRITE(7,'(2E15.6)') UVAL(J),YVAL1(J)
        WRITE(17,'(2E15.6)') UVAL(J),YVAL1(J)
 502    continue
      END IF
      else
c     WRITE(7,'(I6,'' log(y+) vs u+ at x='',f12.5)') NJJ,xx(1) 
      write(17,'(''variables="log(y+)","u+"'')')
      write(17,'(''zone,t="x='',f12.5,''"'')') xx(1)
      DO 507 J=1,NJJ
c     WRITE(7,'(2E15.6)') y1(J),u1(J) 
      WRITE(17,'(2E15.6)') y1(J),u1(J) 
 507  continue
c     WRITE(7,'(I6,'' log(y+) vs u+ at x='',f12.5)') NJJ,xx(2) 
      write(17,'(''zone,t="x='',f12.5,''"'')') xx(2)
      DO 508 J=1,NJJ
c     WRITE(7,'(2E15.6)') y2(J),u2(J) 
      WRITE(17,'(2E15.6)') y2(J),u2(J) 
 508  continue
      if(ispal .eq. 1) then
        nnnh=50
c       write(7,'(i6,'' Spalding theory log(y+) vs u+'')') nnnh
        write(17,'(''zone,t="Spalding theory"'')')
        do 509 j=1,nnnh
c       write(7,'(2e15.6)') yspal(j),uspal(j)
        write(17,'(2e15.6)') yspal(j),uspal(j)
 509    continue
      end if
      end if
c     CALL NFRAME 
C 
C   TEMPERATURE PROFILES
C 
      IF(IQUIT .EQ. 1) GOTO 8888
      NBL=14
c     CALL AXES(1.5,1.5,0.,12.,.4,.1,1.,0.,'TEMP/(TEMP INF)',.28,-15,2) 
      IF(CCORD .LE. 0.) THEN
c     CALL AXES(1.5,1.5,90.,8.,0.,.2,1.,0.,'Y/DELTA',.28,7,2) 
      ELSE
      IF(DELDET .GT. 0) THEN
c     CALL AXES(1.5,1.5,90.,8.,0.,DELDET,1.,0.,'Y/CHORD',.28,7,3) 
      ELSE
c     CALL AXESLOG(1.5,1.5,90.,8.,-1.,1.,1,0,'Y+',.28,2)
      END IF
      END IF
      T1(NJJ+1)=0.4-(1.5*0.1) 
      T1(NJJ+2)=.1
      T2(NJJ+1)=0.4-(1.5*0.1) 
      T2(NJJ+2)=.1
      IF(IBLAS .EQ. 1) THEN 
c     CALL DASHPLT(TVAL,YVAL3,NBL,1,PATTERN,2)
c     CALL CHARACT(2.,6.5,.2,'DASH = BLASIUS',0.,14,0)
      END IF
c     CALL LINPLT(T1,Y1,NJJ,1,0,0,1,1)
c     CALL DASHPLT(T2,Y2,NJJ,1,TATTERN,4) 
c     CALL CHARACT(2.,6.1,.2,'SOLID  = X AT ',0.,14,0)
c     CALL NUMBER(4.5,6.1,.2,XX(1),0.,3)
c     CALL CHARACT(2.,5.7,.2,'DTDASH = X AT ',0.,14,0)
c     CALL NUMBER(4.5,5.7,.2,XX(2),0.,3)
c     WRITE(9,'(I6)') NJJ 
c     DO 700 J=1,NJJ
c700  WRITE(9,'(2E15.8)') T1(J),Y1(J) 
c     WRITE(9,'(I6)') NJJ 
c     DO 701 J=1,NJJ
c701  WRITE(9,'(2E15.8)') T2(J),Y2(J) 
      IF(IBLAS .EQ. 1) THEN 
c     WRITE(9,'(I6)') NBL 
c     DO 702 J=1,NBL
c702  WRITE(9,'(2E15.8)') TVAL(J),YVAL3(J)
      END IF
c     CALL NFRAME 
C 
C   V-VELOCITY PROFILES 
C 
      IF(IQUIT .EQ. 2) GOTO 8888
      NBL=14
      IF(IUD.EQ.1) THEN 
c     CALL AXES(1.5,1.5,0.,12.,-.6,.2,1.,0.,'(V/U INF)*SQRT(RE)', 
c    +.28,-18,2)
      ELSE
c     CALL AXES(1.5,1.5,0.,12.,-.6,.2,1.,0.,'(V/U EGDE)*SQRT(RE)',
c    +.28,-19,2)
      END IF
      IF(CCORD .LE. 0.) THEN
c     CALL AXES(1.5,1.5,90.,8.,0.,.2,1.,0.,'Y/DELTA',.28,7,2) 
      ELSE
      IF(DELDET .GT. 0) THEN
c     CALL AXES(1.5,1.5,90.,8.,0.,DELDET,1.,0.,'Y/CHORD',.28,7,3) 
      ELSE
c     CALL AXESLOG(1.5,1.5,90.,8.,-1.,1.,1,0,'Y+',.28,2)
      END IF
      END IF
      V1(NJJ+1)=-.6-(1.5*.2)
      V1(NJJ+2)=.2
      V2(NJJ+1)=-.6-(1.5*.2)
      V2(NJJ+2)=.2
      IF(IBLAS .EQ. 1) THEN 
c     CALL DASHPLT(VVAL,YVAL2,NBL,1,PATTERN,2)
c     CALL CHARACT(2.,6.5,.2,'DASH = BLASIUS',0.,14,0)
      END IF
c     CALL LINPLT(V1,Y1,NJJ,1,0,0,1,1)
c     CALL DASHPLT(V2,Y2,NJJ,1,TATTERN,4) 
c     CALL CHARACT(2.,6.1,.2,'SOLID  = X AT ',0.,14,0)
c     CALL NUMBER(4.5,6.1,.2,XX(1),0.,3)
c     CALL CHARACT(2.,5.7,.2,'DTDASH = X AT ',0.,14,0)
c     CALL NUMBER(4.5,5.7,.2,XX(2),0.,3)
c     WRITE(8,'(I6)') NJJ 
c     DO 600 J=1,NJJ
c600  WRITE(8,'(2E15.8)') V1(J),Y1(J) 
c     WRITE(8,'(I6)') NJJ 
c     DO 601 J=1,NJJ
c601  WRITE(8,'(2E15.8)') V2(J),Y2(J) 
      IF(IBLAS .EQ. 1) THEN 
c     WRITE(8,'(I6)') NBL 
c     DO 602 J=1,NBL
c602  WRITE(8,'(2E15.8)') VVAL(J),YVAL2(J)
      END IF
c     CALL NFRAME 
 8888 CONTINUE
c     CALL CALPLT(0.,0.,999)
      WRITE(6,'(/,'' TECPLOT DATA OUTPUT TO flatplate.dat'')') 
c     WRITE(6,'('' TAPE8.DAT=Y VS V,  TAPE9.DAT=Y VS T'')') 
      STOP
      END 
