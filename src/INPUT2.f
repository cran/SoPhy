*
*  Authors:   Simunek, J., T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, Martin.Schlather@uni-bayreuth.de 
*  for Computers & Geosciences:
*     The use of the language interface of R: two examples for        
*      modelling water flux and solute transport                       
*                                                                     
*  Copyright (C) 2002 Simunek, J., T. Vogel and M. Th. van Genuchten
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
* 
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*
*
* NOTE: The original file was published as a public domain contribution.
*       The modified file is published under the GNU licence in
*       accordance with Rien van Genuchten, email correspondance August, 2002


* Source file INPUT2.FOR |||||||||||||||||||||||||||||||||||||||||||||||

      subroutine BasInf(KAT,MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,
     !                  SeepF,FluxF,FreeD,DrainF,
     !                  IntVec, NInt, DblVec, NDbl)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision DblVec(NDbl)
      integer IntVec(NInt)

      character*72 Hed
      character*5  LUnit,TUnit,MUnit
      logical lWat,lChem,AtmInF,ShortF,SeepF,FluxF,FreeD,DrainF

      dimension IU(11)
      data IU /50,71,72,75,76,77,78,79,80,81,82/

C      read(30,*)
C      read(30,*)
C      read(30,*) Hed
C      read(30,*)
C      read(30,*) LUnit,TUnit,MUnit
C      read(30,*)
C      read(30,*) KAT
      KAT = IntVec(1)
CC      write(51, 749) KAT
 749  format('KAT',i5)
C      read(30,*)
C      read(30,*) MaxIt,TolTh,TolH
      MaxIt = IntVec(2)
      TolTh = DblVec(1)
      TolH = DblVec(2)
CC      write(51, 748) MaxIt,TolTh,TolH 
 748  format("maxit",i5,"   TolTh",e15.8,"   TolH",e15.8)
C      read(30,*)
C      read(30,*) lWat,lChem,CheckF,ShortF,FluxF,AtmInF,SeepF,FreeD,
C     !           DrainF
      lWat = IntVec(3).ne.0
      lChem = IntVec(4).ne.0
      ShortF = IntVec(5).ne.0
      FluxF = IntVec(6).ne.0
      AtmInF = IntVec(7).ne.0
      SeepF = IntVec(8).ne.0
      FreeD = IntVec(9).ne.0
      DrainF = IntVec(10).ne.0
      
C      do 11 i=1,11
C        write(IU(i),*) Hed
C        write(IU(i),*)
C        write(IU(i),*)'Program SWMS_2D'
c        call getdat(ii,imonth,iday)
c        call gettim(ihours,mins,isecs,ii)
c        write(IU(i),100) iday,imonth,ihours,mins,isecs
C        if(AtmInF) then
C         write(*,*)'Time dependent boundary conditions'
C        else
C          write(*,*)'Time independent boundary conditions'
C        end if
C        if(KAT.eq.0) write(IU(i),110)
C        if(KAT.eq.1) write(IU(i),120)
C        if(KAT.eq.2) write(IU(i),130)
C        write(IU(i),*) 'Units: L = ',LUnit,', T = ',TUnit,', M = ',MUnit
C11    continue
C      write(*,*) '-----------------------------------------------------'
C      write(*,*) '|                                                   |'
C      write(*,*) '|                     SWMS_2D                       |'
C      write(*,*) '|                                                   |'
C      write(*,*) '|     Code for simulating water flow and solute     |'
C      write(*,*) '|       transport in two-dimensional variably       |'
C      write(*,*) '|             saturated porous media                |'
C      write(*,*) '|                                                   |'
C      write(*,*) '|                  version 1.22                     |'
C      write(*,*) '|          Last modified: January, 1994             |'
C      write(*,*) '|                                                   |'
C      write(*,*) '-----------------------------------------------------'
C      write(*,*)
C      write(*,*) Hed
C      if(KAT.eq.0) write(*,110)
C      if(KAT.eq.1) write(*,120)
C      if(KAT.eq.2) write(*,130)
C      write(51,140) MaxIt,TolTh,TolH
C      write(51, 746) lWat,lChem,ShortF, FluxF,AtmInF,
C     ! SeepF,FreeD,DrainF
 746  format('Wat', l3, ' lChem' ,l3,'  ShortF', l3,
     ! '  FluxF',l3, '  AtmInF',l3,'  SeepF',l3,'  FreeD',l3,
     ! '  DrainF',l3)

100   format(' Date: ',i3,'.',i2,'.','    Time: ',i3,':',i2,':',i2)
110   format(' Horizontal plane flow, V = L*L')
120   format(' Axisymmetric flow, V = L*L*L')
130   format(' Vertical plane flow, V = L*L')
140   format(/' Max. number of iterations           ',i4/
     !        ' Absolute water content tolerance [-]',f8.5/ 
     !        ' Absolute pressure head tolerance [L]',f8.5/)
      return
      end

************************************************************************

      subroutine MatIn(NMat,NLay,Par,hTab1,hTabN,
     !                 IntVec, NInt, DblVec, NDbl)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt)
      double precision DblVec(NDbl),  K      
      dimension Par(10,NMat),Qe(10)
      data Qe /1.,.99,.90,.85,.75,.65,.50,.35,.20,.10/

C      write(*,*) 'reading material information'
      Imax=10
C      read(30,*)
C      read(30,*)
C      read(30,*) NMat,NLay,hTab1,hTabN
      NLay = IntVec(11)
      hTab1 = DblVec(3)
      hTabN = DblVec(4)

      hTab1=-dmin1(abs(hTab1),abs(hTabN))
      hTabN=-dmax1(abs(hTab1),abs(hTabN))
C      read (30,*)
C      write(50,110)
C      do 11 M=1,NMat
C        read (30,*)     (Par(i,M),i=1,NPar) 
C        write(50,120) M,(Par(i,M),i=1,NPar)
C11    continue
C      write(51,110)
C      write(51,130)
      do 13 M=1,NMat
C        write(51,*)
        do 12 i=1,Imax
          h=FH(Qe(i),Par(1,M))
          K=FK(h,Par(1,M))
          C=FC(h,Par(1,M))
          Q=FQ(h,Par(1,M))
C          write(51,140) M,Qe(i),Q,h,C,K
12      continue
13    continue

110   format(/' MatNum,        thR    thS    tha    thm      alpha      
     !   n          Ks          Kk          thk'/)
120   format(i5,8x,4f7.3,16e12.3)
130   format(//' MatNum          Qe     Q        h         C         K') 
140   format(i5,8x,2f7.3,f10.3,e10.2,e12.3)
      return
      end

************************************************************************

      subroutine GenMat(NTab,NTabD,NMat,thR,hSat,Par,hTab,ConTab,CapTab,
     !                  ConSat,TheTab,thSat)

      IMPLICIT REAL*8 (A-H,O-Z)
      dimension Par(10,NMat),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          TheTab(NTabD,NMat),hTab(NTab),ConSat(NMat),hSat(NMat),
     !          thR(NMat),thSat(NMat)
      
C      write(*,*) 'generating materials'
      hTab1=hTab(1)
      hTabN=hTab(NTab)
      dlh=(log10(-hTabN)-log10(-hTab1))/(NTab-1)
      do 11 i=1,NTab
        alh=log10(-hTab1)+(i-1)*dlh
        hTab(i)=-10**alh
11    continue
      do 13 M=1,NMat
        Hr       =FH(0.0D0,Par(1,M))
        hSat(M)  =FH(1.0D0,Par(1,M))
        ConSat(M)=FK(0.0D0,Par(1,M))  
        thR(M)   =FQ(Hr ,Par(1,M))
        thSat(M) =FQ(0.0D0,Par(1,M))
        do 12 i=1,NTab
          ConTab(i,M)=FK(hTab(i),Par(1,M))
          CapTab(i,M)=FC(hTab(i),Par(1,M))
          TheTab(i,M)=FQ(hTab(i),Par(1,M))
12      continue
13    continue  
      return
      end

************************************************************************

      subroutine TmIn(tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !                TPrint,t,dtOpt,dtOld,AtmInF,
     !                IntVec, NInt, DblVec, NDbl, MPL)

      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt)
      double precision DblVec(NDbl)

      logical AtmInF
      dimension TPrint(MPL + 1)

C      write(*,*) 'reading time information'
C      read(30,*)
C      read(30,*)
C      read(30,*) dt,dtMin,dtMax,dMul,dMul2,MPL
      dt = DblVec(5)
      dtMin = DblVec(6)
      dtMax = DblVec(7)
      dMul = DblVec(8)
      dMul2 = DblVec(9)
C      write(*,748)dt , dtMin ,  dtMax, DMul , DMul2,MPL
 748  format('dt',f7.3,'   dtMin',f7.3, '   dtMax',f7.3,'   DMul', f7.3,
     ! '   DMul2' ,f7.3,'   MPL',i5)
C      read(30,*)
C      read(30,*) (TPrint(i),i=1,MPL)
C      do 703 i=1,MPL
C        TPrint(i) = TPrintR(i)
C 703  continue
C        write(51,749) (TPrint(i),i=1,MPL)
 749    format('TPr',100f8.2)
      dtOpt=dt
      dtOld=dt
      if(.not.AtmInF) then
        tMax=TPrint(MPL)
        tAtm=tMax
      end if
      TPrint(MPL+1)=tMax
      tOld=tInit
      t=tInit+dt
      return
      end

************************************************************************

      subroutine DrainIn(NDr,NDrD,NElDrD,NumEl,ND,NED,KElDr,EfDim,
     !                   ConAxx,ConAxz,ConAzz, 
     !                   IntVec, NInt, DblVec, NDbl,
     !                   NrErr)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt)
      double precision DblVec(NDbl)

      integer e
      dimension ND(NDrD),NED(NDrD),KElDr(NDrD,NElDrD),EfDim(2,NDrD),
     !          ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl)

C      write(*,*) 'reading drainage information'
C      read(30,*)
C      read(30,*)
C      read(30,*) NDr,DrCorr      
      NDr = IntVec(13)
      DrCorr = DblVec(10)
      if(NDr.gt.NDrD) then
C        write(*,*) 'Dimension in NDrD is exceeded'
         NrErr = 3
         return
        stop
      end if
C      read(30,*)
C      read(30,*) (ND(i),i=1,NDr)
C      read(30,*)
C      read(30,*) (NED(i),i=1,NDr)
       do 11 i=1,NDr
        if(NED(i).gt.NElDrD) then
C           write(*,*) 'Dimension in NElDrD is exceeded'
           NrErr = 4
           return
           stop
        end if
11    continue
C      read(30,*)
C      do 12 i=1,NDr
C        read(30,*) (EfDim(i,j),j=1,2)
C         do 704 j=1,2
C            EfDim(i,j) = EfDimR(i,j)
C 704     continue
C12    continue
C      read(30,*)
C      do 13 i=1,NDr
C        read(30,*) (KElDr(i,j),j=1,NED(i))
C         do 705 j=1,NED(i)
C            KElDr(i,j) = KElDrR(i,j)
C 705     continue         
C13    continue
      do 15 i=1,Ndr
        rho=EfDim(i,2)/EfDim(i,1)
        A=(1.+0.405*rho**(-4))/(1.-0.405*rho**(-4))
        B=(1.+0.163*rho**(-8))/(1.-0.163*rho**(-8))
        C=(1.+0.067*rho**(-12))/(1.-0.067*rho**(-12))
        Red=376.7/(138.*log10(rho)+6.48-2.34*A-0.48*B-0.12*C)/DrCorr
        do 14 j=1,NED(i)
          e=KElDr(i,j)
          ConAxx(e)=ConAxx(e)*Red
          ConAxz(e)=ConAxz(e)*Red
          ConAzz(e)=ConAzz(e)*Red
14      continue
15    continue
      return
      end

************************************************************************

      subroutine NodInf(NumNP,NumEl,NumBP,NumKD,
     !                  NObs,NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,y,
     !                  MatNum,Beta,Axz,Bxz,Dxz,
     !                  IntVec, NInt, nCodeM, d1doub, 
     !                  NrErr
     !                 )
      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt), nCodeM(NumNP, 3)
      double precision d1doub(NumNP, 9)

      dimension Kode(NumNP),Q(NumNP),hOld(NumNP),x(NumNP),y(NumNP),
     !          hNew(NumNP),hTemp(NumNP),MatNum(NumNP),Beta(NumNP),
     !          Axz(NumNP),Bxz(NumNP),Dxz(NumNP),Conc(NumNP)

C      write(*,*) 'reading nodal information'
C      read(32,*)
C      read(32,*)
C      read(32,*) NumNP,NumEl,IJ,NumBP,NObs
      ii = 0
      NObs  = IntVec(14)
C      write(51,745) NumNP, NumEl, IJ, NumBP, NObs
 745  format('NumNP', i5, '  NumEl', i5, '  IJ', i5, '  NumBP', i5, 
     ! '  NObs',i5)
      if(NumNP.gt.NumNP) then
C         write(*,*) 'Dimension in NumNP is exceeded'
         NrErr = 5
         return
        stop
      else if(NObs.gt.NObsD) then
C        write(*,*) 'Dimension in NObsD is exceeded'
        stop
         NrErr = 8
         return
      end if
C      read(32,*)
      NPR=0
      k=0

11    k=k+1
C      read(32,*) n,Kode(n),x(n),y(n),hOld(n),Conc(n),Q(n),MatNum(n),
C     !           Beta(n),Axz(n),Bxz(n),Dxz(n)
      ii = ii + 1
      n = nCodeM(ii, 1)
      Kode(n) = nCodeM(ii, 2)
      MatNum(n) = nCodeM(ii, 3)
      x(n) = d1doub(ii, 1)
      y(n) = d1doub(ii, 2)
      hOld(n) = d1doub(ii, 3)
      Conc(n) = d1doub(ii, 4)
      Q(n) = d1doub(ii, 5)
      Beta(n) = d1doub(ii, 6)
      Axz(n) = d1doub(ii, 7)
      Bxz(n) = d1doub(ii, 8)
      Dxz(n) = d1doub(ii, 9)
 
      if(Kode(n).gt.NumKD) then
C       write(*,*) 'Dimension in NumKD is exceeded'
        NrErr=9
        return
        stop
      end if
      if(n-k) 12,15,13
12    continue
C     write(*,130) n
      NrErr=10
      return
      stop
13    Deno=n-k+1
      DX=(x(n)-x(NPR))/Deno
      DY=(y(n)-y(NPR))/Deno
      DP=(hOld(n)-hOld(NPR))/Deno
      DConc=(Conc(n)-Conc(NPR))/Deno
      DBeta=(Beta(n)-Beta(NPR))/Deno
      DA=(Axz(n)-Axz(NPR))/Deno
      DB=(Bxz(n)-Bxz(NPR))/Deno
      DD=(Dxz(n)-Dxz(NPR))/Deno
14    x(k)=x(k-1)+DX
      y(k)=y(k-1)+DY
      hOld(k)=hOld(k-1)+DP
      Conc(k)=Conc(k-1)+DConc
      Beta(k)=Beta(k-1)+DBeta
      Axz(k)=Axz(k-1)+DA
      Bxz(k)=Bxz(k-1)+DB
      Dxz(k)=Dxz(k-1)+DD
      MatNum(k)=MatNum(k-1)
      Kode(k)=Kode(k-1)
      Q(k)=Q(k-1)
      k=k+1
      if(k.lt.N) goto 14
15    NPR=N

      if(k.lt.NumNP) goto 11

      do 16 n=1,NumNP
        hNew(n)=hOld(n)
        hTemp(n)=hOld(n)
16    continue

110   format(////' NODAL POINT INFORMATION'//'  NODE NO.',6x,'KODE',
     !       7x,'X,R',12x,'Y,Z',11x,'.PSI.',12x,'Q',11x,'Conc'/)
C120   format(2i10,5e15.6,i10,4f7.3)
120   format(2i4,5f8.2,i2,4f7.3)
130   format(' ERROR IN NodInf AT N=',i5)
      return
      end

************************************************************************

      subroutine ElemIn(NumEl,NumNP,KX,LayNum,ConAxx,ConAzz,
     !                  ConAxz,ListNE,MBand,MBandD,lChem,lOrt,
     !                  IntVec,NInt,ConAX,KXR
     !                  )
      IMPLICIT REAL*8 (A-H,O-Z)    
      integer IntVec(NInt), KXR(NumEl, 6)
      double precision ConAX(NumEl, 3)

      logical lChem,lConst,lOrt
      integer e
      dimension KX(NumEl,4),ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),
     !          LayNum(NumEl),ListNE(NumNP)

C     write(*,*) 'reading element information'
      Num=0
C      read(32,*)
C      read(32,*)
      ii = 0
      do 14 e=1,NumEL
        IF (Num-e) 11,14,12
 11     continue
C     read(32,*) Num,(KX(Num,i),i=1,4),ConAxz(Num),ConAxx(Num),
C     !             ConAzz(Num),LayNum(Num)
        ii = ii + 1
        Num = KXR(ii, 1)
        KX(Num, 1) = KXR(ii, 2)
        KX(Num, 2) = KXR(ii, 3)
        KX(Num, 3) = KXR(ii, 4)
        KX(Num, 4) = KXR(ii, 5)        
        ConAxz(Num) = ConAX(ii, 1)
        ConAxx(Num) = ConAX(ii, 2)
        ConAzz(Num) = ConAX(ii, 3)
        LayNum(Num) = KXR(ii, 6)
        if(KX(Num,4).eq.0) KX(Num,4)=KX(Num,3)
        if(Num.eq.e) goto 14
12      do 13 i=1,4
          KX(e,i)=KX(e-1,i)+1
13      continue
        ConAxx(e)=ConAxx(e-1)
        ConAzz(e)=ConAzz(e-1)
        ConAxz(e)=ConAxz(e-1)
        LayNum(e)=LayNum(e-1)
14    continue
      AA=3.141592654/180.
      do 15 e=1,NumEl
        Ang=AA*ConAxz(e)
        CAxx=ConAxx(e)
        CAzz=ConAzz(e)
        ConAxx(e)=CAxx*cos(Ang)*cos(Ang) + CAzz*sin(Ang)*sin(Ang)
        ConAzz(e)=CAxx*sin(Ang)*sin(Ang) + CAzz*cos(Ang)*cos(Ang)
        ConAxz(e)=(CAxx-CAzz)*sin(Ang)*cos(Ang)
15    continue

      do 17 i=1,NumNP
        ListNE(i)=0
17    continue
      do 19 e=1,NumEl
        NCorn=4
        if(KX(e,3).eq.KX(e,4)) NCorn=3
        do 18 n=1,NCorn-2
          i=KX(e,1)
          j=KX(e,n+1)
          k=KX(e,n+2)
          ListNE(i)=ListNE(i)+1
          ListNE(j)=ListNE(j)+1
          ListNE(k)=ListNE(k)+1
18      continue
19    continue

      lOrt=.false.
      lConst=.true.
      MBand=1
      do 21 e=1,NumEl
        NUS=4
        if(KX(e,3).eq.KX(e,4)) NUS=3
        do 20 kk=1,NUS-2
          MB=1
          i=KX(e,1)
          j=KX(e,kk+1)
          k=KX(e,kk+2)
          if(abs(i-j).gt.MB) MB=abs(i-j)
          if(abs(i-k).gt.MB) MB=abs(i-k)
          if(abs(j-k).gt.MB) MB=abs(j-k)
          if(MB.gt.MBand) MBand=MB
          if(e.eq.1.and.kk.eq.1) then
            MB1=MB
          else
            if(MB1.ne.MB) lConst=.false.
          end if
20      continue
21    continue
      MBand=MBand+1
      if(MBand.gt.MBandD.or.(lChem.and.2*MBand-1.gt.MBandD)) lOrt=.true.
 748  format("mband", i6, i6, l5)
C      write(*, 748) MBand, MBandD, lChem
C      if(.not.lConst) IJ=NumNP
      if(MBand.gt.10.or.NumNP.gt.200) lOrt=.true.

110   format (////' Element Information'//' Element    C O R N E R    N 
     !O D E S   ConAxz    ConAxx    ConAzz  LayNum'/)
120   format (i6,i9,3i6,e14.3,2f8.3,i5)
      return
      end

************************************************************************

      subroutine GeomIn(NumKD,NumNP,NumBP,NObs,NObsD,SWidth,Width,Kode,
     !                  KXB,rLen,Node, 
     !                  DblVec,NDbl)
      IMPLICIT REAL*8 (A-H,O-Z)
      double precision DblVec(NDbl)
      
      character*50 Text1
      dimension KXB(NumBP),Width(NumBP),SWidth(NumKD),Kode(NumNP),
     !          Node(NObsD)

C     write(*,*) 'reading geometric information'
C      read(32,*)
C      read(32,*)
C      read(32,*) (KXB(i),i=1,NumBP)
C      read(32,*)
C      read(32,*) (Width(i),i=1,NumBP)
C      do 708 i=1,NumBP
C         Width(i) = WidthR(i)
C 708  continue      
C      write(51,754)  (KXB(i), i=1,NumBP)
 754  format('Nodes',100i5)
C      write(51,744) (Width(i), i=1,NumBP)
 744  format('Width',100f7.3)

C      read(32,*)
C      read(32,*) rLen
      rlen = DblVec(11)
C         write(51,745) rlen
 745     format('LENGTH',f7.3)
      do 11 i=1,NumKD
        SWidth(i)=0.
11    continue
      do 12 i=1,NumBP
        n=KXB(i)
        j=iabs(Kode(n))
        if(j.eq.0) goto 12
        SWidth(j)=SWidth(j)+Width(i)
12    continue
C      if(NObs.gt.0) then
C        read(32,*)
C        read(32,*) (Node(i),i=1,NObs)
C        Text1='    hNew     theta     conc    '
C       write(92,110) (Node(j),j=1,NObs)
C       write(92,120) (Text1,i=1,NObs)
C      end if
110   format (///14x,5(11x,'Node(',i3,')',11x))
120   format (/'     time     ',5(a31))
      return
      end

************************************************************************

      subroutine AtmIn(GWL0L,SinkF,qGWLF,tInit,tMax,Aqh,Bqh,hCritS,
     !                 MaxAL, 
     !                 IntVec, NInt, DblVec, NDbl
     !                )

      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt)
      double precision DblVec(NDbl)

      logical SinkF,qGWLF

C     write(*,*) 'reading atmospheric information'
C      read(31,*)
C      read(31,*)
C      read(31,*)
C      read(31,*)
C      read(31,*) SinkF,qGWLF
      SinkF = IntVec(15).NE.0
      qGWLF = IntVec(16).NE.0
C      read(31,*)
C      read(31,*) GWL0L,Aqh,Bqh
      GWL0L = DblVec(12)
      Aqh = DblVec(13)
      Bqh = DblVec(14)
C      read(31,*)
C      read(31,*) tInit,MaxAL
      tInit = DblVec(15)
C      read(31,*)
C      read(31,*) hCritS
      hCritS = DblVec(16)
C      read(31,*)
C      do 11 i=1,MaxAL-1
C        read(31,*)
C11    continue
C      read(31,*) tMax
      tMax = DblVec(17)
C      rewind 31
C      do 12 i=1,12
C        read(31,*)
C     12    continue 
C      write(51, 740)
C      write(51, 741) SinkF,qGWLF,GWL0L,Aqh,Bqh,tInit,MaxAL,hCritS,tMax
 740  format('   SinkF   qGWLF   GWL0L     Aqh      Bqh   tInit',
     !       '   MaxAL  hCritS   tMax')
 741  format(2l8,4f8.3,i8,2e8.1)
      return
      end

************************************************************************

      subroutine SinkIn(NumEl,NumNP,KAT,KX,x,y,Beta,DblVec,NDbl)
      IMPLICIT REAL*8 (A-H,O-Z)
      double precision  DblVec(NDbl)

      integer e
      dimension Beta(NumNP),KX(NumEl,4),x(NumNP),y(NumNP)

C     write(*,*) 'reading sink information'
C      read(30,*)
C      read(30,*)
C      read(30,*) P0,P2H,P2L,P3,r2H,r2L
CC      P0 = DblVec(18)
CC      P2H = DblVec(19)
CC      P2L = DblVec(20)
CC      P3 = DblVec(21)
CC      r2H = DblVec(22)
CC      r2L = DblVec(23)
C      read(30,*)
C      read(30,*) (POptm(i),i=1,NMat)
C      write(51, 742) P0, P2H, P2L, P3, r2H, r2L
 742  format("P0", f12.3, "  P2H", f12.3, "  P2L", f12.3, "  P3", f12.3,
     !       "  r2H", f12.3, "  r2L", f12.3)
C      write(51, 741) (POptm(i),i=1,NMat)
 741  format(100f12.3)
C see also swms2d.R
C      P0 =-abs(P0)
C      P2L=-abs(P2L)
C      P2H=-abs(P2H)
C      P3 =-abs(P3)
      xMul=1.
      SBeta=0.
      do 12 e=1,NumEl
        NUS=4
        IF(KX(e,3).eq.KX(e,4)) NUS=3
        do 11 k=1,NUS-2
          i=KX(e,1)
          j=KX(e,k+1)
          l=KX(e,k+2)
          CJ=x(i)-x(l)
          CK=x(j)-x(i)
          BJ=y(l)-y(i)
          BK=y(i)-y(j)
          AE=(CK*BJ-CJ*BK)/2.
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          BetaE=(Beta(i)+Beta(j)+Beta(l))/3.
          SBeta=SBeta+xMul*AE*BetaE
11      continue
12    continue
      do 13 i=1,NumNP
        Beta(i)=Beta(i)/SBeta
13    continue
      return
      end

************************************************************************

      subroutine ChemIn(NMat,NumBP,cBound,epsi,tPulse,KodCB,NLevel,
     !                  lUpW,lArtD,PeCr,IntVec, NInt, DblVec, NDbl)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt)
      double precision  DblVec(NDbl)

      logical lUpW,lArtD
      dimension KodCB(NumBP),cBound(6)

C     write(*,*) 'reading transport information'
      NLevel=1
      write(51,110)
C      read(30,*)
C      read(30,*)
C      read(30,*) epsi,lUpW,lArtD,PeCr
      epsi = DblVec(18)
      lUpW = IntVec(17).ne.0
      lArtD= IntVec(18).ne.0
      PeCr = DblVec(19)
      
C      write(51, 750) epsi, lUpW, lArtD, PeCr
 750  format('epsi',1e11.3,'   lUpW/lArtD',2l3,'    PeCr',e11.3)

      PeCr=dmax1(PeCr,0.01D0)
      if(epsi.lt.0.999) NLevel=2      
C      read(30,*) 
       if(lUpW) then
C          write(51,120)
      else
C        write(51,130)
C        if(lArtD) write(51,140) PeCr
      end if
C      write(51,150)
C      do 11 M=1,NMat    
C        read(30,*) (ChPar(j,M),j=1,9)
C       write(51,160) (ChPar(j,M),j=1,10)
C11    continue   

C      read(30,*)
C      read(30,*) (KodCB(i),i=1,NumBP)
C      read(30,*)
C      read(30,*) (cBound(i),i=1,6)
      do 711 i=1,6
         cBound(i)=DblVec(19 + i)
 711  continue
C      write(51,751) (KodCB(i),i=1,NumBP)
 751  format(' KodCD ',100i5)
C      write(51,170) (cBound(i),i=1,6)
C      read(30,*)
C      read(30,*) tPulse
      tPulse = DblVec(26)
C      write(51,180) tPulse
  
110   format (/' Solute transport information'/1x,28('='))
120   format (/' Upstream weighting finite-element method')
130   format (/' Galerkin finite-element method')
140   format (/' Artificial dispersion is added when Peclet number is',
     !         ' higher than',f10.3)
150   format (/'  Bulk.d.   Difus.         Disper.        Adsorp.    ',
     !         'SinkL1    SinkS1    SinkL0    SinkS0')
160   format (10e10.3) 
170   format (/'   Conc1     Conc2     Conc3     Conc4     cSink     ',
     !         'cWell'/8e10.3)
180   format (/' tPulse =   ',f15.3)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
