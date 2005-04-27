*
*  Authors:   Simunek, J., T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, schlath@hsu-hh.de 
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

* NOTE: The original file was published as a public domain contribution.
*       The modified file is published under the GNU licence in
*       accordance with Rien van Genuchten, email correspondance August, 2002


* Source file WATFLOW2.FOR |||||||||||||||||||||||||||||||||||||||||||||

      subroutine WatFlow(NumNP,NumEl,NTab,NTabD,MBand,MBandD,
     !                   NMat,NSeep,NSeepD,NumSPD,NSP,NP,NumBP,ItCum,
     !                   MaxIt,Iter,Kode,KAT,t,dt,dtMin,dtOpt,dtOld,
     !                   tOld,hCritA,hCritS,TolTh,TolH,rLen,Width,rTop,
     !                   vMeanR,hMeanR,AtmInf,SinkF,SeepF,qGWLF,FreeD,
     !                   Par,hTab,ConTab,CapTab,TheTab,hNew,hOld,hTemp,
     !                   thR,thSat,ThNew,ThOld,MatNum,Con,Cap,ConSat,
     !                   Axz,Bxz,Dxz,hSat,A,B,Q,F,x,y,KX,Sink,DS,Beta,
     !                   ConAxx,ConAzz,ConAxz,KXB,Explic,GWL0L,Aqh,Bqh,
     !                   lWat,TLevel,lOrt,DrainF,ND,NDr,NDrD,rRoot,PRL,
     !                   ConO,
     !                   A1,B1,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,
     !                   QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO, 
     !                   NrErr)

      IMPLICIT REAL*8 (A-H,O-Z)
      logical AtmInf,SinkF,SeepF,Explic,ItCrit,lWat,qGWLF,FreeD,lOrt,
     !        DrainF
      double precision A,B,A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      integer TLevel
      dimension A(MBandD,NumNP),B(NumNP),Kode(NumNP),Q(NumNP),F(NumNP),
     !          hNew(NumNP),hTemp(NumNP),hOld(NumNP),ConSat(NMat),
     !          hTab(NTab),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          Con(NumNP),Cap(NumNP),x(NumNP),y(NumNP),MatNum(NumNP),
     !          KX(NumEl,4),KXB(NumBP),Par(10,NMat),Width(NumBP),
     !          ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),hSat(NMat),
     !          NP(NSeepD,NumSPD),NSP(NSeepD),DS(NumNP),Beta(NumNP),
     !          Axz(NumNP),Bxz(NumNP),Dxz(NumNP),Sink(NumNP),thR(NMat),
     !          thSat(NMat),TheTab(NTabD,NMat),ThNew(NumNP),ND(NDrD),
     !          ThOld(NumNP),ConO(NumNP),
     !          A1(MBandD,NumNP),B1(NumNP),RES(NumNP),IAD(MBandD,NumNP),
     !          IADN(NumNP),IADD(NumNP),VRV(NumNP),RQI(NumNP,MNorth),
     !          RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNP,MNorth),
     !          PRL(NumNP, 7)

      logical ldebug   
      ldebug = .false.
C      ldebug = .true.
C      write(*,800) 1, (hNew(iii), iii=1,42)
 800  format(i8,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",
     !1e14.7,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",
     !1e14.7,"\n")

      if (ldebug) write(*,*) "W.A"
      if(lWat.and.TLevel.ne.1) then
        do 10 i=1,NumNP
          hOld(i) =hNew(i)
          ThOld(i)=ThNew(i)
          ConO(i)=Con(i)
          if(Kode(i).lt.1) then
            hTemp(i)=hNew(i)+(hNew(i)-hOld(i))*dt/dtOld
            hNew(i) =hTemp(i)
          else
            hTemp(i)=hNew(i)
          end if
10      continue
      end if

11    continue

      if (ldebug) write(*,*) "W.B"
      Iter=0
      Explic=.false.

* --- Beginning of iteration loop --------------------------------------

12    continue
      if (ldebug) write(*,*) "W.C"
      if(SinkF) !.and..not.lWat.and.Iter.ne.0)
     !  call SetSnk(NumNP,hNew,rRoot,Sink,PRL,Beta,rLen)
      if (ldebug) write(*,*) "W.D"
      call SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !            MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,hTemp,
     !            Explic,TheTab,thSat,thR,ThNew)
C      write(*,800) 2,(hNew(iii), iii=1,42)
      if (ldebug) write(*,*) "W.E"
      call Reset (A,B,Kode,Q,hNew,F,Con,Cap,x,y,KX,NumNP,NumEl,
     !            MBand,MBandD,KAT,dt,Iter,SinkF,Sink,DS,Beta,rLen,
     !            vMeanR,hMeanR,ConAxx,ConAzz,ConAxz,ThNew,ThOld,lWat,
     !            PRL,lOrt,IAD,IADN,IADD,B1)
C      write(*,800) 3,(hNew(iii), iii=1,42)
C      write(*,800) 31,(B1(iii), iii=1,42)
      if (ldebug) write(*,*) "W.F"
      call Shift (NumNP,NumBP,NSeepD,NumSPD,NSeep,NSP,NP,hNew,hOld,Q,
     !            Kode,rTop,Width,KXB,hCritA,hCritS,SeepF,AtmInF,Explic,
     !            qGWLF,FreeD,GWL0L,Aqh,Bqh,Con,ConO,DrainF,ND,NDr,NDrD,
     !            lWat,Iter)
C      write(*,800) 4,(hNew(iii), iii=1,42)
      if (ldebug) write(*,*) "W.G"
      call Dirich(A,B,Kode,hNew,NumNP,MBand,MBandD,lOrt,IADD)
      if (ldebug) write(*,*) "W.H"
      if(lOrt) then
        call ILU (A,NumNP,MBandD,IAD,IADN,IADD,A1)
        North=0
C        write(*,*) "WATERFLOW orthomin"
      if (ldebug) write(*,*) "W.|"
        call OrthoMin(A,B1,B,NumNP,MBandD,IAD,IADN,IADD,A1,VRV,
     !                RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,RCNVRG,ACNVRG,
     !                North,MNorth,MaxItO, NrErr)
        if (NrErr.ne.0) return
C      write(*,800) 31,(B1(iii), iii=1,42)
      else
      if (ldebug) write(*,*) "W.J"
        call Solve (A,B,NumNP,MBand,MBandD)
      end if

C      write(*,800) 5,(hNew(iii), iii=1,42)
C      write(*,800) 51,(B(iii), iii=1,42)
C      write(*,800) 52,(B1(iii), iii=1,42)
C      write(*,801) lOrt
 801  format(l8)
      do 13 i=1,NumNP
        hTemp(i)=hNew(i)
        if(lOrt) B(i)=B1(i)
        hNew(i)=sngl(B(i))
        if(abs(hNew(i)).gt.1.d+10) hNew(i)=sign(1.d+10,hNew(i))
13    continue
      Iter =Iter+1
      ItCum=ItCum+1
      if(Explic) goto 18

      if (ldebug) write(*,*) "W.K"
C      write(*,800) 6,(hNew(iii), iii=1,42)
*     Test for convergence
      ItCrit=.true.
      do 14 i=1,NumNP
        m=MatNum(i)
        EpsTh=0.
        EpsH=0.
        if(hTemp(i).lt.hSat(m).and.hNew(i).lt.hSat(m)) then
          Th=ThNew(i)+Cap(i)*(hNew(i)-hTemp(i))/(ThSat(m)-ThR(m))/Dxz(i)
          EpsTh=abs(ThNew(i)-Th)
        else
          EpsH=abs(hNew(i)-hTemp(i))
        end if
        if(EpsTh.gt.TolTh.or.EpsH.gt.TolH) then
          ItCrit=.false.
          goto 15
        end if
14    continue
15    continue
      if (ldebug) write(*,*) "W.L"

C      write(*,800) 7,(hNew(iii), iii=1,42)
      if(.not.ItCrit) then
        if(Iter.lt.MaxIt.or.(.not.lWat.and.Iter.lt.5*MaxIt)) then
          goto 12
        else if(dt.le.dtMin) then
          Explic=.true.
          do 16 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
16        continue
          goto 12
        else if(.not.lWat) then
C         write(*,*) ' No steady state solution found'
          NrErr = 11
          return
C          stop
        else
          do 17 i=1,NumNP
            hNew(i) =hOld(i)
            hTemp(i)=hOld(i)
17        continue
          dt=dmax1(dt/3.,dtMin)
          dtOpt=dt
          t=tOld+dt
          goto 11
        end if
      end if
* --- End of iteration loop --------------------------------------------

      if (ldebug) write(*,*) "W.M"
18    continue

      if(.not.lWat.and.TLevel.eq.1) then
C       write(*,101) Iter
        do 19 i=1,NumNP
          hOld(i) =hNew(i)
          ThOld(i)=ThNew(i)
19      continue
      end if

101   format(' Steady state was reached after',i6,' iterations.')
      return
      end

************************************************************************

      subroutine Reset(A,B,Kode,Q,hNew,F,Con,Cap,x,y,KX,NumNP,NumEl,
     !                 MBand,MBandD,KAT,dt,Iter,SinkF,Sink,DS,
     !                 Beta,rLen,vMeanR,hMeanR,ConAxx,ConAzz,ConAxz,
     !                 ThNew,ThOld,lWat,PRL,lOrt,IAD,IADN,IADD,B1)

      IMPLICIT REAL*8 (A-H,O-Z)
      logical SinkF,lWat,lOrt
      double precision A,B,B1
      dimension A(MBandD,NumNP),B(NumNP),Q(NumNP),hNew(NumNP),F(NumNP),
     !          Con(NumNP),Cap(NumNP),x(NumNP),y(NumNP),KX(NumEl,4),
     !          ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),Kode(NumNP),
     !          Sink(NumNP),DS(NumNP),Beta(NumNP),ThNew(NumNP),
     !          ThOld(NumNP),E(3,3),iLoc(3),Bi(3),Ci(3),
     !          B1(NumNP),IAD(MBandD,NumNP),IADN(NumNP),IADD(NumNP),
     !          PRL(NumNP, 7)

*     Initialisation
      xMul=1.
      AreaR =0.
      if(Iter.eq.0) then
        vMeanR=0.
        hMeanR=0.
      end if
      do 12 i=1,NumNP
        B(i)=0.d0
        if(lOrt) B1(i)=hNew(i)
        F(i)=0.
        if(Iter.eq.0) DS(i)=0.
        do 11 j=1,MBandD
          A(j,i)=0.d0
11      continue
12    continue

*     Loop on elements
      do 16 n=1,NumEl
        CondI=ConAxx(n)
        CondJ=ConAzz(n)
        CondK=ConAxz(n)
        NUS=4
        if(KX(n,3).eq.KX(n,4)) NUS=3

*       Loop on subelements
        do 15 k=1,NUS-2
          i=KX(n,1)
          j=KX(n,k+1)
          l=KX(n,k+2)
C          write(*,777) i,j,k
 777      format(3i15)
          iLoc(1)=1
          iLoc(2)=k+1
          iLoc(3)=k+2
          Ci(1)=x(l)-x(j)
          Ci(2)=x(i)-x(l)
          Ci(3)=x(j)-x(i)
          Bi(1)=y(j)-y(l)
          Bi(2)=y(l)-y(i)
          Bi(3)=y(i)-y(j)
          AE=(Ci(3)*Bi(2)-Ci(2)*Bi(3))/2.
          ConE=(Con(i)+Con(j)+Con(l))/3.
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.
          AMul=xMul*ConE/4./AE
          BMul=xMul*ConE/2.
          FMul=xMul*AE/12.
          BetaE=(Beta(i)+Beta(j)+Beta(l))/3.
          if(SinkF.and.BetaE.gt.0..and.Iter.eq.0) then
            SinkE=(Sink(i)+Sink(j)+Sink(l))/3.
C PRL(i,4) is P3
            if(hNew(i).gt.PRL(i,4)) DS(i)=DS(i)+FMul*(3.*SinkE+Sink(i))
            if(hNew(j).gt.PRL(j,4)) DS(j)=DS(j)+FMul*(3.*SinkE+Sink(j))
            if(hNew(l).gt.PRL(l,4)) DS(l)=DS(l)+FMul*(3.*SinkE+Sink(l))
            hNewE=(hNew(i)+hNew(j)+hNew(l))/3.
            vMeanR=vMeanR+xMul*AE*SinkE/rLen
            hMeanR=hMeanR+xMul*hNewE*AE
            AreaR=AreaR+xMul*AE
          end if
          do 14 i=1,3
            iG=KX(n,iLoc(i))
C          write(*,778) 111, ig
            F(iG)=F(iG)+FMul*4.
            if(KAT.ge.1) B(iG)=B(iG)+BMul*(CondK*Bi(i)+CondJ*Ci(i))
            do 13 j=1,3
              jG=KX(n,iLoc(j))
C          write(*,778) 222, JG
 778      format(2i15)

              E(i,j)=CondI*Bi(i)*Bi(j)+CondK*(Bi(i)*Ci(j)+Ci(i)*Bi(j))+
     !               CondJ*Ci(i)*Ci(j)
              if(lOrt) then
                call Find(iG,jG,kk,NumNP,MBandD,IAD,IADN)
C                write(*,749)  i,j,kk,iG,A(kk,iG), AMul, E(i,j)
 749            format("wa",4i6,3e20.11)
                A(kk,iG)=A(kk,iG)+AMul*E(i,j)
              else
                iB=iG-jG+1
C          write(*,778) 333, ib
                if(iB.ge.1) A(iB,jG)=A(iB,jG)+AMul*E(i,j)
              end if
13          continue
14        continue
15      continue

16    continue
      if(Iter.eq.0.and.AreaR.gt.0.) hMeanR=hMeanR/AreaR
*     Determine boundary fluxes
      do 19 n=1,NumNP
        if(Kode(n).lt.1) goto 19
        QN=sngl(B(n))+DS(n)+F(n)*(ThNew(n)-ThOld(n))/dt
        if(lOrt) then
          do 17 j=1,IADN(n)
            QN=QN+sngl(A(j,n))*hNew(IAD(j,n))
17        continue
        else
          QN=QN+sngl(A(1,n))*hNew(n)
          do 18 j=2,MBand
            k=n-j+1
            if(k.ge.1) then
              QN=QN+sngl(A(j,k))*hNew(k)
            end if
            k=n+j-1
            if(K.le.NumNP) then
              QN=QN+sngl(A(j,n))*hNew(k)
            end if
18        continue
        end if
        Q(n)=QN
19    continue

*     Complete construction of RHS vector and form effective matrix
      do 20 i=1,NumNP
        if(.not.lWat) F(i)=0.
        j=1
        if(lOrt) j=IADD(i)
        A(j,i)=A(j,i)+F(i)*Cap(i)/dt
        B(i)=F(i)*Cap(i)*hNew(i)/dt-F(i)*(ThNew(i)-ThOld(i))/dt+
     !       Q(i)-B(i)-DS(i)
20    continue
      return
      end

************************************************************************

      subroutine Dirich(A,B,Kode,hNew,NumNP,MBand,MBandD,lOrt,IADD)

      IMPLICIT REAL*8 (A-H,O-Z)
      logical lOrt
      double precision A,B
      dimension A(MBandD,NumNP),B(NumNP),Kode(NumNP),hNew(NumNP),
     !          IADD(NumNP)

      do 12 n=1,NumNP
        if(Kode(n).lt.1) goto 12
        if(lOrt) then
          A(IADD(n),n)=10.d30
          B(n)=10.d30*hNew(n)
        else
          do 11 m=2,MBand
            k=n-m+1
            if(k.gt.0) then
              B(k)=B(k)-A(m,k)*hNew(n)
              A(m,k)=0.d0
            end if
            l=n+m-1
            if(NumNP-l.ge.0) then
              B(l)=B(l)-A(m,n)*hNew(n)
            end if
            A(m,n)=0.d0
11        continue
          A(1,n)=1.d0
          B(n)=hNew(n)
        end if
12    continue
      return
      end

************************************************************************

      subroutine Shift(NumNP,NumBP,NSeepD,NumSPD,NSeep,NSP,NP,hNew,hOld,
     !                 Q,Kode,EI,Width,KXB,hCritA,hCritS,SeepF,AtmInF,
     !                 Explic,qGWLF,FreeD,GWL0L,Aqh,Bqh,Con,ConO,DrainF,
     !                 ND,NDr,NDrD,lWat,Iter)

      IMPLICIT REAL*8 (A-H,O-Z)
      logical SeepF,AtmInF,Explic,qGWLF,FreeD,DrainF,lWat
      dimension Kode(NumNP),Q(NumNP),hNew(NumNP),hOld(NumNP),
     !          Width(NumBP),KXB(NumBP),NP(NSeepD,NumSPD),NSP(NSeepD),
     !          Con(NumNP),ConO(NumNP),ND(NDrD)

      logical ldebug
      ldebug = .false.
      ldebug = .true.

*     Modify conditions on seepage faces
      if(SeepF) then
        do 12 i=1,NSeep
c          iCheck=0
          NS=NSP(i)
          do 11 j=1,NS
            n=NP(i,j)
C            write(*,710) i,NSeep,j,NS,n
C 710        format(5i15)
            if(Kode(n).eq.-2) then
              if(hNew(n).lt.0.) then
c                iCheck=1
              else
                Kode(n)=2
                hNew(n)=0.
              end if
            else
c              if(iCheck.gt.0.or.Q(n).ge.0.) then
              if(Q(n).ge.0.) then
                Kode(n)=-2
                Q(n)=0.
c                iCheck=1
              end if
            end if
11        continue
12      continue
      end if

*     Modify conditions in drainage node
      if(DrainF) then
        do 13 i=1,NDr
          n=ND(i)
          if(Kode(n).eq.-5) then
            if(hNew(n).ge.0.) then
              Kode(n)=5
              hNew(n)=0.
            end if
          else
            if(Q(n).ge.0.) then
              Kode(n)=-5
              Q(n)=0.
            end if
          end if
13      continue
      end if

*     Modify potential surface flux boundaries
      if(AtmInF) then
        do 14 i=1,NumBP
          n=KXB(i)
          k=Kode(n)
          if(Explic.and.iabs(k).eq.4) then
            Kode(n)=-iabs(k)
            goto 14
          end if

*         Critical surface pressure on ...
          if(K.eq.4) then
            if(abs(Q(n)).gt.abs(-EI*Width(i)).or.Q(n)*(-EI).le.0) then
              Kode(n)=-4
              Q(n)=-EI*Width(i)
            end if
            goto 14
          end if

*         Surface flux on ...
          if(K.eq.-4.and.Iter.ne.0) then
            if(hNew(n).le.hCritA) then
              Kode(n)=4
              hNew(n)=hCritA
              goto 14
            end if
            if(hNew(n).ge.hCritS) then
              Kode(n)=4
              hNew(n)=hCritS
            end if
          end if

*         Time variable flux
          if(K.eq.-3) then
            if(lWat) then
              if(qGWLF) Q(n)=-Width(i)*Fqh(hOld(n)-GWL0L,Aqh,Bqh)
            else
              if(qGWLF) Q(n)=-Width(i)*Fqh(hNew(n)-GWL0L,Aqh,Bqh)
            end if
          end if
14      continue
      end if

*     Free Drainage
      if(FreeD) then      
        do 15 i=1,NumBP
          n=KXB(i)
          k=Kode(n)
          if(lWat) then
            if(k.eq.-3) Q(n)=-Width(i)*ConO(n)
          else
            if(k.eq.-3) Q(n)=-Width(i)*Con(n)
          end if
15      continue
      end if
      return
      end

************************************************************************

      subroutine SetMat(NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,
     !                  hOld,MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,
     !                  hTemp,Explic,TheTab,thSat,thR,theta)

      IMPLICIT REAL*8 (A-H,O-Z)
      logical Explic
      dimension hTab(NTab),ConTab(NTabD,NMat),CapTab(NTabD,NMat),
     !          hNew(NumNP),hOld(NumNP),MatNum(NumNP),Par(10,NMat),
     !          Con(NumNP),Cap(NumNP),ConSat(NMat),hSat(NMat),
     !          Axz(NumNP),Bxz(NumNP),Dxz(NumNP),hTemp(NumNP),
     !          TheTab(NTabD,NMat),thSat(NMat),thR(NMat),theta(NumNP)

      alh1=dlog10(-hTab(1))
      dlh =(dlog10(-hTab(NTab))-alh1)/(NTab-1)
C     write(*,770) NumNP
 770  format("numnp", i10)
      do 11 i=1,NumNP
        M=MatNum(i)
C       write(*,778) i, M
 778    format(2i10)
 779    format(e15.8)
C       write(*,779) hSat(M)
C       write(*,779) hTemp(i)
C       write(*,779) Axz(i)
        hi1=dmin1(hSat(M),hTemp(i)/Axz(i))
        hi2=dmin1(hSat(M), hNew(i)/Axz(i))
        if(Explic) hi2=hi1
        hiM=0.1*hi1+0.9*hi2
        if(hi1.ge.hSat(M).and.hi2.ge.hSat(M)) then
          Ci=ConSat(M)
        else if(hiM.ge.hTab(NTab).and.hiM.le.hTab(1)) then
          iT=int((dlog10(-hiM)-alh1)/dlh)+1
          S1=(ConTab(iT+1,M)-ConTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=ConTab(iT,M)+S1*(hiM-hTab(iT))
        else
          Ci=FK(hiM,Par(1,M))
        end if
         Con(i)=Bxz(i)*Ci
        hi1=hOld(i)/Axz(i)
        hi2=hNew(i)/Axz(i)
        if(Explic) hi2=hi1
        if(hi2.ge.hSat(M)) then
          Ci=0
          Ti=thSat(M)
        else if(hi2.ge.hTab(NTab).and.hi2.le.hTab(1)) then
          iT=int((dlog10(-hi2)-alh1)/dlh)+1
          S1=(CapTab(iT+1,M)-CapTab(iT,M))/(hTab(iT+1)-hTab(iT))
          S2=(TheTab(iT+1,M)-TheTab(iT,M))/(hTab(iT+1)-hTab(iT))
          Ci=CapTab(iT,M)+S1*(hi2-hTab(iT))
          Ti=TheTab(iT,M)+S2*(hi2-hTab(iT))
        else
          Ci=FC(hi2,Par(1,M))
          Ti=FQ(hi2,Par(1,M))
        end if
        Cap(i)=Ci*Dxz(i)/Axz(i)
        theta(i)=thR(M)+(Ti-thR(M))*Dxz(i)
11    continue
      return
      end

************************************************************************

      subroutine Solve(A,B,NumNP,MBand,MBandD)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision A,B,C
      dimension A(MBandD,NumNP),B(NumNP)

*     Reduction
      do 13 n=1,NumNP
        do 12 m=2,MBand
          if(dabs(A(m,n)).lt.1.d-30) goto 12
          C=A(m,n)/A(1,n)
          i=n+m-1
          if(i.gt.NumNP) goto 12
          j=0
          do 11 k=m,MBand
            j=j+1
            A(j,i)=A(j,i)-C*A(k,n)
11        continue
          A(m,n)=C
          B(i)=B(i)-A(m,n)*B(n)
12      continue
        B(n)=B(n)/A(1,n)
13    continue

*     Back substitution
      n=NumNP
14    do 15 k=2,MBand
        l=n+k-1
        if(l.gt.NumNP) goto 16
        B(n)=B(n)-A(k,n)*B(l)
15    continue
16    n=n-1
      if(n.gt.0) goto 14
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
