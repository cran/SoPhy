*
*  Authors:   Simunek, J., T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, martin.schlather@math.uni-goettingen.de  (2004 -- 2006)
* 
*  Copyright (C) 2002 Simunek, J., T. Vogel and M. Th. van Genuchten

* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later  version.
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


*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                      *
*     SWMS_2D  - Numerical model of two-dimensional flow and solute    *
*                transport in a variably saturated porous medium       *
*                Conjugate gradient solver for symmetric matrix        *
*                ORTHOMIN solver for asymmetric matrix                 *
*                version 1.22                                          *
*                                                                      *
*     Updated by J.Simunek (1994)                                      *
*     Based on model SWMS_2D (Simunek et al., 1992)                    *
*                                                                      *
*                                         Last modified: January, 1996 *
*                                                                      *
*||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*


C    error messages
C    0  no error
C    1 'Dimension in NSeepD is exceeded' (unused)
C    2 'Dimension in NumSPD is exceeded' (unused)
C    3 'Dimension in NDrD is exceeded'
C    4 'Dimension in NElDrD is exceeded'
C    5  jj error
C    6  iteration > breakpoint && iteration > mxLoop[IntVec(21)]
C    7  iteration > breakpoint 
C    8 'Dimension in NObsD is exceeded'
C    9 'Dimension in NumKD is exceeded'
C   10 unknown error in NodInf
C   11 'No steady state solution found'
C   12 "too many iterations (orthomin)"
C   13  iteration > mxLoop[IntVec(21)]
C   14  
C   15  
C   16  
C   17  
C   18  
C   19  

      subroutine swms2d(
     !                 IntVec,
C              1.KAT, 2.MaxIt, 3.lWat, 4.lChem, 5.ShortF,
C              6.FluxF, 7.AtmInF, 8.SeepF, 9.FreeD, 10.DrainF,
C              11.NLay, 12.NSeep, 13.NDr,  14.NObs, 
C              15.SinkF, 16.qGWLF, 17.lUpW, 18.lArtD, 19.iiprnt.incr 
C              20.start, 21.mxLoop, 22.brkpnt, 23.brk.tprint 24.ii 25.iiprnt
     !                 DblVec,
C              1.TolTh, 2.TolH, 3.hTab1, 4.hTabN, 5.dt, 6.dtMin,
C              7.dtMax, 8.dMul, 9.dMul2, 10.DrCorr, 11.rlen, 12.GWL0L,
C              13.Aqh, 14.Bqh, 15.tInit, 16.hCritS, 17.tMax, 18.epsi, 
C              19.PeCr 20.-25.cBound, 26.tPulse 27.time!

     !                  NMat, Par,
     !                  MPL, TPrint,
     !                  NSP, NP,
     !                  ND, NED, EfDim, KElDr,
     !                  NumNP, NumEl, NumBP,       
C             xR, yR, hOldR, ConcR, QR,BetaR, AxzR, BxzR, DxzR
     !                  nCodeM, d1doub,
C             ConAxx, ConAxz, ConAzz
     !                  KXR, ConAX,
     !                  KXB, Width, Node,
     !                  MaxAL, 
     !                  ChPar, KodCB, 
C             tAtm,Prec,cPrec,rSoil,rRoot,hCritA,rGWL,GWL,crt,cht
     !                  atmX,
     !                  PRL, 
C    output parameters and variables 
     !                  TlSolR, TlSolC, balaR, balaC, 
     !                  hQThF, TlSolO, atmOut, balanc, boundr,
     !                  NrErr,
     !                  intpck, dblpck
     !                  )


      parameter (MBandD=20,
     !           NSeepD=2,
     !           NumSPD=50,
     !           NDrD=2,
     !           NElDrD=8,
     !           NTabD=100,
     !           NumKD=6,
     !           NObsD=4,
     !           MNorth=4,
     !           NInt = 25, NDbl = 27)

      IMPLICIT REAL*8 (A-H,O-Z)
      integer IntVec(NInt)
      double precision DblVec(NDbl)

      integer NMat

      double precision TPrint(MPL + 1)

      integer NSP, NP
*      dimension NSP(NSeepD),NP(NSeepD,NumSPD)

      integer ND, NED
*     dimension ND(NDrD), NED(NDrD) 
      
      integer NumNP, NumEl, NumBP
      
      integer nCodeM(NumNP, 3)
      double precision d1doub(NumNP, 9)
*     dimension Kode(NumNP), MatNum(NumNP)

      integer KXR(NumEl, 6)
      double precision ConAX(NumEl, 3)
*     dimension KX(NumEl,4)
      
      integer KXB, Node
*     dimension KXB(NumBP), Node(NObsD)
      
      integer MaxAl 

      integer KodCB
      double precision ChPar
*     dimension KodCB(NumBP),
      
      double precision atmX(MaxAl, 10)
      
      double precision PRL(NumNP, 7)

      integer TlSolC, TlSolR, balaC, balaR
      double precision hQThF(6, NumNP, MPL + 2), TlSolO(TlSolR,TlSolC), 
     !   atmOut(MaxAL, 9), balanc(balaR, balaC), 
     !   boundr(NumBP, 11, balaR)
      
      integer NrErr, intpck(*) 
      double precision dblpck(*)
 
C     end of new definitions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      double precision A,B, prii
c      double precision RTime1,RTime2
      double precision A1,B1,VRV,RES,RQI,RQ,QQ,QI,RQIDOT,ECNVRG,
     !                 RCNVRG,ACNVRG
      logical lWat,lChem,SinkF,qGWLF,AtmInF,ShortF,SeepF,FluxF,
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt,lMinStep, ldebug
      integer PLevel,ALevel,TLevel

      integer ii, jj, atmkk, nPLvl, TLkk, nbalnc

      dimension prii(5),A(MBandD,NumNP),B(NumNP),Kode(NumNP),Q(NumNP),
     !  hNew(NumNP),hTemp(NumNP),hOld(NumNP),ConSat(NMat),F(NumNP),
     !  hTab(NTabD),ConTab(NTabD,NMat),CapTab(NTabD,NMat),Con(NumNP),
     !  Cap(NumNP),x(NumNP),y(NumNP),MatNum(NumNP),LayNum(NumEl),
     !  KX(NumEl,4),KXB(NumBP),Par(10,NMat),Width(NumBP),
     !  ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),SWidth(NumKD),
     !  NP(NSeepD,NumSPD),NSP(NSeepD),hSat(NMat),WatIn(NumEl),
     !  Axz(NumNP),Bxz(NumNP),Dxz(NumNP),thR(NMat),thSat(NMat),
     !  TheTab(NTabD,NMat),ThNew(NumNP),ThOld(NumNP),ListNE(NumNP),
     !  Sink(NumNP),Beta(NumNP),DS(NumNP),CumQ(NumKD),
     !  vMean(NumKD),hMean(NumKD),KodCB(NumBP),Qc(NumNP),Vx(NumNP),
     !  Vz(NumNP),ChPar(10,NMat),Dispxx(NumNP),Dispzz(NumNP),
     !  Dispxz(NumNP),cBound(6),Ac(NumNP),Fc(NumNP),SolIn(NumEl),
     !  Conc(NumNP),SMean(NumKD),ChemS(NumKD),WeTab(3,2*NumEl),
     !  Gc(NumNP),ND(NDrD),NED(NDrD),EfDim(2,NDrD),KElDr(NDrD,NElDrD),
     !  ConO(NumNP),Node(NObsD),
     !  B1(NumNP),IAD(MBandD,NumNP),IADN(NumNP),IADD(NumNP),
     !  A1(MBandD,NumNP),RES(NumNP),VRV(NumNP),RQI(NumNP,MNorth),
     !  RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNP,MNorth)


C   #### to avoid "uninitialized" error messages
      iiprnt = 0

C !! start or continue from a previous run? 
      if (IntVec(20).eq.0) goto 2  

      t = DblVec(27)
      atmkk = 0
      TLkk = 0
      nPLvl = 1
      nbalnc = 1
      SinkF = .false.
      qGWLF = .false.
      tInit = 0
      NTab = 100
      ItCum = 0
      Iter = 0
      TLevel = 1
      ALevel = 1
      PLevel = 1
      do 997 i=1,NumKD
         CumQ(i) =  0.0
         ChemS(i) = 0.0
 997  continue
      do 998 i=1,NumNP
        Sink(i) = 0.0
 998  continue
      CumQrT = 0.0
      CumQrR = 0.0
      CumQvR = 0.0
      rRoot = 0.0
      rTop = 0.0
      CumCh0 = 0.0
      CumCh1 = 0.0
      CumChR = 0.0
      dtMaxC = 1.d+30
      wCumA = 0.0
      cCumA = 0.0
      Explic = .false.
      lMinStep = .false.
      ECNVRG = 1.0d-8
      ACNVRG = 1.0d-8
      RCNVRG = 1.0d-8
      MaxItO = 200
      

* --- Reading of the input files and initial calculations --------------
      ldebug = .false.
C      ldebug = .true.

      if (ldebug) write(*,*) "BasInf"
      call BasInf (KAT,MaxIt,TolTh,TolH,lWat,lChem,AtmInF,ShortF,SeepF,
     !             FluxF,FreeD,DrainF,
     !             IntVec,NInt,DblVec,NDbl)
      if (ldebug) write(*,*) "NodInf"
      call NodInf (NumNP,NumEl,NumBP,NumKD,NObs,
     !             NObsD,Kode,Q,Conc,hNew,hOld,hTemp,x,y,MatNum,Beta,
     !             Axz,Bxz,Dxz,
     !             IntVec,NInt,nCodeM,d1doub,
     !             NrErr)
      
C      write(*,800) (hNew(iii), iii=1,42)
 800  format(5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",
     !1e14.7,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",5e14.7,"\n",
     !1e14.7,"\n")

      if (NrErr.NE.0) goto 999
      if (ldebug) write(*,*) "ElemIn"
      call ElemIn (NumEl,NumNP,KX,LayNum,ConAxx,ConAzz,ConAxz,
     !             ListNE,MBand,MBandD,lChem,lOrt,
     !             IntVec,NInt,ConAX,KXR
     !             )
      if (ldebug) write(*,*) "GeomIn"
      call GeomIn (NumKD,NumNP,NumBP,NObs,NObsD,SWidth,Width,Kode,KXB,
     !             rLen,Node,
     !             DblVec,NDbl)
      call IADMake(KX,NumNP,NumEl,MBandD,IAD,IADN,IADD)
      if (ldebug) write(*,*) "MatIn"
      call MatIn  (NMat,NLay,Par,hTab(1),hTab(NTab), 
     !             IntVec, NInt,DblVec,NDbl)
      call GenMat (NTab,NTabD,NMat,thR,hSat,Par,hTab,ConTab,CapTab,
     !             ConSat,TheTab,thSat)
      if (ldebug) write(*,*) "SetMat"
      call SetMat (NumNP,NTab,NTabD,NMat,hTab,ConTab,CapTab,hNew,hOld,
     !             MatNum,Par,Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,hTemp,
     !             Explic,TheTab,thSat,thR,ThOld)
      if(AtmInF) then
         if (ldebug) write(*,*) "AtmIn"
        call AtmIn (GWL0L,SinkF,qGWLF,tInit,tMax,Aqh,Bqh,hCritS,MaxAL,
     !              IntVec, NInt, DblVec, NDbl
     !              )
        call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,Q,
     !              NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt,lMinStep,
     !              atmX, MaxAL, atmkk
     !             )
      end if
      if (ldebug) write(*,*) "TmIn"
      call TmIn    (tInit,tMax,tAtm,tOld,dt,dtMax,dMul,dMul2,dtMin,
     !              TPrint,t,dtOpt,dtOld,AtmInF,
     !              IntVec, NInt, DblVec, NDbl, MPL)
      dtInit=dt
      if(SinkF) then
       if (ldebug) write(*,*) "sinkin"
       call SinkIn(NumEl,NumNP,KAT,KX,x,y, Beta, DblVec, NDbl)
        call SetSnk(NumNP,hNew,rRoot,Sink,PRL,Beta,rLen)
      end if     
      NSeep = IntVec(12)

      if(DrainF)
     !  call DrainIn(NDr,NDrD,NElDrD,NumEl,ND,NED,KElDr,EfDim,ConAxx,
     !               ConAxz,ConAzz, 
     !               IntVec, NInt, DblVec, NDbl, 
     !               NrErr)
      if (NrErr.NE.0) goto 999
      if(lChem) then
      if (ldebug) write(*,*) "chemin"
        call ChemIn(NMat,NumBP,cBound,epsi,tPulse,KodCB,NLevel,
     !              lUpW,lArtD,PeCr,IntVec, NInt, DblVec, NDbl)
       if(lWat)
     !  call ChInit(NumNP,NumEl,NMat,x,y,KX,MatNum,NLevel,Con,
     !              hNew,Sink,cBound(5),Vx,Vz,ConAxx,ConAzz,ConAxz,
     !              Dispxx,Dispzz,Dispxz,ChPar,ThOld,thSat,Conc,Fc,Gc,
     !              ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,Courant,KAT,
     !              lArtD,PeCr,ConO)
       call cOut  (NumNP,Conc,x,y,tInit, hQThF, MPL, nPLvl)
      else
         lUpW = .false.
         lArtD = .false.
      end if

      if (ldebug) write(*,*) "hout"
      call hOut  (hNew,x,y,NumNP,tInit, hQThF, MPL, nPLvl)
      call thOut (ThOld,x,y,NumNP,tInit, hQThF, MPL, nPLvl)
      if (ldebug) write(*,*) "subreg"
      call SubReg(NumEl,NumNP,NMat,hNew,ThOld,ThOld,x,y,MatNum,
     !            LayNum,KX,KAT,tInit,dt,NLay,0,lWat,lChem,Conc,ChPar,
     !            wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,SolIn,
     !            balanc, balaR, balaC, nbalnc)
CC      if(NObs.gt.0) call ObsNod(tInit,NumNP,NObs,NObsD,Node,hNew,ThOld,
CC     !                          Conc, TlSolO, TlSolR, TlSolC, NumKD, 
CC     !                          TLkk)

C      write(*,*)'beginning of numerical solution'
c      call getdat(i,i,iday)
c      call gettim(ihours,mins,isecs,i)
c      Rtime1=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs

      iiprnt = IntVec(19)
      goto 1

 2    continue

      ldebug = .false.
C      ldebug = .true.

      iiprnt = IntVec(25)
      call xunpck(intpck, dblpck,
     !  A, B, Q, hNew, htemp, hOld ,ConSat, F,
     !  hTab, ConTab, CapTab, Con,  Cap, x, y, LayNum,
     !  Par, ConAxx, ConAzz, ConAxz, SWidth,
     !  hSat, WatIn,  Axz, Bxz, Dxz, thR, thSat,
     !  TheTab, ThNew, ThOld, ListNE,  Sink, Beta, DS, CumQ,
     !  vMean, hMean,  Qc, Vx,  Vz, ChPar, Dispxx, Dispzz,
     !  Dispxz, cBound, Ac, Fc, SolIn,  Conc, SMean, ChemS, WeTab,
     !  Gc, ConO,
     !  B1, IAD, IADN, IADD,  A1, RES, VRV, RQI,  RQ, QQ, RQIDOT, QI,
     !
     !  SinkF, qGWLF, tInit, NTab, ItCum, Iter, TLevel, ALevel, PLevel,
     !  CumQrT, CumQrR, CumQvR, rRoot, rTop, CumCh0, CumCh1, CumChR,
     !  dtMaxC, wCumA, cCumA, Explic, lMinStep, ECNVRG, ACNVRG, RCNVRG,
     !  MaxItO, 
     !  lWat, lChem, AtmInF, ShortF, SeepF, FluxF, lUpW, FreeD,
     !  DrainF, lArtD, lOrt, 
     !  atmkk, nPLvl, TLkk, nbalnc, 
     !  t, KAT, MaxIt, TolTh, TolH, MBand, rLen, NLay, GWL0L, tMax, 
     !  Aqh, Bqh, hCritS, MaxAL, tAtm, hCritA, cPrec, cht, crt,
     !  tOld, dt, dtMax, dMul, dMul2, dtMin, dtOpt, dtOld, 
     !  dtInit, NSeep, epsi, tPulse, NLevel, PeCr, Peclet, Courant,
     !  MPL, wCumT, cCumT, wVolI, cVolI,
     !  vMeanR, hMeanR, NDr, hMeanT, hMeanG, NObs, 
     !  Kode, MatNum, KX, 
     !
     !  MBandD, NumNP, NMat, NTabD, NumEl, NumBP, NumKD, NSeepD, 
     !  NumSPD, NDrD, NElDrD, NObsD, MNorth
     !  )

* --- Beginning of time loop -------------------------------------------

 1    continue
      ii = IntVec(24)
      jj = 0
11    continue
      
C      ldebug = .true.

      if (ii.ge.iiprnt) then
         if (iiprnt.eq.IntVec(19)) 
     !      call dblepr("\ncurrent time; percent of allowed number 
     !                    of iterations",  
     !                  53, prii, 0)
         iiprnt = iiprnt + IntVec(19)
         prii(1) = t
         prii(2) = ii * 100.0 / IntVec(21)
         call dblepr(" ", 0, prii, 2)
      end if

      if (ldebug) write(*,739) ii, t, IntVec(21)
 739  format('  ii=',i7, 1e15.8,'  max loop=',i7)  

      if (ii.ge.IntVec(21)) then
         if (jj.eq.IntVec(22)) then 
            NrErr = 6
         else
            NrErr = 13
         end if
         goto 999
      else
         if (jj.eq.IntVec(22)) then
            NrErr = 7
            goto 999
         end if
      end if
    
      ii = ii + 1
      jj = jj + 1

*     Calculate water flow
      if (ldebug) write(*,*) "WatFlow"
      if(lWat.or.TLevel.eq.1)
     !  call WatFlow(NumNP,NumEl,NTab,NTabD,MBand,MBandD,NMat,
     !              NSeep,NSeepD,NumSPD,NSP,NP,NumBP,ItCum,MaxIt,Iter,
     !              Kode,KAT,t,dt,dtMin,dtOpt,dtOld,tOld,hCritA,hCritS,
     !              TolTh,TolH,rLen,Width,rTop,vMeanR,hMeanR,AtmInf,
     !              SinkF,SeepF,qGWLF,FreeD,Par,hTab,ConTab,CapTab,
     !              TheTab,hNew,hOld,hTemp,thR,thSat,ThNew,ThOld,MatNum,
     !              Con,Cap,ConSat,Axz,Bxz,Dxz,hSat,A,B,Q,F,x,y,KX,Sink,
     !              DS,Beta,ConAxx,ConAzz,ConAxz,KXB,Explic,GWL0L,Aqh,
     !              Bqh,lWat,TLevel,lOrt,DrainF,ND,NDr,NDrD,rRoot,PRL,
     !              ConO,
     !              A1,B1,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO,
     !              NrErr)
      if (ldebug) write(*,*) "WatEndFlow"
      if (NrErr.NE.0) goto 999
      if(.not.lWat.and.TLevel.eq.1) then
        if(lChem) then
      if (ldebug) write(*,*) "ChInit"
          call ChInit(NumNP,NumEl,NMat,x,y,KX,MatNum,NLevel,Con,
     !                hNew,Sink,cBound(5),Vx,Vz,ConAxx,ConAzz,ConAxz,
     !                Dispxx,Dispzz,Dispxz,ChPar,ThNew,thSat,Conc,Fc,Gc,
     !                ListNE,lUpW,WeTab,dt,dtMaxC,Peclet,Courant,KAT,
     !                lArtD,PeCr,ConO)
        else
      if (ldebug) write(*,*) "Veloc"
          call Veloc(KAT,NumNP,NumEl,hNew,x,y,KX,ListNE,Con,
     !               ConAxx,ConAzz,ConAxz,Vx,Vz)
        end if
        Iter=1
      end if

*     Calculate solute transport
      if (ldebug) write(*,*) "Solute"
      if(lChem)
     !  call Solute(NumNP,NumEl,MBand,MBandD,NMat,t,Kode,A,B,Q,
     !              hNew,hOld,F,x,y,KX,KAT,dt,DS,Sink,MatNum,Con,ConO,
     !              ConAxx,ConAzz,ConAxz,Vx,Vz,Dispxx,Dispzz,Dispxz,
     !              ChPar,ThNew,ThOld,thSat,Ac,Fc,Gc,Qc,Conc,ListNE,
     !              cBound,tPulse,NumBP,KodCB,KXB,NLevel,cPrec,crt,cht,
     !              lWat,lUpW,WeTab,epsi,CumCh0,CumCh1,CumChR,dtMaxC,
     !              Peclet,Courant,lArtD,PeCr,lOrt,
     !              A1,B1,IAD,IADN,IADD,VRV,RES,RQI,RQ,QQ,QI,
     !              RQIDOT,ECNVRG,RCNVRG,ACNVRG,MNorth,MaxItO, NrErr)
      if (NrErr.NE.0) goto 999

*     T-Level information
      if ((.not.lWat.and.TLevel.eq.1).or.
     !    (.not.ShortF.or.abs(TPrint(PLevel)-t).lt.0.001*dt)) then
         TLkk = TLkk + 1
         if (TLkk.gt.TlSolR) then
            write(*,741) TLkk, TlSolR, lWat, TLevel, PLevel, 
     !                   TPrint(PLevel), t
 741        format('TLkk =',i5,' > TlSolR =', i5,' --  lWat',l8, 
     !             '   Tlevel',i15, '   Plevel', i15, 2f15.7)
            stop
         end if 
      end if
      if (ldebug) write(*,*) "TLInf"
      call TLInf   (NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,dt,
     !              TLevel,ShortF,TPrint(PLevel),Iter,ItCum,rTop,rRoot,
     !              vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,CumQrT,
     !              CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,rLen,
     !              Peclet,Courant,wCumT,wCumA, TlSolO, TlSolR, TlSolC, 
     !              TLkk)
      if (ldebug) write(*,*) "SolInf"
      if(lChem)
     !  call SolInf(NumNP,Kode,Qc,t,dt,TLevel,ShortF,TPrint(PLevel),
     !              NumKD,SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,cCumT,
     !              lWat, TlSolO, TlSolR, TlSolC, TLkk)
      if (ldebug) write(*,*) "ObsNod"
      if(NObs.gt.0.and.
     !  (.not.ShortF.or.abs(TPrint(PLevel)-t).lt.0.001*dt))
     !     call ObsNod(t,NumNP,NObs,NObsD,Node,hNew,ThNew,Conc,
     !                 TlSolO, TlSolR, TlSolC, NumKD, TLkk)

*     P-Level information
      if(abs(TPrint(PLevel)-t).lt.0.001*dt.or.
     !     (.not.lWat.and..not.lChem)) then   
        if (lWat.or.lChem.or.(PLevel.eq.1)) then
           nPLvl = nPLvl + 1
           if  (nPLvl.gt.MPL + 2) then
             NrErr = 14
             goto 999
           end if
        end if
        if(lWat.or.(.not.lWat.and.PLevel.eq.1)) then
          if (ldebug) write(*,*) "hout"
          call hOut (hNew,x,y,NumNP,t,hQThF, MPL, nPLvl)
          if (ldebug) write(*,*) "thout"
          call thOut(ThNew,x,y,NumNP,t, hQThF, MPL, nPLvl)
          if(FluxF) then
            if (ldebug) write(*,*) "flxout, qout"
            if(.not.lChem)
     !        call Veloc(KAT,NumNP,NumEl,hNew,x,y,KX,ListNE,Con,
     !                   ConAxx,ConAzz,ConAxz,Vx,Vz)
            call FlxOut(Vx,Vz,x,y,NumNP,t, hQThF, MPL, nPLvl)
            call QOut  (Q,x,y,NumNP,t, hQThF, MPL, nPLvl)
          end if
        end if
        if (ldebug) write(*,*) "subreg"
        nbalnc = nbalnc + 1
        if (nbalnc > balaR) then
           NrErr = 15
           goto 999
        end if
        call SubReg(NumEl,NumNP,NMat,hNew,ThOld,ThNew,x,y,MatNum,
     !              LayNum,KX,KAT,t,dt,NLay,PLevel,lWat,lChem,Conc,
     !              ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,
     !              SolIn, balanc, balaR, balaC, nbalnc)
        if (ldebug) write(*,*) "bouout"
        call BouOut(NumNP,NumBP,t,hNew,ThNew,Q,Width,KXB,Kode,x,y,Conc,
     !              boundr, balaR, nbalnc)
        if(lChem) then
          if (ldebug) write(*,*) "cOut"
          call cOut(NumNP,Conc,x,y,t,hQThF, MPL, nPLvl)
        end if
        PLevel=PLevel+1
      else 
         if (jj.eq.IntVec(22)) then
            nnn = nPLvl + 1
            if  (nnn.gt.MPL + 2) then
               NrErr = 5
               goto 999
            end if
            call hOut (hNew,x,y,NumNP,t,hQThF, MPL, nnn)
            call thOut(ThNew,x,y,NumNP,t, hQThF, MPL, nnn)
            if(FluxF) then
               if(.not.lChem)
     !              call Veloc(KAT,NumNP,NumEl,hNew,x,y,KX,ListNE,Con,
     !              ConAxx,ConAzz,ConAxz,Vx,Vz)
               call FlxOut(Vx,Vz,x,y,NumNP,t, hQThF, MPL, nnn)
               call QOut  (Q,x,y,NumNP,t, hQThF, MPL, nnn)
            endif
            if (lChem) call cOut(NumNP, Conc, x, y, t, hQThF, MPL, nnn)
           IntVec(23) = nnn
         end if
      end if
      if (ldebug) write(*,*) "end P-Level"
       
*     A-level information
      if(abs(t-tAtm).le.0.001*dt.and.AtmInF) then
        if (ldebug) write(*,*) "alinf"
        if(lWat.and.(ALevel.le.MaxAL))
     !    call ALInf (t,CumQ,hMeanT,hMeanR,hMeanG,ALevel,CumQrT,CumQrR,
     !                CumQvR,NumKD, atmOut, MaxAL)
        if(ALevel.lt.MaxAL) then
          if (ldebug) write(*,*) "setatm"
          call SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,hNew,
     !                Q,NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt,lMinStep,
     !                  atmX, MaxAL, atmkk)
          ALevel=ALevel+1
        end if
      end if


*     Root extraction
*      if(SinkF)
*     !  call SetSnk(NumNP,hNew,rRoot,Sink,PRL,Beta,rLen)
*     Time governing
      if(abs(t-tMax).le.0.001*dt.or.(.not.lWat.and..not.lChem)) then
         if (ldebug) write(*,*) "the end"
c        call getdat(i,i,iday)
c        call gettim(ihours,mins,isecs,i)
c        Rtime2=iday*24.*60.*60.+ihours*60.*60.+mins*60.+isecs
c        write(70,*)
c        write(70,*) 'Real time [sec]',Rtime2-RTime1
c        write( *,*) 'Real time [sec]',Rtime2-RTime1
         goto 999
C         stop
      end if
      tOld=t
      dtOld=dt
      if (ldebug)      write(*,*) "tmcont"
      call TmCont(dt,dtMax,dtOpt,dMul,dMul2,dtMin,Iter,TPrint(PLevel),
     !            tAtm,t,tMax,dtMaxC,lMinStep,dtInit)
      TLevel=TLevel+1
      t=t+dt
      goto 11

C     final settings
 999  continue
      DblVec(27) = t
      IntVec(24) = ii
      IntVec(25) = iiprnt
C NP, NSP, ND, NED, KXB, Node, KodCB, PRL, TPrint
      call xpack(intpck, dblpck,
     !  A, B, Q, hNew, htemp, hOld ,ConSat, F,
     !  hTab, ConTab, CapTab, Con,  Cap, x, y, LayNum,
     !  Par,  ConAxx, ConAzz, ConAxz, SWidth,
     !  hSat, WatIn,  Axz, Bxz, Dxz, thR, thSat,
     !  TheTab, ThNew, ThOld, ListNE,  Sink, Beta, DS, CumQ,
     !  vMean, hMean,  Qc, Vx,  Vz, ChPar, Dispxx, Dispzz,
     !  Dispxz, cBound, Ac, Fc, SolIn,  Conc, SMean, ChemS, WeTab,
     !  Gc, ConO,
     !  B1, IAD, IADN, IADD,  A1, RES, VRV, RQI,  RQ, QQ, RQIDOT, QI,
     !
     !  SinkF, qGWLF, tInit, NTab, ItCum, Iter, TLevel, ALevel, PLevel,
     !  CumQrT, CumQrR, CumQvR, rRoot, rTop, CumCh0, CumCh1, CumChR,
     !  dtMaxC, wCumA, cCumA, Explic, lMinStep, ECNVRG, ACNVRG, RCNVRG,
     !  MaxItO, 
     !  lWat, lChem, AtmInF, ShortF, SeepF, FluxF, lUpW, FreeD,
     !  DrainF, lArtD, lOrt, 
     !  atmkk, nPLvl, TLkk, nbalnc, 
     !  t, KAT, MaxIt, TolTh, TolH, MBand, rLen, NLay, GWL0L, tMax, 
     !  Aqh, Bqh, hCritS, MaxAL, tAtm, hCritA, cPrec, cht, crt,
     !  tOld, dt, dtMax, dMul, dMul2, dtMin, dtOpt, dtOld, 
     !  dtInit, NSeep, epsi, tPulse, NLevel, PeCr, Peclet, Courant,
     !  MPL, wCumT, cCumT, wVolI, cVolI,
     !  vMeanR, hMeanR, NDr, hMeanT, hMeanG, NObs,
     !  Kode, MatNum, KX, 
     !
     !  MBandD, NumNP, NMat, NTabD, NumEl, NumBP, NumKD, NSeepD, 
     !  NumSPD, NDrD, NElDrD, NObsD, MNorth
     !  )
      return

* --- end of time loop -------------------------------------------------

      end
 
   
  
     

*|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
