*                                                                     *
*     an analytic solution of water transport                         *
*                                                                     *
*     Martin Schlather and Bernd Huwe (2003)                          *
*                                                                     *
*     The use of the language interface of R: two examples for        *
*     modelling water flux and solute transport                       *
*                                                                     *
*                            Computers & Geosciences                  *
*                                                                     *
* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      subroutine xpack(intpck, dblpck,
     !  A, B, Q, hNew, htemp, hOld ,ConSat, F,
     !  hTab, ConTab, CapTab, Con, Cap, x, y, LayNum,
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
     !               )
   
      IMPLICIT REAL*8 (A-H,O-Z)
      
      integer intpck(*)
      double precision dblpck(*)

      dimension A(MBandD * NumNP),B(NumNP), Q(NumNP),
     !  hNew(NumNP),hTemp(NumNP),hOld(NumNP),ConSat(NMat),F(NumNP),
     !  hTab(NTabD),ConTab(NTabD * NMat),CapTab(NTabD*NMat),Con(NumNP),
     !  Cap(NumNP),x(NumNP),y(NumNP),LayNum(NumEl),
     !  Par(10 * NMat),
     !  ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),SWidth(NumKD),
     !  hSat(NMat),WatIn(NumEl),
     !  Axz(NumNP),Bxz(NumNP),Dxz(NumNP),thR(NMat),thSat(NMat),
     !  TheTab(NTabD * NMat),ThNew(NumNP),ThOld(NumNP),ListNE(NumNP),
     !  Sink(NumNP),Beta(NumNP),DS(NumNP),CumQ(NumKD),
     !  vMean(NumKD),hMean(NumKD),Qc(NumNP),Vx(NumNP),
     !  Vz(NumNP),ChPar(10 * NMat),Dispxx(NumNP),Dispzz(NumNP),
     !  Dispxz(NumNP),cBound(6),Ac(NumNP),Fc(NumNP),SolIn(NumEl),
     !  Conc(NumNP),SMean(NumKD),ChemS(NumKD),WeTab(3 * 2*NumEl),
     !  Gc(NumNP),
     !  ConO(NumNP),
     !  B1(NumNP),IAD(MBandD * NumNP),IADN(NumNP),IADD(NumNP),
     !  A1(MBandD * NumNP),RES(NumNP),VRV(NumNP),RQI(NumNP * MNorth),
     !  RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNP * MNorth),
     !  Kode(NumNP),MatNum(NumNP),KX(NumEl * 4)
      logical lWat,lChem,SinkF,qGWLF,AtmInF,ShortF,SeepF,FluxF,
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt,lMinStep, ldebug
      integer PLevel,ALevel,TLevel, atmkk, nPLvl, TLkk, nbalnc

      idbl = 0
      iint = 0
      ndbl = MBandD * NumNP
      call cpdbl(A, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(B, dblpck, ndbl, idbl)
      call cpdbl(Q, dblpck, ndbl, idbl)
      call cpdbl(hNew, dblpck, ndbl, idbl)
      call cpdbl(hTemp, dblpck, ndbl, idbl)
      call cpdbl(hOld, dblpck, ndbl, idbl)
      ndbl = NMat
      call cpdbl(ConSat, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(F, dblpck, ndbl, idbl)
      ndbl = NTabD
      call cpdbl(hTab, dblpck, ndbl, idbl)
      ndbl = NTabD * NMat
      call cpdbl(ConTab, dblpck, ndbl, idbl)
      call cpdbl(CapTab, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(Con, dblpck, ndbl, idbl)
      call cpdbl(Cap, dblpck, ndbl, idbl)
      call cpdbl(x, dblpck, ndbl, idbl)
      call cpdbl(y, dblpck, ndbl, idbl)
      ndbl = 10 * NMat
      call cpdbl(Par, dblpck, ndbl, idbl)
      ndbl = NumEl
      call cpdbl(ConAxx, dblpck, ndbl, idbl)
      call cpdbl(ConAzz, dblpck, ndbl, idbl)
      call cpdbl(ConAxz, dblpck, ndbl, idbl)
      ndbl = NumKD
      call cpdbl(SWidth, dblpck, ndbl, idbl)
      ndbl = NMat
      call cpdbl(hSat, dblpck, ndbl, idbl)
      ndbl = NumEl
      call cpdbl(WatIn, dblpck, ndbl, idbl)
      ndbl = NumBP
      call cpdbl(Axz, dblpck, ndbl, idbl)
      call cpdbl(Bxz, dblpck, ndbl, idbl)
      call cpdbl(Dxz, dblpck, ndbl, idbl)
      ndbl = NMat
      call cpdbl(thR, dblpck, ndbl, idbl)
      call cpdbl(thSat, dblpck, ndbl, idbl)
      ndbl = NTabD * NMat
      call cpdbl(TheTab, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(ThNew, dblpck, ndbl, idbl)
      call cpdbl(ThOld, dblpck, ndbl, idbl)
      call cpdbl(Sink, dblpck, ndbl, idbl)
      call cpdbl(Beta, dblpck, ndbl, idbl)
      call cpdbl(DS, dblpck, ndbl, idbl)
      ndbl = NumKD
      call cpdbl(CumQ, dblpck, ndbl, idbl)
      call cpdbl(vMean, dblpck, ndbl, idbl)
      call cpdbl(hMean, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(Qc, dblpck, ndbl, idbl)
      call cpdbl(Vx, dblpck, ndbl, idbl)
      call cpdbl(Vz, dblpck, ndbl, idbl)
      ndbl = 10 * NMat
      call cpdbl(ChPar, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(Dispxx, dblpck, ndbl, idbl)
      call cpdbl(Dispzz, dblpck, ndbl, idbl)
      call cpdbl(Dispxz, dblpck, ndbl, idbl)
      ndbl = 6
      call cpdbl(cBound, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(Ac, dblpck, ndbl, idbl)
      call cpdbl(Fc, dblpck, ndbl, idbl)
      ndbl = NumEl
      call cpdbl(SolIn, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(Conc, dblpck, ndbl, idbl)
      ndbl = NumKD
      call cpdbl(SMean, dblpck, ndbl, idbl)
      call cpdbl(ChemS, dblpck, ndbl, idbl)
      ndbl = 3 * 2 * NumEl
      call cpdbl(WeTab, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(Gc, dblpck, ndbl, idbl)
      call cpdbl(ConO, dblpck, ndbl, idbl)
      call cpdbl(B1, dblpck, ndbl, idbl)
      ndbl = MBandD * NumNP
      call cpdbl(A1, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(RES, dblpck, ndbl, idbl)
      call cpdbl(VRV, dblpck, ndbl, idbl)
      ndbl = NumNP * MNorth
      call cpdbl(RQI, dblpck, ndbl, idbl)
      ndbl = NumNP
      call cpdbl(RQ, dblpck, ndbl, idbl)
      call cpdbl(QQ, dblpck, ndbl, idbl)
      ndbl = MNorth
      call cpdbl(RQIDOT, dblpck, ndbl, idbl)
      ndbl = NumNP * MNorth 
      call cpdbl(QI, dblpck, ndbl, idbl)
      ndbl = 1
      call cpdbl(tInit, dblpck, ndbl, idbl)
      call cpdbl(CumQrT, dblpck, ndbl, idbl)
      call cpdbl(CumQrR, dblpck, ndbl, idbl)
      call cpdbl(CumQvR, dblpck, ndbl, idbl)
      call cpdbl(rRoot, dblpck, ndbl, idbl)
      call cpdbl(rTop, dblpck, ndbl, idbl)
      call cpdbl(CumCh0, dblpck, ndbl, idbl)
      call cpdbl(CumCh1, dblpck, ndbl, idbl)
      call cpdbl(CumChR, dblpck, ndbl, idbl)
      call cpdbl(dtMaxC, dblpck, ndbl, idbl)
      call cpdbl(wCumA, dblpck, ndbl, idbl)
      call cpdbl(cCumA, dblpck, ndbl, idbl)
      call cpdbl(ECNVRG, dblpck, ndbl, idbl)
      call cpdbl(ACNVRG, dblpck, ndbl, idbl)
      call cpdbl(RCNVRG, dblpck, ndbl, idbl)
      call cpdbl(t, dblpck, ndbl, idbl)
      call cpdbl(TolTh, dblpck, ndbl, idbl)
      call cpdbl(TolH, dblpck, ndbl, idbl)
      call cpdbl(rLen, dblpck, ndbl, idbl)
      call cpdbl(GWL0L, dblpck, ndbl, idbl)
      call cpdbl(tMax, dblpck, ndbl, idbl)
      call cpdbl(Aqh, dblpck, ndbl, idbl)
      call cpdbl(Bqh, dblpck, ndbl, idbl)
      call cpdbl(hCritS, dblpck, ndbl, idbl)
      call cpdbl(tAtm, dblpck, ndbl, idbl)
      call cpdbl(hCritA, dblpck, ndbl, idbl)
      call cpdbl(cPrec, dblpck, ndbl, idbl)
      call cpdbl(cht, dblpck, ndbl, idbl)
      call cpdbl(crt, dblpck, ndbl, idbl)
      call cpdbl(tOld, dblpck, ndbl, idbl)
      call cpdbl(dt, dblpck, ndbl, idbl)
      call cpdbl(dtMax, dblpck, ndbl, idbl)
      call cpdbl(dMul, dblpck, ndbl, idbl)
      call cpdbl(dMul2, dblpck, ndbl, idbl)
      call cpdbl(dtMin, dblpck, ndbl, idbl)
      call cpdbl(dtOpt, dblpck, ndbl, idbl)
      call cpdbl(dtOld, dblpck, ndbl, idbl)
      call cpdbl(dtInit, dblpck, ndbl, idbl)
      call cpdbl(epsi, dblpck, ndbl, idbl)
      call cpdbl(tPulse, dblpck, ndbl, idbl)
      call cpdbl(PeCr, dblpck, ndbl, idbl)
      call cpdbl(Peclet, dblpck, ndbl, idbl)
      call cpdbl(Courant, dblpck, ndbl, idbl)
      call cpdbl(wCumT, dblpck, ndbl, idbl)
      call cpdbl(cCumT, dblpck, ndbl, idbl)
      call cpdbl(wVolI, dblpck, ndbl, idbl)
      call cpdbl(cVolI, dblpck, ndbl, idbl)
      call cpdbl(vMeanR, dblpck, ndbl, idbl)
      call cpdbl(hMeanR, dblpck, ndbl, idbl)
      call cpdbl(hMeanT, dblpck, ndbl, idbl)
      call cpdbl(hMeanG, dblpck, ndbl, idbl)
      xxx = 999999.25
      call cpdbl(xxx, dblpck, ndbl, idbl)
C      write(*, 777) idbl, dblpck(idbl)
 777  format(i10, f15.7)
C      stop
C      call cpdbl(, dblpck, ndbl, idbl)

CCCCCCCCCCCCCCCCCCCCCCCCC
      nint = NumEl
      call cpint(LayNum, intpck, nint, iint)
      nint = NumNP
      call cpint(ListNE, intpck, nint, iint)
      nint = MBandD * NumNP
      call cpint(IAD, intpck, nint, iint)
      nint = NumNP
      call cpint(IADN, intpck, nint, iint)
      call cpint(IADD, intpck, nint, iint)
      nint = 1
      call cpint(NTab, intpck, nint, iint)
      call cpint(ItCum, intpck, nint, iint)
      call cpint(Iter, intpck, nint, iint)
      call cpint(TLevel, intpck, nint, iint)
      call cpint(ALevel, intpck, nint, iint)
      call cpint(PLevel, intpck, nint, iint)
      call cpint(MaxItO, intpck, nint, iint)
      call cpint(atmkk, intpck, nint, iint)
      call cpint(nPLvl, intpck, nint, iint)
      call cpint(TLkk, intpck, nint, iint)
      call cpint(nbalnc, intpck, nint, iint)
      call cpint(KAT, intpck, nint, iint)
      call cpint(MaxIt, intpck, nint, iint)
      call cpint(MBand, intpck, nint, iint)
      call cpint(NLay, intpck, nint, iint)
      call cpint(MaxAL, intpck, nint, iint)
      call cpint(NSeep, intpck, nint, iint)
      call cpint(NLevel, intpck, nint, iint)
      call cpint(MPL, intpck, nint, iint)
      call cpint(NDr, intpck, nint, iint)
      call cpint(NObs, intpck, nint, iint)
      nint = NumNP
      call cpint(Kode, intpck, nint, iint)
      call cpint(MatNum, intpck, nint, iint)
      nint = NumEl * 4
      call cpint(KX, intpck, nint, iint)
C      call cpint(, intpck, nint, iint)

CCCCCCCCCCCCCCCCCCCCCCCCC
      nint = 1
      call cplog(SinkF, intpck, nint, iint)
      call cplog(qGWLF, intpck, nint, iint)
      call cplog(Explic, intpck, nint, iint)
      call cplog(lMinStep, intpck, nint, iint)
      call cplog(lWat, intpck, nint, iint)
      call cplog(lChem, intpck, nint, iint)
      call cplog(AtmInF, intpck, nint, iint)
      call cplog(ShortF, intpck, nint, iint)
      call cplog(SeepF, intpck, nint, iint)
      call cplog(FluxF, intpck, nint, iint)
      call cplog(lUpW, intpck, nint, iint)
      call cplog(FreeD, intpck, nint, iint)
      call cplog(DrainF, intpck, nint, iint)
      call cplog(lArtD, intpck, nint, iint)
      call cplog(lOrt, intpck, nint, iint)
C      call cplog(, intpck, nint, iint)

      nxxx = 999999
      call cpint(nxxx, intpck, nint, iint)
      return
      end


      subroutine cpdbl(from, to, nfrom, ito)
      double precision from(*), to(*)
      integer nfrom, ito
       do 11 i=1,nfrom
        to(ito + i) = from(i)
 11   continue
      ito = ito + nfrom 
      return
      end


      subroutine cpint(from, to, nfrom, ito)
      integer from(*), to(*)
      integer nfrom, ito
       do 11 i=1,nfrom
        to(ito + i) = from(i)
 11   continue
      ito = ito + nfrom 
      return
      end


      subroutine cplog(from, to, nfrom, ito)
      logical from(*)
      integer to(*)
      integer nfrom, ito
      do 11 i=1,nfrom
         if (from(i)) then 
            to(ito + i) = 1
         else
            to(ito + i) = 0
         end if
 11   continue
      ito = ito + nfrom 
      return
      end


      subroutine redbl(from, to, nfrom, ito)
      double precision from(*), to(*)
      integer nfrom, ito
       do 11 i=1,nfrom
        from(i) = to(ito + i) 
 11   continue
      ito = ito + nfrom 
      return
      end


      subroutine reint(from, to, nfrom, ito)
      integer from(*), to(*)
      integer nfrom, ito
       do 11 i=1,nfrom
         from(i) = to(ito + i)
 11   continue
      ito = ito + nfrom 
      return
      end


      subroutine relog(from, to, nfrom, ito)
      logical from(*)
      integer to(*)
      integer nfrom, ito
      do 11 i=1,nfrom
         from(i) =  to(ito + i).eq.1
 11   continue
      ito = ito + nfrom 
      return
      end


      subroutine xunpck(intpck, dblpck, 
     !  A, B, Q, hNew, htemp, hOld ,ConSat, F,
     !  hTab, ConTab, CapTab, Con, Cap, x, y, LayNum,
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
     !               )
   
      IMPLICIT REAL*8 (A-H,O-Z)
      
      integer intpck(*)
      double precision dblpck(*)

      dimension A(MBandD * NumNP),B(NumNP), Q(NumNP),
     !  hNew(NumNP),hTemp(NumNP),hOld(NumNP),ConSat(NMat),F(NumNP),
     !  hTab(NTabD),ConTab(NTabD * NMat),CapTab(NTabD*NMat),Con(NumNP),
     !  Cap(NumNP),x(NumNP),y(NumNP),LayNum(NumEl),
     !  Par(10 * NMat),
     !  ConAxx(NumEl),ConAzz(NumEl),ConAxz(NumEl),SWidth(NumKD),
     !  hSat(NMat),WatIn(NumEl),
     !  Axz(NumNP),Bxz(NumNP),Dxz(NumNP),thR(NMat),thSat(NMat),
     !  TheTab(NTabD * NMat),ThNew(NumNP),ThOld(NumNP),ListNE(NumNP),
     !  Sink(NumNP),Beta(NumNP),DS(NumNP),CumQ(NumKD),
     !  vMean(NumKD),hMean(NumKD),Qc(NumNP),Vx(NumNP),
     !  Vz(NumNP),ChPar(10 * NMat),Dispxx(NumNP),Dispzz(NumNP),
     !  Dispxz(NumNP),cBound(6),Ac(NumNP),Fc(NumNP),SolIn(NumEl),
     !  Conc(NumNP),SMean(NumKD),ChemS(NumKD),WeTab(3 * 2*NumEl),
     !  Gc(NumNP),
     !  ConO(NumNP),
     !  B1(NumNP),IAD(MBandD * NumNP),IADN(NumNP),IADD(NumNP),
     !  A1(MBandD * NumNP),RES(NumNP),VRV(NumNP),RQI(NumNP * MNorth),
     !  RQ(NumNP),QQ(NumNP),RQIDOT(MNorth),QI(NumNP * MNorth),
     !  Kode(NumNP),MatNum(NumNP),KX(NumEl*4)

      logical lWat,lChem,SinkF,qGWLF,AtmInF,ShortF,SeepF,FluxF,
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt,lMinStep, ldebug
      integer PLevel,ALevel,TLevel, atmkk, nPLvl, TLkk, nbalnc

      idbl = 0
      iint = 0
      ndbl = MBandD * NumNP
      call redbl(A, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(B, dblpck, ndbl, idbl)
      call redbl(Q, dblpck, ndbl, idbl)
      call redbl(hNew, dblpck, ndbl, idbl)
      call redbl(hTemp, dblpck, ndbl, idbl)
      call redbl(hOld, dblpck, ndbl, idbl)
      ndbl = NMat
      call redbl(ConSat, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(F, dblpck, ndbl, idbl)
      ndbl = NTabD
      call redbl(hTab, dblpck, ndbl, idbl)
      ndbl = NTabD * NMat
      call redbl(ConTab, dblpck, ndbl, idbl)
      call redbl(CapTab, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(Con, dblpck, ndbl, idbl)
      call redbl(Cap, dblpck, ndbl, idbl)
      call redbl(x, dblpck, ndbl, idbl)
      call redbl(y, dblpck, ndbl, idbl)
      ndbl = 10 * NMat
      call redbl(Par, dblpck, ndbl, idbl)
      ndbl = NumEl
      call redbl(ConAxx, dblpck, ndbl, idbl)
      call redbl(ConAzz, dblpck, ndbl, idbl)
      call redbl(ConAxz, dblpck, ndbl, idbl)
      ndbl = NumKD
      call redbl(SWidth, dblpck, ndbl, idbl)
      ndbl = NMat
      call redbl(hSat, dblpck, ndbl, idbl)
      ndbl = NumEl
      call redbl(WatIn, dblpck, ndbl, idbl)
      ndbl = NumBP
      call redbl(Axz, dblpck, ndbl, idbl)
      call redbl(Bxz, dblpck, ndbl, idbl)
      call redbl(Dxz, dblpck, ndbl, idbl)
      ndbl = NMat
      call redbl(thR, dblpck, ndbl, idbl)
      call redbl(thSat, dblpck, ndbl, idbl)
      ndbl = NTabD * NMat
      call redbl(TheTab, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(ThNew, dblpck, ndbl, idbl)
      call redbl(ThOld, dblpck, ndbl, idbl)
      call redbl(Sink, dblpck, ndbl, idbl)
      call redbl(Beta, dblpck, ndbl, idbl)
      call redbl(DS, dblpck, ndbl, idbl)
      ndbl = NumKD
      call redbl(CumQ, dblpck, ndbl, idbl)
      call redbl(vMean, dblpck, ndbl, idbl)
      call redbl(hMean, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(Qc, dblpck, ndbl, idbl)
      call redbl(Vx, dblpck, ndbl, idbl)
      call redbl(Vz, dblpck, ndbl, idbl)
      ndbl = 10 * NMat
      call redbl(ChPar, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(Dispxx, dblpck, ndbl, idbl)
      call redbl(Dispzz, dblpck, ndbl, idbl)
      call redbl(Dispxz, dblpck, ndbl, idbl)
      ndbl = 6
      call redbl(cBound, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(Ac, dblpck, ndbl, idbl)
      call redbl(Fc, dblpck, ndbl, idbl)
      ndbl = NumEl
      call redbl(SolIn, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(Conc, dblpck, ndbl, idbl)
      ndbl = NumKD
      call redbl(SMean, dblpck, ndbl, idbl)
      call redbl(ChemS, dblpck, ndbl, idbl)
      ndbl = 3 * 2 * NumEl
      call redbl(WeTab, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(Gc, dblpck, ndbl, idbl)
      call redbl(ConO, dblpck, ndbl, idbl)
      call redbl(B1, dblpck, ndbl, idbl)
      ndbl = MBandD * NumNP
      call redbl(A1, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(RES, dblpck, ndbl, idbl)
      call redbl(VRV, dblpck, ndbl, idbl)
      ndbl = NumNP * MNorth
      call redbl(RQI, dblpck, ndbl, idbl)
      ndbl = NumNP
      call redbl(RQ, dblpck, ndbl, idbl)
      call redbl(QQ, dblpck, ndbl, idbl)
      ndbl = MNorth
      call redbl(RQIDOT, dblpck, ndbl, idbl)
      ndbl = NumNP * MNorth 
      call redbl(QI, dblpck, ndbl, idbl)
      ndbl = 1
      call redbl(tInit, dblpck, ndbl, idbl)
      call redbl(CumQrT, dblpck, ndbl, idbl)
      call redbl(CumQrR, dblpck, ndbl, idbl)
      call redbl(CumQvR, dblpck, ndbl, idbl)
      call redbl(rRoot, dblpck, ndbl, idbl)
      call redbl(rTop, dblpck, ndbl, idbl)
      call redbl(CumCh0, dblpck, ndbl, idbl)
      call redbl(CumCh1, dblpck, ndbl, idbl)
      call redbl(CumChR, dblpck, ndbl, idbl)
      call redbl(dtMaxC, dblpck, ndbl, idbl)
      call redbl(wCumA, dblpck, ndbl, idbl)
      call redbl(cCumA, dblpck, ndbl, idbl)
      call redbl(ECNVRG, dblpck, ndbl, idbl)
      call redbl(ACNVRG, dblpck, ndbl, idbl)
      call redbl(RCNVRG, dblpck, ndbl, idbl)
      call redbl(t, dblpck, ndbl, idbl)
      call redbl(TolTh, dblpck, ndbl, idbl)
      call redbl(TolH, dblpck, ndbl, idbl)
      call redbl(rLen, dblpck, ndbl, idbl)
      call redbl(GWL0L, dblpck, ndbl, idbl)
      call redbl(tMax, dblpck, ndbl, idbl)
      call redbl(Aqh, dblpck, ndbl, idbl)
      call redbl(Bqh, dblpck, ndbl, idbl)
      call redbl(hCritS, dblpck, ndbl, idbl)
      call redbl(tAtm, dblpck, ndbl, idbl)
      call redbl(hCritA, dblpck, ndbl, idbl)
      call redbl(cPrec, dblpck, ndbl, idbl)
      call redbl(cht, dblpck, ndbl, idbl)
      call redbl(crt, dblpck, ndbl, idbl)
      call redbl(tOld, dblpck, ndbl, idbl)
      call redbl(dt, dblpck, ndbl, idbl)
      call redbl(dtMax, dblpck, ndbl, idbl)
      call redbl(dMul, dblpck, ndbl, idbl)
      call redbl(dMul2, dblpck, ndbl, idbl)
      call redbl(dtMin, dblpck, ndbl, idbl)
      call redbl(dtOpt, dblpck, ndbl, idbl)
      call redbl(dtOld, dblpck, ndbl, idbl)
      call redbl(dtInit, dblpck, ndbl, idbl)
      call redbl(epsi, dblpck, ndbl, idbl)
      call redbl(tPulse, dblpck, ndbl, idbl)
      call redbl(PeCr, dblpck, ndbl, idbl)
      call redbl(Peclet, dblpck, ndbl, idbl)
      call redbl(Courant, dblpck, ndbl, idbl)
      call redbl(wCumT, dblpck, ndbl, idbl)
      call redbl(cCumT, dblpck, ndbl, idbl)
      call redbl(wVolI, dblpck, ndbl, idbl)
      call redbl(cVolI, dblpck, ndbl, idbl)
      call redbl(vMeanR, dblpck, ndbl, idbl)
      call redbl(hMeanR, dblpck, ndbl, idbl)
      call redbl(hMeanT, dblpck, ndbl, idbl)
      call redbl(hMeanG, dblpck, ndbl, idbl)

      call redbl(xxx, dblpck, ndbl, idbl)  
      if (xxx.ne.999999.25) then
         write(*,*) "xxx.ne.999999.25"
         stop
      end if
C      stop
C      call redbl(, dblpck, ndbl, idbl)

CCCCCCCCCCCCCCCCCCCCCCCCC
      nint = NumEl
      call reint(LayNum, intpck, nint, iint)
      nint = NumNP
      call reint(ListNE, intpck, nint, iint)
      nint = MBandD * NumNP
      call reint(IAD, intpck, nint, iint)
      nint = NumNP
      call reint(IADN, intpck, nint, iint)
      call reint(IADD, intpck, nint, iint)
      nint = 1
      call reint(NTab, intpck, nint, iint)
      call reint(ItCum, intpck, nint, iint)
      call reint(Iter, intpck, nint, iint)
      call reint(TLevel, intpck, nint, iint)
      call reint(ALevel, intpck, nint, iint)
      call reint(PLevel, intpck, nint, iint)
      call reint(MaxItO, intpck, nint, iint)
      call reint(atmkk, intpck, nint, iint)
      call reint(nPLvl, intpck, nint, iint)
      call reint(TLkk, intpck, nint, iint)
      call reint(nbalnc, intpck, nint, iint)
      call reint(KAT, intpck, nint, iint)
      call reint(MaxIt, intpck, nint, iint)
      call reint(MBand, intpck, nint, iint)
      call reint(NLay, intpck, nint, iint)
      call reint(MaxAL, intpck, nint, iint)
      call reint(NSeep, intpck, nint, iint)
      call reint(NLevel, intpck, nint, iint)
      call reint(MPL, intpck, nint, iint)
      call reint(NDr, intpck, nint, iint)
      call reint(NObs, intpck, nint, iint)
      nint = NumNP
      call reint(Kode, intpck, nint, iint)
      call reint(MatNum, intpck, nint, iint)
      nint = NumEl * 4
      call reint(KX, intpck, nint, iint)
C      call reint(, intpck, nint, iint)

CCCCCCCCCCCCCCCCCCCCCCCCC
      nint = 1
      call relog(SinkF, intpck, nint, iint)
      call relog(qGWLF, intpck, nint, iint)
      call relog(Explic, intpck, nint, iint)
      call relog(lMinStep, intpck, nint, iint)
      call relog(lWat, intpck, nint, iint)
      call relog(lChem, intpck, nint, iint)
      call relog(AtmInF, intpck, nint, iint)
      call relog(ShortF, intpck, nint, iint)
      call relog(SeepF, intpck, nint, iint)
      call relog(FluxF, intpck, nint, iint)
      call relog(lUpW, intpck, nint, iint)
      call relog(FreeD, intpck, nint, iint)
      call relog(DrainF, intpck, nint, iint)
      call relog(lArtD, intpck, nint, iint)
      call relog(lOrt, intpck, nint, iint)
C      call relog(, intpck, nint, iint)


      call reint(nxxx, intpck, nint, iint)
      if (nxxx.ne.999999) then
         write(*,*) "nxxx.ne.999999"
         stop
      end if

      return
      end









