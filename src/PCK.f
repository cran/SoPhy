* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                     *
*     an analytic solution of water transport                         *
*                                                                     *
*     Martin Schlather and Bernd Huwe (2002, 2003)                    *
*                                                                     *
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
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt,lMinStep
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
C      ndbl = 1
      call Ecpdbl(tInit, dblpck, idbl)
      call Ecpdbl(CumQrT, dblpck, idbl)
      call Ecpdbl(CumQrR, dblpck, idbl)
      call Ecpdbl(CumQvR, dblpck, idbl)
      call Ecpdbl(rRoot, dblpck, idbl)
      call Ecpdbl(rTop, dblpck, idbl)
      call Ecpdbl(CumCh0, dblpck, idbl)
      call Ecpdbl(CumCh1, dblpck, idbl)
      call Ecpdbl(CumChR, dblpck, idbl)
      call Ecpdbl(dtMaxC, dblpck, idbl)
      call Ecpdbl(wCumA, dblpck, idbl)
      call Ecpdbl(cCumA, dblpck, idbl)
      call Ecpdbl(ECNVRG, dblpck, idbl)
      call Ecpdbl(ACNVRG, dblpck, idbl)
      call Ecpdbl(RCNVRG, dblpck, idbl)
      call Ecpdbl(t, dblpck, idbl)
      call Ecpdbl(TolTh, dblpck, idbl)
      call Ecpdbl(TolH, dblpck, idbl)
      call Ecpdbl(rLen, dblpck, idbl)
      call Ecpdbl(GWL0L, dblpck, idbl)
      call Ecpdbl(tMax, dblpck, idbl)
      call Ecpdbl(Aqh, dblpck, idbl)
      call Ecpdbl(Bqh, dblpck, idbl)
      call Ecpdbl(hCritS, dblpck, idbl)
      call Ecpdbl(tAtm, dblpck, idbl)
      call Ecpdbl(hCritA, dblpck, idbl)
      call Ecpdbl(cPrec, dblpck, idbl)
      call Ecpdbl(cht, dblpck, idbl)
      call Ecpdbl(crt, dblpck, idbl)
      call Ecpdbl(tOld, dblpck, idbl)
      call Ecpdbl(dt, dblpck, idbl)
      call Ecpdbl(dtMax, dblpck, idbl)
      call Ecpdbl(dMul, dblpck, idbl)
      call Ecpdbl(dMul2, dblpck, idbl)
      call Ecpdbl(dtMin, dblpck, idbl)
      call Ecpdbl(dtOpt, dblpck, idbl)
      call Ecpdbl(dtOld, dblpck, idbl)
      call Ecpdbl(dtInit, dblpck, idbl)
      call Ecpdbl(epsi, dblpck, idbl)
      call Ecpdbl(tPulse, dblpck, idbl)
      call Ecpdbl(PeCr, dblpck, idbl)
      call Ecpdbl(Peclet, dblpck, idbl)
      call Ecpdbl(Courant, dblpck, idbl)
      call Ecpdbl(wCumT, dblpck, idbl)
      call Ecpdbl(cCumT, dblpck, idbl)
      call Ecpdbl(wVolI, dblpck, idbl)
      call Ecpdbl(cVolI, dblpck, idbl)
      call Ecpdbl(vMeanR, dblpck, idbl)
      call Ecpdbl(hMeanR, dblpck, idbl)
      call Ecpdbl(hMeanT, dblpck, idbl)
      call Ecpdbl(hMeanG, dblpck, idbl)
      xxx = 999999.25
      call Ecpdbl(xxx, dblpck, idbl)
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
C      nint = 1

      call Ecpint(NTab, intpck, iint)
      call Ecpint(ItCum, intpck, iint)
      call Ecpint(Iter, intpck, iint)
      call Ecpint(TLevel, intpck, iint)
      call Ecpint(ALevel, intpck, iint)
      call Ecpint(PLevel, intpck, iint)
      call Ecpint(MaxItO, intpck, iint)
      call Ecpint(atmkk, intpck, iint)
      call Ecpint(nPLvl, intpck, iint)
      call Ecpint(TLkk, intpck, iint)
      call Ecpint(nbalnc, intpck, iint)
      call Ecpint(KAT, intpck, iint)
      call Ecpint(MaxIt, intpck, iint)
      call Ecpint(MBand, intpck, iint)
      call Ecpint(NLay, intpck, iint)
      call Ecpint(MaxAL, intpck, iint)
      call Ecpint(NSeep, intpck, iint)
      call Ecpint(NLevel, intpck, iint)
      call Ecpint(MPL, intpck, iint)
      call Ecpint(NDr, intpck, iint)
      call Ecpint(NObs, intpck, iint)
      nint = NumNP
      call cpint(Kode, intpck, nint, iint)
      call cpint(MatNum, intpck, nint, iint)
      nint = NumEl * 4
      call cpint(KX, intpck, nint, iint)
C      call cpint(, intpck, nint, iint)

CCCCCCCCCCCCCCCCCCCCCCCCC
C      nint = 1
      call Ecplog(SinkF, intpck, iint)
      call Ecplog(qGWLF, intpck, iint)
      call Ecplog(Explic, intpck, iint)
      call Ecplog(lMinStep, intpck, iint)
      call Ecplog(lWat, intpck, iint)
      call Ecplog(lChem, intpck, iint)
      call Ecplog(AtmInF, intpck, iint)
      call Ecplog(ShortF, intpck, iint)
      call Ecplog(SeepF, intpck, iint)
      call Ecplog(FluxF, intpck, iint)
      call Ecplog(lUpW, intpck, iint)
      call Ecplog(FreeD, intpck, iint)
      call Ecplog(DrainF, intpck, iint)
      call Ecplog(lArtD, intpck, iint)
      call Ecplog(lOrt, intpck, iint)
C      call Ecplog(, intpck, iint)

      nxxx = 999999
      call Ecpint(nxxx, intpck, iint)

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


C      subroutine cplog(from, to, nfrom, ito)
C      logical from(*)
C      integer to(*)
C      integer nfrom, ito
C      do 11 i=1,nfrom
C         if (from(i)) then 
C            to(ito + i) = 1
C         else
C            to(ito + i) = 0
C         end if
C 11   continue
C      ito = ito + nfrom 
C      return
C      end

      subroutine Ecpdbl(from, to, ito)
      double precision from, to(*)
      integer ito
      ito = ito + 1
      to(ito) = from
      return
      end


      subroutine Ecpint(from, to, ito)
      integer from, to(*)
      integer ito
      ito = ito + 1
      to(ito) = from
      return
      end

      subroutine Ecplog(from, to, ito)
      logical from
      integer to(*)
      integer ito
      ito = ito + 1
      if (from) then 
         to(ito) = 1
      else
         to(ito) = 0
      end if
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


C      subroutine relog(from, to, nfrom, ito)
C      logical from(*)
C      integer to(*)
C      integer nfrom, ito
C      do 11 i=1,nfrom
C         from(i) =  to(ito + i).eq.1
C 11   continue
C      ito = ito + nfrom 
C      return
C      end

      subroutine Eredbl(from, to, ito)
      double precision from, to(*)
      integer ito
      ito = ito + 1
      from = to(ito) 
      return
      end


      subroutine Ereint(from, to, ito)
      integer from, to(*)
      integer ito
      ito = ito + 1
      from = to(ito)
      return
      end


      subroutine Erelog(from, to, ito)
      logical from
      integer to(*)
      integer ito
      ito = ito + 1
      from =  to(ito).eq.1
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
     !        Explic,lUpW,FreeD,DrainF,lArtD,lOrt,lMinStep
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
C      ndbl = 1
      call Eredbl(tInit, dblpck, idbl)
      call Eredbl(CumQrT, dblpck, idbl)
      call Eredbl(CumQrR, dblpck, idbl)
      call Eredbl(CumQvR, dblpck, idbl)
      call Eredbl(rRoot, dblpck, idbl)
      call Eredbl(rTop, dblpck, idbl)
      call Eredbl(CumCh0, dblpck, idbl)
      call Eredbl(CumCh1, dblpck, idbl)
      call Eredbl(CumChR, dblpck, idbl)
      call Eredbl(dtMaxC, dblpck, idbl)
      call Eredbl(wCumA, dblpck, idbl)
      call Eredbl(cCumA, dblpck, idbl)
      call Eredbl(ECNVRG, dblpck, idbl)
      call Eredbl(ACNVRG, dblpck, idbl)
      call Eredbl(RCNVRG, dblpck, idbl)
      call Eredbl(t, dblpck, idbl)
      call Eredbl(TolTh, dblpck, idbl)
      call Eredbl(TolH, dblpck, idbl)
      call Eredbl(rLen, dblpck, idbl)
      call Eredbl(GWL0L, dblpck, idbl)
      call Eredbl(tMax, dblpck, idbl)
      call Eredbl(Aqh, dblpck, idbl)
      call Eredbl(Bqh, dblpck, idbl)
      call Eredbl(hCritS, dblpck, idbl)
      call Eredbl(tAtm, dblpck, idbl)
      call Eredbl(hCritA, dblpck, idbl)
      call Eredbl(cPrec, dblpck, idbl)
      call Eredbl(cht, dblpck, idbl)
      call Eredbl(crt, dblpck, idbl)
      call Eredbl(tOld, dblpck, idbl)
      call Eredbl(dt, dblpck, idbl)
      call Eredbl(dtMax, dblpck, idbl)
      call Eredbl(dMul, dblpck, idbl)
      call Eredbl(dMul2, dblpck, idbl)
      call Eredbl(dtMin, dblpck, idbl)
      call Eredbl(dtOpt, dblpck, idbl)
      call Eredbl(dtOld, dblpck, idbl)
      call Eredbl(dtInit, dblpck, idbl)
      call Eredbl(epsi, dblpck, idbl)
      call Eredbl(tPulse, dblpck, idbl)
      call Eredbl(PeCr, dblpck, idbl)
      call Eredbl(Peclet, dblpck, idbl)
      call Eredbl(Courant, dblpck, idbl)
      call Eredbl(wCumT, dblpck, idbl)
      call Eredbl(cCumT, dblpck, idbl)
      call Eredbl(wVolI, dblpck, idbl)
      call Eredbl(cVolI, dblpck, idbl)
      call Eredbl(vMeanR, dblpck, idbl)
      call Eredbl(hMeanR, dblpck, idbl)
      call Eredbl(hMeanT, dblpck, idbl)
      call Eredbl(hMeanG, dblpck, idbl)

      call Eredbl(xxx, dblpck, idbl)  
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
C      nint = 1
      call Ereint(NTab, intpck, iint)
      call Ereint(ItCum, intpck, iint)
      call Ereint(Iter, intpck, iint)
      call Ereint(TLevel, intpck, iint)
      call Ereint(ALevel, intpck, iint)
      call Ereint(PLevel, intpck, iint)
      call Ereint(MaxItO, intpck, iint)
      call Ereint(atmkk, intpck, iint)
      call Ereint(nPLvl, intpck, iint)
      call Ereint(TLkk, intpck, iint)
      call Ereint(nbalnc, intpck, iint)
      call Ereint(KAT, intpck, iint)
      call Ereint(MaxIt, intpck, iint)
      call Ereint(MBand, intpck, iint)
      call Ereint(NLay, intpck, iint)
      call Ereint(MaxAL, intpck, iint)
      call Ereint(NSeep, intpck, iint)
      call Ereint(NLevel, intpck, iint)
      call Ereint(MPL, intpck, iint)
      call Ereint(NDr, intpck, iint)
      call Ereint(NObs, intpck, iint)
      nint = NumNP
      call reint(Kode, intpck, nint, iint)
      call reint(MatNum, intpck, nint, iint)
      nint = NumEl * 4
      call reint(KX, intpck, nint, iint)
C      call reint(, intpck, nint, iint)

CCCCCCCCCCCCCCCCCCCCCCCCC
C      nint = 1
      call Erelog(SinkF, intpck, iint)
      call Erelog(qGWLF, intpck, iint)
      call Erelog(Explic, intpck, iint)
      call Erelog(lMinStep, intpck, iint)
      call Erelog(lWat, intpck, iint)
      call Erelog(lChem, intpck, iint)
      call Erelog(AtmInF, intpck, iint)
      call Erelog(ShortF, intpck, iint)
      call Erelog(SeepF, intpck, iint)
      call Erelog(FluxF, intpck, iint)
      call Erelog(lUpW, intpck, iint)
      call Erelog(FreeD, intpck, iint)
      call Erelog(DrainF, intpck, iint)
      call Erelog(lArtD, intpck, iint)
      call Erelog(lOrt, intpck, iint)
C      call Erelog(, intpck, iint)


      call Ereint(nxxx, intpck, iint)
      if (nxxx.ne.999999) then
         write(*,*) "nxxx.ne.999999"
         stop
      end if

      return
      end









