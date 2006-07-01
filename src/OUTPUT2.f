*
*  Authors:   Simunek, J., T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, martin.schlather@math.uni-goettingen.de  (2004 -- 2006)
* 
*  Copyright (C) 2002 Simunek, J., T. Vogel and M. Th. van Genuchten, 
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



* Source file OUTPUT2.FOR ||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TLInf(NumNP,NumBP,Kode,Q,hNew,CumQ,Width,SWidth,KXB,t,
     !                 dt,TLevel,ShortF,TPrint,Iter,ItCum,rTop,rRoot,
     !                 vMeanR,hMeanT,hMeanR,hMeanG,AtmInF,SinkF,CumQrT,
     !                 CumQrR,CumQvR,NumKD,hMean,vMean,lWat,lChem,rLen,
     !                 Peclet,Courant,wCumT,wCumA, 
     !                 TlSolO, TlSolR, TlSolC, TLkk)

      IMPLICIT REAL*8 (A-H,O-Z)
      integer TlSolR,TlSolC, Tlkk
      double precision TlSolO(TlSolR,TlSolC)

      integer TLevel
      logical ShortF,SinkF,AtmInF,lWat,lChem
      dimension Q(NumNP),Kode(NumNP),KXB(NumBP),SWidth(NumKD),
     !          hNew(NumNP),Width(NumBP),CumQ(NumKD),hMean(NumKD),
     !          vMean(NumKD)

      if(lWat.or.TLevel.eq.1) then
        do 11 i=1,NumKD
          vMean(i)=0.
          hMean(i)=0.
11      continue
        do 12 i=1,NumBP
          n=KXB(i)
          j=iabs(Kode(n))
          if(j.eq.0) goto 12
          hMean(j)=hMean(j)+hNew(n)*Width(i)/SWidth(j)
         if(j.eq.4) vMean(j)=vMean(j)-Q(n)/SWidth(j)
12      continue
        hMeanT=hMean(4)
        hMeanG=hMean(3)
        do 13 i=1,NumNP
          j=iabs(Kode(i))
          wCumA=wCumA+dabs(Q(i))*dt
          if(j.ne.0.and.j.ne.4) then
            vMean(j)=vMean(j)-Q(i)
          end if
13      continue
        if(.not.lWat.and.TLevel.eq.1) then
           TlSolO(TLkk, 1) = t
           TlSolO(TLkk, 2) = rTop
           TlSolO(TLkk, 3) = rRoot
           TlSolO(TLkk, 4) = vMean(4)
           TlSolO(TLkk, 5) = vMeanR
           TlSolO(TLkk, 6) = vMean(3)
           TlSolO(TLkk, 7) = vMean(1)
           TlSolO(TLkk, 8) = vMean(2)
           do 700 ii=5,NumKD
              TlSolO(TLkk, 4+ii) = vMean(ii)
 700       continue
           jj = 4 + NumKD 
           TlSolO(TLkk, jj+1) = hMean(4)
           TlSolO(TLkk, jj+2) = hMeanR
           TlSolO(TLkk, jj+3) = hMean(3)
           TlSolO(TLkk, jj+4) = hMean(1)
           TlSolO(TLkk, jj+5) = hMean(2)
           jj = 5 + NumKD 
           do 701 ii=5,NumKD
              TlSolO(TLkk, jj + ii) = hMean(ii)
 701       continue
           jj = 5 + 2 * NumKD           
        end if
      end if
      if(lWat) then
        wCumA=wCumA+dabs(vMeanR*dt*rLen)
        wCumT=CumQvR
        do 14 j=1,NumKD
          if(j.eq.4) then
            CumQ(j)=CumQ(j)+vMean(j)*dt*SWidth(4)
          else
            CumQ(j)=CumQ(j)+vMean(j)*dt
          end if
          wCumT=wCumT+CumQ(j)
14      continue
        CumQrT=CumQrT+rTop  *dt*Swidth(4)
        CumQrR=CumQrR+rRoot *dt*rLen
        CumQvR=CumQvR+vMeanR*dt*rLen
      end if
      if(.not.ShortF.or.dabs(TPrint-t).lt.0.001*dt) then 
         if(lWat) then
           TlSolO(TLkk, 1) = t
           TlSolO(TLkk, 2) = rTop
           TlSolO(TLkk, 3) = rRoot
           TlSolO(TLkk, 4) = vMean(4)
           TlSolO(TLkk, 5) = vMeanR
           TlSolO(TLkk, 6) = vMean(3)
           TlSolO(TLkk, 7) = vMean(1)
           TlSolO(TLkk, 8) = vMean(2)
           jj = 4
           do 702 ii=5,NumKD
              TlSolO(TLkk, jj+ii) = vMean(ii)
 702       continue
           jj = 4 + NumKD 
           TlSolO(TLkk, jj+1) = hMean(4)
           TlSolO(TLkk, jj+2) = hMeanR
           TlSolO(TLkk, jj+3) = hMean(3)
           TlSolO(TLkk, jj+4) = hMean(1)
           TlSolO(TLkk, jj+5) = hMean(2)
           jj = 5 + NumKD 
           do 703 ii=5,NumKD
              TlSolO(TLkk, jj + ii) = hMean(ii)
 703       continue
           jj = 5 + 2 * NumKD  
           TlSolO(TLkk, jj+1) = CumQrT
           TlSolO(TLkk, jj+2) = CumQrR
           TlSolO(TLkk, jj+3) = CumQ(4)
           TlSolO(TLkk, jj+4) = CumQvR
           TlSolO(TLkk, jj+5) = CumQ(3)
           TlSolO(TLkk, jj+6) = CumQ(1)
           TlSolO(TLkk, jj+7) = CumQ(2)
           jj = 8 + 2 * NumKD
           do 704 ii=5,NumKD
              TlSolO(TLkk, jj + ii) = hMean(ii)
 704       continue   
           jj = 8 + 3 * NumKD
        end if
        jj = 3 * (NumKD + 1) + 5
        if(lChem) then 
          if(lWat) then
             TlSolO(TLkk, jj+1) = dt
             TlSolO(TLkk, jj+2) = Iter
             TlSolO(TLkk, jj+3) = ItCum
             TlSolO(TLkk, jj+4) = Peclet
             TlSolO(TLkk, jj+5) = Courant
          else
            TlSolO(TLkk, jj+1) = dt
            TlSolO(TLkk, jj+4) = Peclet
            TlSolO(TLkk, jj+5) = Courant
          end if
        else
           TlSolO(TLkk, jj+1) = dt
           TlSolO(TLkk, jj+2) = Iter
           TlSolO(TLkk, jj+3) = ItCum
        end if
      end if
      return
      end

************************************************************************

      subroutine ALInf(t,CumQ,hMeanT,hMeanR,hMeanG,ALevel,CumQrT,CumQrR,
     !                 CumQvR,NumKD, atmOut, MaxAL)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision atmOut(MaxAL, 9)

      integer ALevel
      dimension CumQ(NumKD)
      
      atmOut(ALevel, 1) = t
      atmOut(ALevel, 2) = CumQrT
      atmOut(ALevel, 3) = CumQrR
      atmOut(ALevel, 4) = CumQ(4)
      atmOut(ALevel, 5) = CumQvR
      atmOut(ALevel, 6) = CumQ(3)
      atmOut(ALevel, 7) = hMeanT
      atmOut(ALevel, 8) = hMeanR
      atmOut(ALevel, 9) = hMeanG 

      return
      end

************************************************************************

      subroutine SubReg(NumEl,NumNP,NMat,hNew,ThO,ThN,x,y,MatNum,
     !                  LayNum,KX,KAT,t,dt,NLay,PLevel,lWat,lChem,Conc,
     !                  ChPar,wCumA,wCumT,cCumA,cCumT,wVolI,cVolI,WatIn,
     !                  SolIn, balanc, balaR, balaC, nbalnc)
      IMPLICIT REAL*8 (A-H,O-Z)
      integer balaR, balaC, nbalnc
      double precision balanc(balaR, balaC)
 
      logical lWat,lChem, ldebug
      integer PLevel,e
      dimension hNew(NumNP),x(NumNP),y(NumNP),MatNum(NumNP),Conc(NumNP),
     !          KX(NumEl,4),ChPar(10,NMat),LayNum(NumEl),ThO(NumNP),
     !          ThN(NumNP),WatIn(NumEl),SolIn(NumEl),Area(10),hMean(10),
     !          SubVol(10),SubCha(10),cMean(10),ConSub(10)

      ldebug = .false.
C      ldebug = .true.
      if (ldebug) write(*,*) "subreg: start"
      xMul=1.
      ATot=0.
C   ### cancelling of the if statements if of no harm, but
C   ### ensures initialised values 
C      if(lWat.or.PLevel.le.1) then
        Volume=0.
        Change=0.
        hTot=0.
        DeltW=0.
C      end if
C      if(lChem) then
        cTot=0.
        ConVol=0.
        DeltC=0.
C      end if
      do 11 i=1,NLay
        Area(i)=0.
        if(lWat.or.PLevel.le.1) then
          SubVol(i)=0.
          SubCha(i)=0.
          hMean(i)=0.
        end if
        if(lChem) then
          ConSub(i)=0.
          cMean(i)=0.
        end if
11    continue
      if (ldebug) write(*,*) "subreg: 11"
      do 13 e=1,NumEl
        Lay=LayNum(e)
        wEl=0.
        cEl=0.
        NUS=4
        if(KX(e,3).eq.KX(e,4)) NUS=3
        do 12 k=1,NUS-2
           if (ldebug) write(*,*) "subreg: 12"
          i=KX(e,1)
          j=KX(e,k+1)
          l=KX(e,k+2)
          Mi=MatNum(i)
          Mj=MatNum(j)
          Mk=MatNum(l)
          Cj=x(i)-x(l)
          Ck=x(j)-x(i)
          Bj=y(l)-y(i)
          Bk=y(i)-y(j)
C          write(*,*) e,i,j,l,Mi,Mj,Mk,Lay, Cj,Ck,Bj,Bk
C 777      format(8i10,  4f15.8)
          if(KAT.eq.1) xMul=2.*3.1416*(x(i)+x(j)+x(l))/3.D0
          AE=xMul*(Ck*Bj-Cj*Bk)/2.D0
          Area(Lay)=Area(Lay)+AE
          if(lWat.or.PLevel.le.1) then
            hE=(hNew(i)+hNew(j)+hNew(l))/3.D0
            VNewE=AE*(thN(i)+thN(j)+thN(l))/3.D0
            VOldE=AE*(thO(i)+thO(j)+thO(l))/3.D0
            Volume=Volume+VNewE
            wEl=wEl+VNewE
            Change=Change+(VNewE-VOldE)/dt
            SubVol(Lay)=SubVol(Lay)+VNewE
            SubCha(Lay)=SubCha(Lay)+(VNewE-VOldE)/dt
            hTot=hTot+hE*AE
            hMean(Lay)=hMean(Lay)+hE*AE
          end if
          if (ldebug) write(*,*) "subreg: 12,lchem"
          if(lChem) then
            cE=(Conc(i)+Conc(j)+Conc(l))/3.D0
            cNewE=AE*((thN(i)+ChPar(1,Mi)*ChPar(5,Mi))*Conc(i)+
     !                (thN(j)+ChPar(1,Mj)*ChPar(5,Mj))*Conc(j)+
     !                (thN(l)+ChPar(1,Mk)*ChPar(5,Mk))*Conc(l))/3.D0
            ConVol=ConVol+cNewE
            cEl=cEl+cNewE
            ConSub(Lay)=ConSub(Lay)+cNewE
            cTot=cTot+cE*AE
            cMean(Lay)=cMean(Lay)+cE*AE
          end if
          if(k.eq.NUS-2) then
            if(PLevel.eq.0) then
              if(lWat) WatIn(e)=wEl
              if(lChem) SolIn(e)=cEl
            else
              if(lWat) DeltW=DeltW+dabs(WatIn(e)-wEl)
              if(lChem) DeltC=DeltC+dabs(SolIn(e)-cEl)
            end if
          end if
12      continue
13    continue
      if (ldebug) write(*,*) "subreg: 13"
      do 14 Lay=1,NLay
        if(lWat.or.PLevel.le.1) hMean(Lay)=hMean(Lay)/Area(Lay)
        if(lChem) cMean(Lay)=cMean(Lay)/Area(Lay)
        ATot=ATot+Area(Lay)
14    continue
      if(lWat.or.PLevel.le.1) hTot=hTot/ATot
      if(lChem) cTot=cTot/ATot      
      balanc(nbalnc, 1) = t
      balanc(nbalnc, 6) = ATot
      jj = 6
      do 705 ii=1,NLay
         balanc(nbalnc, jj + ii) = Area(ii)
 705  continue
      if(lWat.or.PLevel.le.1) then
         jj = 6 + NLay      
         balanc(nbalnc, jj+1) = Volume 
         jj = 7 + NLay      
         do 706 ii=1,NLay
            balanc(nbalnc, jj + ii) = SubVol(ii)
 706     continue
         jj = 7 + 2 * NLay      
         balanc(nbalnc, jj+1) = Change 
         jj = 8 + 2 * NLay
         do 707 ii=1,NLay
            balanc(nbalnc, jj + ii) = SubCha(ii)
 707     continue
         jj = 8 + 3 * NLay      
         balanc(nbalnc, jj+1) = hTot
         jj = 9 + 3 * NLay      
         do 708 ii=1,NLay
            balanc(nbalnc, jj + ii) = hMean(ii)
 708     continue         
      end if
CCCC falls lchem, dann notwendig TPrint point, see SWMS2d
      if (ldebug) write(*,*) "subreg: chem"
      if(lChem) then
         jj = 9 + 4 * NLay      
         balanc(nbalnc, jj+1) = ConVol
         jj = 10 + 4 * NLay   
         do 709 ii=1,NLay
            balanc(nbalnc, jj + ii) = ConSub(ii)
 709     continue         
         jj = 10 + 5 * NLay 
         balanc(nbalnc, jj+1) = cTot
         jj = 11 + 5 * NLay 
         do 710 ii=1,NLay
            balanc(nbalnc, jj + ii) = cMean(ii)
 710     continue   
         jj = 11 + 6 * NLay 
      end if
      
*     Mass balance calculation
      if(PLevel.eq.0) then
        wVolI=Volume
        if(lChem) cVolI=ConVol
      else
        jj = 1    
        if(lWat) then
          wBalT=Volume-wVolI+wCumT 
          balanc(nbalnc, jj+1) = wBalT
          ww=dmax1(DeltW,wCumA)
          if(ww.ge.1.d-25) then
            wBalR=dabs(wBalT)/ww*100.D0
            balanc(nbalnc, jj+2) = wBalR
          end if
        end if
        if(lChem) then
          cBalT=ConVol-cVolI+cCumT
          balanc(nbalnc, jj+3) = cBalT
          cc=dmax1(DeltC,cCumA)
          if(cc.ge.1.d-25) then
            cBalR=dabs(cBalT)/cc*100.D0
            balanc(nbalnc, jj+4) = cBalR
          end if
        end if
      end if

      return
      end

************************************************************************

      subroutine BouOut(NumNP,NumBP,t,hNew,theta,Q,Width,KXB,Kode,x,y,
     !                  Conc, boundr, balaR, nbalnc)

      IMPLICIT REAL*8 (A-H,O-Z)
      integer balaR, nbalnc
      double precision boundr(NumBP, 11, balaR)
      
      dimension hNew(NumNP),Q(NumNP),Width(NumBP),theta(NumNP),
     !          KXB(NumBP),Kode(NumNP),x(NumNP),y(NumNP),Conc(NumNP)

      ii=0
      do 12 i=1,NumNP
        if(Kode(i).ne.0) then
          do 11 j=1,NumBP
            n=KXB(j)
            if(n.eq.i) then
              ii=ii+1
              v=-Q(i)/Width(j)
              boundr(ii ,1, nbalnc) = ii
              boundr(ii ,2 ,nbalnc) = i
              boundr(ii ,3 ,nbalnc) = x(i)
              boundr(ii ,4 ,nbalnc) = y(i)
              boundr(ii ,5 ,nbalnc) = Kode(i)
              boundr(ii ,6 ,nbalnc) = Q(i)
              boundr(ii ,7, nbalnc) = v
              boundr(ii ,8, nbalnc) = hNew(i)
              boundr(ii ,9, nbalnc) = theta(i)
              boundr(ii ,10, nbalnc) = Conc(i)
              boundr(ii, 11, nbalnc) = t
              goto 12
            end if
11        continue
          ii=ii+1
          boundr(ii ,1, nbalnc) = ii
          boundr(ii ,2, nbalnc) = i
          boundr(ii ,3, nbalnc) = x(i)
          boundr(ii ,4, nbalnc) = y(i)
          boundr(ii ,5, nbalnc) = Kode(i)
          boundr(ii ,6, nbalnc) = Q(i)
C          boundr(ii ,7, nbalnc) = v
          boundr(ii ,8, nbalnc) = hNew(i)
          boundr(ii ,9, nbalnc) = theta(i)
          boundr(ii ,10, nbalnc) = Conc(i)
          boundr(ii, 11, nbalnc) = t
        end if
12    continue

      return
      end

************************************************************************

      subroutine SolInf(NumNP,Kode,Qc,t,dt,TLevel,ShortF,TPrint,NumKD,
     !                  SMean,ChemS,CumCh0,CumCh1,CumChR,cCumA,cCumT,
     !                  lWat, TlSolO, TlSolR, TlSolC, TLkk)
   
      IMPLICIT REAL*8 (A-H,O-Z)
      integer TlSolR, TlSolC, TLkk
      double precision TlSolO(TlSolR,TlSolC)

      integer TLevel
      logical ShortF,lWat
      dimension Qc(NumNP),Kode(NumNP),ChemS(NumKD),SMean(NumKD)
      
      do 11 i=1,NumKD
        SMean(i)=0.
11    continue
      do 12 i=1,NumNP
        j=iabs(Kode(i))
        if(j.ne.0) then
          SMean(j)=Smean(j)-Qc(i)
        end if
12    continue
      cCumA=dabs(CumCh0)+dabs(CumCh1)+dabs(CumChR)
      cCumT=CumCh0+CumCh1+CumChR
      do 13 j=1,NumKD
        ChemS(j)=ChemS(j)+SMean(j)*dt
        cCumT=cCumT+ChemS(j)
        cCumA=cCumA+dabs(ChemS(j))
13    continue
      jj = 10 + 3 * (NumKD + 1)
      if(.not.ShortF.or.dabs(TPrint-t).lt.0.001*dt) then
         TlSolO(TLkk, jj+1) = CumCh0
         TlSolO(TLkk, jj+2) = CumCh1
         TlSolO(TLkk, jj+3) = CumChR 
         jj = jj + 3
         do 711 ii=1,NumKD
            TlSolO(TLkk, jj + ii) = ChemS(ii)
 711     continue
         jj = jj + NumKD
         do 712 ii=1,NumKD
            TlSolO(TLkk, jj + ii) = SMean(ii)
 712     continue
      end if

      return
      end

***********************************************************************

      subroutine ObsNod(t,NumNP,NObs,NObsD,Node,hNew,ThNew,Conc, 
     !                 TlSolO, TlSolR, TlSolC, NumKD, TLkk)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      integer TlSolR, TlSolC, TLkk
      double precision TlSolO(TlSolR,TlSolC)

      dimension Node(NObsD),hNew(NumNP),ThNew(NumNP),Conc(NumNP)

      jj = 10 + 3 * (NumKD + 1) + 3 + 2 * NumKD
      do 713 ii=1,NObs
         TlSolO(TLkk, jj + ii) = hNew(Node(ii))
 713  continue
      jj = jj + NObs
      do 714 ii=1,NObs
         TlSolO(TLkk, jj + ii) = ThNew(Node(ii))
 714  continue
      jj = jj + NObs
      do 715 ii=1,NObs
         TlSolO(TLkk, jj + ii) = Conc(Node(ii))
 715  continue

      return
      end
      
************************************************************************

      subroutine hOut(hNew,x,y,NumNP,t, hQThF, MPL, kk)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision hQThF(6, NumNP, MPL + 2)

      dimension hNew(NumNP),x(NumNP),y(NumNP)

      do 716 ii=1,NumNP
         hQThF(1, ii, kk) = hNew(ii)
 716  continue

      return
      end

***********************************************************************

      subroutine QOut(Q,x,y,NumNP,t, hQThF, MPL, kk)
   
      IMPLICIT REAL*8 (A-H,O-Z)
      double precision hQThF(6, NumNP, MPL + 2)

      dimension Q(NumNP),x(NumNP),y(NumNP)

      do 717 ii=1,NumNP
         hQThF(2, ii, kk) = Q(ii)
 717  continue

      return
      end

************************************************************************

      subroutine thOut(theta,x,y,NumNP,t,hQThF, MPL, kk)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision hQThF(6, NumNP, MPL + 2)

      dimension theta(NumNP),x(NumNP),y(NumNP)

      do 718 ii=1,NumNP
         hQThF(3, ii, kk) = theta(ii)
 718  continue

      return
      end

************************************************************************

      subroutine FlxOut(Vx,Vz,x,y,NumNP,t, hQThF, MPL, kk)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision hQThF(6, NumNP,  MPL + 2)

      dimension x(NumNP),y(NumNP),Vx(NumNP),Vz(NumNP)

      do 719 ii=1,NumNP
         hQThF(4, ii, kk) = Vx(ii)
         hQThF(5, ii, kk) = Vz(ii)
 719  continue

      return
      end

************************************************************************

      subroutine cOut(NumNP,Conc,x,y,t, hQThF, MPL, kk)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision hQThF(6, NumNP, MPL + 2)

      dimension Conc(NumNP),x(NumNP),y(NumNP)

      do 720 ii=1,NumNP
         hQThF(6, ii, kk) = Conc(ii)
 720  continue

      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||





