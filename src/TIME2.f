*
*  Authors:   J. Simunek, T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, Martin.Schlather@uni-bayreuth.de 
*  for Computers & Geosciences:
*     The use of the language interface of R: two examples for        
*      modelling water flux and solute transport                       
*
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

* NOTE: The original file was published as a public domain contribution.
*       The modified file is published under the GNU licence in
*       accordance with Rien van Genuchten, email correspondance August, 2002


* Source file TIME2.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine TmCont(dt,dtMaxW,dtOpt,dMul,dMul2,dtMin,Iter,tPrint,
     !                  tAtm,t,tMax,dtMaxC,lMinStep,dtInit)
      IMPLICIT REAL*8 (A-H,O-Z)
      logical lMinStep

      if(lMinStep) then
        dtMax=dmin1(dtMaxW,dtMaxC,dtInit)
        lMinStep=.false.
      else
        dtMax=dmin1(dtMaxW,dtMaxC)
      end if
      tFix=dmin1(tPrint,tAtm,tMax)
      if(Iter.le.3.and.(tFix-t).ge.dMul*dtOpt) 
     !  dtOpt=dmin1(dtMax,dMul*dtOpt)
      if(Iter.ge.7)
     !  dtOpt=dmax1(dtMin,dMul2*dtOpt)
      dt=dmin1(dtOpt,tFix-t)
      dt=dmin1((tFix-t)/anint((tFix-t)/dt),dtMax)
      if(tFix-t.ne.dt.and.dt.gt.(tFix-t)/2.D0) dt=(tFix-t)/2.D0
      return
      end

************************************************************************

      subroutine SetAtm(tAtm,rTop,rRoot,hCritA,Width,KXB,NumBP,Kode,
     !                  hNew,Q,NumNP,GWL0L,qGWLF,FreeD,cPrec,cht,crt,
     !                  lMinStep, 
     !                  atmX, MaxAL, atmkk)
      
      IMPLICIT REAL*8 (A-H,O-Z)
      double precision atmX(MaxAL, 10)
      integer atmkk

      logical qGWLF,FreeD,lMinStep
      dimension Width(NumBP),KXB(NumBP),Kode(NumNP),hNew(NumNP),Q(NumNP)

      atmkk = atmkk + 1

C      write(*,777) atmkk
 777  format("atmkk ",i10)
      
      if (atmkk.gt.MaxAL) then
         write(*,777) atmkk
         write(*,*) 'atmkk > MaxAL'
         stop
      end if 
      rTopOld=rTop
C      read (31,*) tAt`m,Prec,cPrec,rSoil,rRoot,hCritA,rGWL,GWL,crt,cht
      tAtm = atmX(atmkk, 1)
      Prec = atmX(atmkk, 2)
      cPrec = atmX(atmkk, 3)
      rSoil = atmX(atmkk, 4)
      rRoot = atmX(atmkk, 5)
      hCritA = atmX(atmkk, 6)
      rGWL = atmX(atmkk, 7)
      GWL = atmX(atmkk, 8)
      crt = atmX(atmkk, 9)
      cht = atmX(atmkk, 10)
C      if (atmkk.eq.1) write(*,739)      
C      write(*,740) tAtm,Prec,cPrec,rSoil,rRoot,hCritA,rGWL,GWL,crt,cht
C      if (atmkk.eq.1) write(51,739)      
C      write(51,740) tAtm,Prec,cPrec,rSoil,rRoot,hCritA,rGWL,GWL,crt,cht
      
 739  format('         tAtm         Prec        cPrec        rSoil',
     !       '        rRoot       hCritA         rGWL          GWL',
     !       '          crt          cht')
 740  format(10f13.3)
      Prec=abs(Prec)
      rSoil=abs(rSoil)
      rRoot=abs(rRoot)
      hCritA=-abs(hCritA)
      hGWL=GWL+GWL0L
      rTop=rSoil-Prec
      do 11 i=1,NumBP
        n=KXB(i)
        K=Kode(n)
        if(K.eq.4.or.K.eq.-4) then
          Kode(n)=-4
          Q(n)=-Width(i)*rTop
          goto 11
        end if
        if(K.eq. 3) then
          if(dabs(hNew(n)-hGWL).gt.1.D-8) lMinStep=.true.
          hNew(n)=hGWL
        end if
        if(K.eq.-3.and..not.qGWLF.and..not.FreeD) then
          if(Width(i).gt.0.D0) rGWLOld=-Q(n)/Width(i)
          if(dabs(rGWLOld-rGWL).gt.1.D-8) lMinStep=.true.
          Q(n)=-Width(i)*rGWL
        end if
11    continue
      if((Prec-rSoil).gt.0.D0) then
        cPrec=Prec/(Prec-rSoil)*cPrec
      else
        cPrec=0.
      end if
      if(dabs(rTop-rTopOld).gt.1.D-8.and.dabs(rTop).gt.0.D0) 
     !  lMinStep=.true.
      return
      end
      
************************************************************************

      double precision function Fqh(GWL,Aqh,Bqh)
      IMPLICIT REAL*8 (A-H,O-Z)
      Fqh=-Aqh*exp(Bqh*dabs(GWL))
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
