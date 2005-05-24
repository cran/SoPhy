*
*  Authors:   Simunek, J., T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, schlath@hsu-hh.de  (2004 -- 2005)
* 
*  Copyright (C) 2002 Simunek, J., T. Vogel and M. Th. van Genuchten

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


* Source file MATERIAL.FOR |||||||||||||||||||||||||||||||||||||||||||||

      double precision function FK(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m,Ks,Kr,Kk, h,Par(10)
      integer PPar

      BPar=.5d0
      PPar=2
      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      Ks=Par(7)
      Kk=Par(8)
      Qk=Par(9)
      m=1.d0-1.d0/n
      HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Qeek=dmin1((Qk-Qa)/(Qm-Qa),Qees)
      Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
      Hk=-1.d0/Alfa*(Qeek**(-1.d0/m)-1.d0)**(1.d0/n)
      if(dble(h).lt.Hk) then
        Qee=(1.d0+(-Alfa*HH)**n)**(-m)
        Qe =(Qm-Qa)/(Qs-Qa)*Qee
        Qek=(Qm-Qa)/(Qs-Qa)*Qeek
        FFQ =1.d0-(1.d0-Qee **(1.d0/m))**m
        FFQk=1.d0-(1.d0-Qeek**(1.d0/m))**m
        if(FFQ.le.0.d0) FFQ=m*Qee**(1.d0/m)
        Kr=(Qe/Qek)**Bpar*(FFQ/FFQk)**PPar*Kk/Ks
        FK=sngl(max(Ks*Kr,1.d-37))
        return
      elseif (dble(h).lt.Hs) then
        Kr=(1.d0-Kk/Ks)/(Hs-Hk)*(dble(h)-Hs)+1.d0
        FK=sngl(Ks*Kr)
      else
        FK=sngl(Ks)
      end if
      return
      end

************************************************************************

      double precision function FC(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m, h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.d0-1.d0/n
      HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
      if(dble(h).lt.Hs) then
        C1=(1.d0+(-Alfa*HH)**n)**(-m-1.d0)
        C2=(Qm-Qa)*m*n*(Alfa**n)*(-HH)**(n-1.d0)*C1
        FC=sngl(max(C2,1.d-37))
        return
      else
        FC=0.0
      end if
      return
      end

************************************************************************

      double precision function FQ(h,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m, h,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      m=1.d0-1.d0/n
      HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)
      HH=max(dble(h),HMin)
      Qees=dmin1((Qs-Qa)/(Qm-Qa),.999999999999999d0)
      Hs=-1.d0/Alfa*(Qees**(-1.d0/m)-1.d0)**(1.d0/n)
      if(dble(h).lt.Hs) then
        Qee=(1.d0+(-Alfa*HH)**n)**(-m)
        FQ=sngl(max(Qa+(Qm-Qa)*Qee,1.d-37))
        return
      else
        FQ=sngl(Qs)
      end if
      return
      end

************************************************************************

      double precision function FH(Qe,Par)

      implicit double precision(A-H,O-Z)
      double precision n,m, Qe,Par(9)

      Qr=Par(1)
      Qs=Par(2)
      Qa=Par(3)
      Qm=Par(4)
      Alfa=Par(5)
      n=Par(6)
      
      m=1.d0-1.d0/n
      HMin=-1.d300**(1.d0/n)/max(Alfa,1.d0)

      QeeM=(1.d0+(-Alfa*HMin)**n)**(-m)
      Qee=dmin1(dmax1(Qe*(Qs-Qa)/(Qm-Qa),QeeM),.999999999999999d0)
      FH=sngl(max(-1.d0/Alfa*(Qee**(-1.d0/m)-1.d0)**(1.d0/n),-1.d37))
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
