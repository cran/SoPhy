*
*  Authors:   Simunek, J., T. Vogel and M. Th. van Genuchten.
*
*  modified by 
*  Martin Schlather, schlath@hsu-hh.de (2004 -- 2005)
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


* Source file SINK2.FOR ||||||||||||||||||||||||||||||||||||||||||||||||

      subroutine SetSnk(NumNP,hNew,TPot,Sink,PRL,Beta,Length)

      IMPLICIT REAL*8 (A-H,O-Z)
      double precision Length
      dimension hNew(NumNP),Beta(NumNP), Sink(NumNP), PRL(NumNP, 7)

      do 11 i=1,NumNP
        if(Beta(i).gt.0.) then
          Alfa=FAlfa(TPot,hNew(i),PRL(i,1),PRL(i,7),PRL(i,2),PRL(i,3),
     !               PRL(i,4),PRL(i,5),PRL(i,6))
          Sink(i)=Alfa*Beta(i)*Length*TPot
        end if
11    continue
      return
      end 
   
************************************************************************
 
      double precision function FAlfa(TPot,h,P0,P1,P2H,P2L,P3,r2H,r2L)
      IMPLICIT REAL*8 (A-H,O-Z)

      if((TPot.ge.r2L).and.(TPot.le.r2H)) then
        P2=P2H+(r2H-TPot)/(r2H-r2L)*(P2L-P2H)
      else
       if(TPot.gt.r2H) then
          P2=P2H
       else 
CC        if(TPot.lt.r2L) 
         P2=P2L
       endif 
      end if
      FAlfa=0.0
      if((h.gt.P3).and.(h.lt.P2)) FAlfa=(h-P3)/(P2-P3)
      if((h.ge.P2).and.(h.le.P1)) FAlfa=1.0
      if((h.gt.P1).and.(h.lt.P0)) FAlfa=(h-P0)/(P1-P0)
      return
      end

* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
