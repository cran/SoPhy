* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*
*                                                                     *
*     an analytic solution of water transport                         *
*                                                                     *
*     Martin Schlather and Bernd Huwe (2002, 2003)                    *
*                                                                     *
*                                                                     *
* ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||*

      subroutine ade(C0, D, v, t, nt, z, nz, conc)

      double precision C0, D, z(nz), t(nt), conc(nz, nt), s

      do 10 i=1,nt
         s = 2 * sqrt(D * t(i))
         do 20 j=1,nz
           conc(j,i) = .5 * C0 * (ERFC( (z(j) - v * t(i)) / s) +
     !        EXP(v * z(j) / D) * (ERFC( (z(j) + v * t(i)) / s)))
 20      continue
 10   continue
      return
      end
