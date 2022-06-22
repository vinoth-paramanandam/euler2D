module flux
  use constant
  use declaration
  use iso_fortran_env

  implicit none

  contains

    subroutine fluxcalc(rhof, uf, vf, pf, normalx, normaly, f)
      implicit none

      real(dp), intent(in) :: rhof, uf, vf, pf, normalx, normaly
      real(dp), dimension(4), intent(out) :: f
      real(dp) :: un

      un = uf*normalx + vf*normaly
      f(1) = rhof*un
      f(2) = rhof*un*uf + pf*normalx
      f(3) = rhof*un*vf + pf*normaly
      f(4) = (pf/g2 + half*rhof*(uf*uf + vf*vf) + pf)*un
    end subroutine fluxcalc

    subroutine lf(ql, qr, normx, normy, dlen, fluxvalue)
      implicit none

      real(dp), dimension(4), intent(in) :: ql, qr
      real(dp), intent(in) :: normx, normy, dlen
      real(dp), dimension(4), intent(out) :: fluxvalue
      real(dp), dimension(4) :: fl, fr
      real(dp) :: rhol, ul, vl, pl, rhor, ur, vr, pr
      real(dp) :: qnl, qnr, asoundl, asoundr, qro, aro, sl, sr

      rhol = ql(1)
      ul = ql(2)/rhol
      vl = ql(3)/rhol
      pl = g2*(ql(4) - half*rhol*(ul*ul + vl*vl))

      rhor = qr(1)
      ur = qr(2)/rhor
      vr = qr(3)/rhor
      pr = g2*(qr(4) - half*rhor*(ur*ur + vr*vr))

      call fluxcalc(rhol, ul, vl, pl, normx, normy, fl)
      call fluxcalc(rhor, ur, vr, pr, normx, normy, fr)

      qnl = ul*normx + vl*normy
      qnr = ur*normx + vr*normy

      asoundl = dsqrt(g1*pl/rhol)
      asoundr = dsqrt(g1*pr/rhor)

      qro = (dsqrt(rhol)*ul + dsqrt(rhor)*ur)/ &
           (dsqrt(rhol) + dsqrt(rhor)) !roeavg(rhol, rhor, qnl, qnr)
      aro = (dsqrt(rhol)*asoundl + dsqrt(rhor)*asoundr)/ &
           (dsqrt(rhol) + dsqrt(rhor))!roeavg(rhol, rhor, asoundl, asoundr)

      sr = dmax1(qro+aro, qnr+asoundr, zero)

      fluxvalue = half*(fl + fr - dabs(sr)*(qr - ql))*dlen
    end subroutine lf

    subroutine ausmpup(ql, qr, normx, normy, dlen, fluxvalue)
      implicit none

      real(dp), dimension(4), intent(in) :: ql, qr
      real(dp), intent(in) :: normx, normy, dlen
      real(dp), dimension(4), intent(out) :: fluxvalue
      real(dp), dimension(4) :: fl, fr
      real(dp) :: rhol, ul, vl, pl, tl, hl, asoundl
      real(dp) :: rhor, ur, vr, pr, tr, hr, asoundr
      real(dp) :: qnl, qnr, astarsqrl, astarsqrr
      real(dp) :: ahatl, ahatr, ahalf, eml, emr
      real(dp) :: emhatsqr, eminfsqr, emzerosqr, emzero
      real(dp) :: fa, rohalf, em1p, em1m, em2m, em2p
      real(dp) :: alpha_a, em4p, p5p, em4m, p5m, em2pp, em2mp
      real(dp) :: emhalf, emdothalf, phalf, flowm
      real(dp) :: akp = quarter, aku = threefourths, sigma = one
      real(dp) :: beta = one/eight

      rhol = ql(1)
      ul = ql(2)/rhol
      vl = ql(3)/rhol
      pl = g2*(ql(4) - half*rhol*(ul*ul + vl*vl))

      rhor = qr(1)
      ur = qr(2)/rhor
      vr = qr(3)/rhor
      pr = g2*(qr(4) - half*rhor*(ur*ur + vr*vr))

      tl = pl/(r*rhol)
      tr = pr/(r*rhor)

      hl = cpval*tl + half*(ul*ul + vl*vl)
      hr = cpval*tr + half*(ur*ur + vr*vr)

      qnl = ul*normx + vl*normy
      qnr = ur*normx + vr*normy

      astarsqrl = two*g2*hl/g3
      astarsqrr = two*g2*hr/g3

      ahatl = astarsqrl/dmax1(dsqrt(astarsqrl), qnl)
      ahatr = astarsqrr/dmax1(dsqrt(astarsqrr), -qnr)

      ahalf = dmin1(ahatl, ahatr)

      eml = qnl/ahalf
      emr = qnr/ahalf

      emhatsqr = (eml*eml + emr*emr)*half
      eminfsqr = nine

      emzerosqr = dmin1(one, dmax1(emhatsqr, eminfsqr))
      emzero = dsqrt(emzerosqr)

      fa = emzero*(two - emzero)
      rohalf = (rhol + rhor)*half

      em1p = half*(eml + dabs(eml))
      em1m = half*(emr - dabs(emr))

      em2p = quarter*(eml + one)**2
      em2m = -quarter*(emr - one)**2

      alpha_a = (three/sixteen)*(five*fa*fa - four)

      if (dabs(eml) >= one) then
         em4p = em1p
         p5p = em1p/eml
      else
         em2mp = -quarter*(eml - one)**2
         em4p = em2p*(one - sixteen*beta*em2mp)
         p5p = em2p*((two - eml) - sixteen*alpha_a*eml*em2mp)
      end if

      if (dabs(emr) >= one) then
         em4m = em1m
         p5m = em1m/emr
      else
         em2pp = quarter*(emr + one)**2
         em4m = em2m*(one + sixteen*beta*em2pp)
         p5m = em2m*((-two - emr) + sixteen*alpha_a*emr*em2pp)
      end if

      emhalf = em4p + em4m -(akp/fa)*(pr - pl)/ &
           (rohalf*ahalf**2)*dmax1((one - sigma*emhatsqr), zero)

      if (emhalf > zero) then
         emdothalf = ahalf*emhalf*rhol
      else
         emdothalf = ahalf*emhalf*rhor
      end if

      phalf = p5p*pl + p5m*pr - aku*p5p*(rhol + rhor)* &
           (fa*ahalf)*(qnr - qnl)
      flowm = emdothalf*dlen

      if (emdothalf > zero) then
         fluxvalue(1) = flowm
         fluxvalue(2) = flowm*ul + phalf*normx*dlen
         fluxvalue(3) = flowm*vl + phalf*normy*dlen
         fluxvalue(4) = flowm*hl
      else
         fluxvalue(1) = flowm
         fluxvalue(2) = flowm*ur + phalf*normx*dlen
         fluxvalue(3) = flowm*vr + phalf*normy*dlen
         fluxvalue(4) = flowm*hr
      end if
    end subroutine ausmpup

    subroutine isecondMUSCL
      implicit none

      integer(i8) :: nb
      integer(i16) :: i, j

      real(dp) :: dxlen, dylen, dlen, normalx, normaly
      real(dp), dimension(4) :: fluxval, vecl, vecr
      real(dp), dimension(4) :: ar, al, br, bl, delR, delL

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 0, nx(nb)
               dxlen = -(x(i+1, j+1, nb) - x(i+1, j, nb))
               dylen = (y(i+1, j+1, nb) - y(i+1, j, nb))
               dlen = dsqrt(dxlen*dxlen + dylen*dylen)
               normalx = dylen/dlen
               normaly = dxlen/dlen

               ar = q(:, i+2, j, nb) - q(:, i+1, j, nb)
               br = q(:, i+1, j, nb) - q(:, i, j, nb)

               al = q(:, i+1, j, nb) - q(:, i, j, nb)
               bl = q(:, i, j, nb) - q(:, i-1, j, nb)

               delR = (ar*(br*br + e_tvd) + br*(ar*ar + e_tvd))/&
                    (ar*ar + br*br + two*e_tvd)
               delL =  (al*(bl*bl + e_tvd) + bl*(al*al + e_tvd))/&
                    (al*al + bl*bl + two*e_tvd)

               vecl = q(:, i, j, nb) + half*delL
               vecr = q(:, i+1, j, nb) - half*delR
               ! vecl = q(:, i, j, nb)
               ! vecr = q(:, i+1, j, nb)

               !call lf(vecl, vecr, normalx, normaly, dlen, fluxval)
               call ausmpup(vecl, vecr, normalx, normaly, dlen, fluxval)

               l_res(:, i, j, nb) = l_res(:, i, j, nb) + fluxval
               l_res(:, i+1, j, nb) = l_res(:, i+1, j, nb) - fluxval
            end do
         end do
      end do
    end subroutine isecondMUSCL

    subroutine jsecondMUSCL
      implicit none

      integer(i8) :: nb
      integer(i16) :: i, j

      real(dp) :: dxlen, dylen, dlen, normalx, normaly
      real(dp), dimension(4) :: fluxval, vecl, vecr
      real(dp), dimension(4) :: ar, al, br, bl, delR, delL

      do nb = 1, nblocks
         do j = 0, ny(nb)
            do i = 1, nx(nb)
               dxlen = -(x(i, j+1, nb) - x(i+1, j+1, nb))
               dylen = (y(i, j+1, nb) - y(i+1, j+1, nb))
               dlen = dsqrt(dxlen*dxlen + dylen*dylen)
               normalx = dylen/dlen
               normaly = dxlen/dlen

               ar = q(:, i, j+2, nb) - q(:, i, j+1, nb)
               br = q(:, i, j+1, nb) - q(:, i, j, nb)

               al = q(:, i, j+1, nb) - q(:, i, j, nb)
               bl = q(:, i, j, nb) - q(:, i, j-1, nb)

               delR = (ar*(br*br + e_tvd) + br*(ar*ar + e_tvd))/&
                    (ar*ar + br*br + two*e_tvd)
               delL =  (al*(bl*bl + e_tvd) + bl*(al*al + e_tvd))/&
                    (al*al + bl*bl + two*e_tvd)

               vecl = q(:, i, j, nb) + half*delL
               vecr = q(:, i, j+1, nb) - half*delR
               ! vecl = q(:, i, j, nb)
               ! vecr = q(:, i, j+1, nb)

               !call lf(vecl, vecr, normalx, normaly, dlen, fluxval)
               call ausmpup(vecl, vecr, normalx, normaly, dlen, fluxval)

               l_res(:, i, j, nb) = l_res(:, i, j, nb) + fluxval
               l_res(:, i, j+1, nb) = l_res(:, i, j+1, nb) - fluxval
            end do
         end do
      end do
    end subroutine jsecondMUSCL

  end module flux
