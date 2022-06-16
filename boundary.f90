module boundary
  use constant
  use declaration
  !use grid

  implicit none

  contains

    subroutine nozzleboundary
      implicit none

      integer(i8) :: nb, i, j
      real(dp) :: dxlen, dylen, dlen, normalx, normaly
      real(dp) :: un, ut, pamb, uamb, vamb, tamb, qnorm, m, rref, cref
      real(dp) :: pdom, udom, vdom, rdom, pprev, rprev, uprev, vprev
      real(dp) :: rin, uin, vin, pin, pb, rb, vb, ub, ramb
      real(dp) :: temp, dpres, tinf

      nb = 1
      ! bottom boundary conditions
      do j = -2, 0
         do i = 1, nx(nb)
            dxlen = (-x(i+1, 1, nb) + x(i, 1, nb))
            dylen = (y(i+1, 1, nb) - y(i, 1, nb))
            dlen = dsqrt(dxlen*dxlen + dylen*dylen)
            normalx = dylen/dlen
            normaly = dxlen/dlen

            rho(i, j, nb) = rho(i, 1, nb)
            u(i, j, nb) = u(i, 1, nb)
            v(i, j, nb) = v(i, 1, nb)
            p(i, j, nb) = p(i, 1, nb)
            t(i, j, nb) = p(i, 1, nb)

            un = u(i, j, nb)*normalx + v(i, j, nb)*normaly
            ut = v(i, j, nb)*normalx - u(i, j, nb)*normaly

            un = -un
            ut = ut

            u(i, j, nb) = un*normalx - ut*normaly
            v(i, j, nb) = ut*normalx + un*normaly

            c(i, j, nb) = dsqrt(g1*p(i, j, nb)/rho(i, j, nb))
            h(i, j, nb) = (p(i, j, nb)/g2 + p(i, j, nb) + &
            & half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)) &
            & /rho(i, j, nb)

            q(1, i, j, nb) = rho(i, j, nb)
            q(2, i, j, nb) = rho(i, j, nb)*u(i, j, nb)
            q(3, i, j, nb) = rho(i, j, nb)*v(i, j, nb)
            q(4, i, j, nb) = p(i, j, nb)/g2 + &
                 half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)
         end do
      end do

      ! top boundary conditions
      do j = ny(nb)+1, ny(nb)+gc
         do i = 1, nx(nb)
            dxlen = (-x(i+1, ny(nb), nb) + x(i, ny(nb), nb))
            dylen = (y(i+1, ny(nb), nb) - y(i, ny(nb), nb))
            dlen = dsqrt(dxlen*dxlen + dylen*dylen)
            normalx = dylen/dlen
            normaly = dxlen/dlen

            pamb = 101325.0
            tamb = 300.0
            ramb = pamb/(r*tamb)
            uamb = 0.0d0
            vamb = 0.0d0

            rho(i, j, nb) = rho(i, ny(nb), nb)
            u(i, j, nb) = u(i, ny(nb), nb)
            v(i, j, nb) = v(i, ny(nb), nb)
            p(i, j, nb) = p(i, ny(nb), nb)
            t(i, j, nb) = t(i, ny(nb), nb)

            pdom = p(i, j, nb)
            rdom = rho(i, j, nb)
            udom = u(i, j, nb)
            vdom = v(i, j, nb)

            pprev = p(i, ny(nb)-1, nb)
            rprev = rho(i, ny(nb)-1, nb)
            uprev = u(i, ny(nb)-1, nb)
            vprev = v(i, ny(nb)-1, nb)

            rref = rho(i, j, nb)
            cref = c(i, j, nb)

            qnorm = (udom*normalx + vdom*normaly)
            m = dabs(qnorm)/cref

            if (m < one) then
               if (qnorm < zero) then
                  pb = half*(pdom + pamb - rref*cref*(&
                     & normalx*(uamb - udom) + normaly*(vamb - vdom)))
                  rb = ramb + (pb - pamb)/cref**2
                  ub = uamb - normalx*(pamb - pb)/(rref*cref)
                  vb = vamb - normaly*(pamb - pb)/(rref*cref)
               else
                  pb = pamb
                  rb = rdom + (pb - pamb)/cref**2
                  ub = udom + normalx*(pdom - pamb)/(rref*cref)
                  vb = vdom + normaly*(pdom - pamb)/(rref*cref)
               end if
               p(i, j, nb) = two*pb - pdom
               u(i, j, nb) = two*ub - udom
               v(i, j, nb) = two*vb - vdom
               rho(i, j, nb) = two*rb - rdom
            else
               p(i, j, nb) = two*pdom - pprev
               rho(i, j, nb) = two*rdom - rprev
               u(i, j, nb) = two*udom - uprev
               v(i, j, nb) = two*vdom - vprev
            end if

            c(i, j, nb) = dsqrt(g1*p(i, j, nb)/rho(i, j, nb))
            h(i, j, nb) = (p(i, j, nb)/g2 + p(i, j, nb) + &
            & half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)) &
            & /rho(i, j, nb)
            t(i, j, nb) = p(i, j, nb)/(r*rho(i, j, nb))

            q(1, i, j, nb) = rho(i, j, nb)
            q(2, i, j, nb) = rho(i, j, nb)*u(i, j, nb)
            q(3, i, j, nb) = rho(i, j, nb)*v(i, j, nb)
            q(4, i, j, nb) = p(i, j, nb)/g2 + &
                 half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)
         end do
      end do

      ! right boundary conditions
      do j = 1, ny(nb)
         do i = nx(nb)+1, nx(nb)+gc
            dxlen = (-x(nx(nb), j, nb) + x(nx(nb), j+1, nb))
            dylen = (y(nx(nb), j, nb) - y(nx(nb), j+1, nb))
            dlen = dsqrt(dxlen*dxlen + dylen*dylen)
            normalx = dylen/dlen
            normaly = dxlen/dlen

            rho(i, j, nb) = rho(nx(nb), j, nb)
            p(i, j, nb) = p(nx(nb), j, nb)
            u(i, j, nb) = u(nx(nb), j, nb)
            v(i, j, nb) = v(nx(nb), j, nb)

            c(i, j, nb) = dsqrt(g1*p(i, j, nb)/rho(i, j, nb))
            h(i, j, nb) = (p(i, j, nb)/g2 + p(i, j, nb) + &
            & half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)) &
            & /rho(i, j, nb)
            t(i, j, nb) = p(i, j, nb)/(r*rho(i, j, nb))

            q(1, i, j, nb) = rho(i, j, nb)
            q(2, i, j, nb) = rho(i, j, nb)*u(i, j, nb)
            q(3, i, j, nb) = rho(i, j, nb)*v(i, j, nb)
            q(4, i, j, nb) = p(i, j, nb)/g2 + &
                 half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)
         end do
      end do

      ! left boundary conditions
      ! nozzle exit conditions
      temp = (1 + g2*half*mach_no*mach_no)
      tinf = 300.0d0/temp
      dpres = temp**3.50d0
      pin = npr*101325.0d0/dpres
      rin = pin/(r*tinf)
      uin = mach_no*dsqrt(g1*pin/rin)
      vin = zero
      do j = 1, ny(nb)
         do i = -2, 0
            dxlen = (-x(1, j, nb) + x(1, j+1, nb))
            dylen = (y(1, j, nb) - y(1, j+1, nb))
            dlen = dsqrt(dxlen*dxlen + dylen*dylen)
            normalx = dylen/dlen
            normaly = dxlen/dlen

            if ((y(1, j, nb)-quarter) .le. zero) then
               rho(i, j, nb) = rin
               u(i, j, nb) = uin
               v(i, j, nb) = vin
               p(i, j, nb) = pin
            else
               rho(i, j, nb) = rho(1, j, nb)
               p(i, j, nb) = p(1, j, nb)
               u(i, j, nb) = u(1, j, nb)
               v(i, j, nb) = v(1, j, nb)

               un = u(i, j, nb)*normalx + v(i, j, nb)*normaly
               ut = v(i, j, nb)*normalx - u(i, j, nb)*normaly

               un = -un
               ut = ut

               u(i, j, nb) = un*normalx - ut*normaly
               v(i, j, nb) = ut*normalx + un*normaly
            end if

            c(i, j, nb) = dsqrt(g1*p(i, j, nb)/rho(i, j, nb))
            h(i, j, nb) = (p(i, j, nb)/g2 + p(i, j, nb) + &
            & half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)) &
            & /rho(i, j, nb)

            q(1, i, j, nb) = rho(i, j, nb)
            q(2, i, j, nb) = rho(i, j, nb)*u(i, j, nb)
            q(3, i, j, nb) = rho(i, j, nb)*v(i, j, nb)
            q(4, i, j, nb) = p(i, j, nb)/g2 + &
                 half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)
         end do
      end do

    end subroutine nozzleboundary

 end module boundary
