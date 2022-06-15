module boundary.f90
  use constants
  use declaration
  use grid

  implicit none

  contains

    subroutine nozzleboundary
      implicit none

      integer(i8) :: nb, i, j
      real(dp) :: dxlen, dylen, dlen, normalx, normaly
      real(dp) :: un, ut, pamb, uamb, vamb, tamb

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
            t(i, j, nb) = p(i, ny(nb), nb)
            
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
         end do
      end do

      ! left boundary conditions
      do j = 1, ny(nb)
         do i = -2, 0
            dxlen = (-x(1, j, nb) + x(1, j+1, nb))
            dylen = (y(1, j, nb) - y(1, j+1, nb))
            dlen = dsqrt(dxlen*dxlen + dylen*dylen)
            normalx = dylen/dlen
            normaly = dxlen/dlen

         end do
      end do

    end subroutine nozzleboundary

 end module boundary
