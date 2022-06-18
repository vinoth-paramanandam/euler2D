module misc
  use constant
  use declaration

  implicit none

  contains

    subroutine primitive2conservative
      implicit none

      integer(i8) :: nb
      integer(i16) :: i, j

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               q(1, i, j, nb) = rho(i, j, nb)
               q(2, i, j, nb) = rho(i, j, nb)*u(i, j, nb)
               q(3, i, j, nb) = rho(i, j, nb)*v(i, j, nb)
               q(4, i, j, nb) = cvval*t(i, j, nb) + half*rho(i, j, nb)*&
                    & (u(i, j, nb)**2 + v(i, j, nb)**2)
            end do
         end do
      end do

    end subroutine primitive2conservative

    subroutine conservative2primitive
      implicit none

      integer(i8) :: nb
      integer(i16) :: i, j

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               rho(i, j, nb) = q(1, i, j, nb)
               u(i, j, nb) = q(2, i, j, nb)/rho(i, j, nb)
               v(i, j, nb) = q(3, i, j, nb)/rho(i, j, nb)
               p(i, j, nb) = (q(4, i, j, nb) - half*rho(i, j, nb)*&
                    & (u(i, j, nb)**2 + v(i, j, nb)**2))*g2
               c(i, j, nb) = dsqrt(g1*p(i, j, nb)/rho(i, j, nb))
               t(i, j, nb) = p(i, j, nb)/(r*rho(i, j, nb))
               h(i, j, nb) = (p(i, j, nb)/g2 + p(i, j, nb) + &
            & half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)) &
            & /rho(i, j, nb)
            end do
         end do
      end do
    end subroutine conservative2primitive

    subroutine primitive_calc
      implicit none

      integer(i8) :: nb
      integer(i16) :: i, j

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               c(i, j, nb) = dsqrt(g1*p(i, j, nb)/rho(i, j, nb))
               t(i, j, nb) = p(i, j, nb)/(r*rho(i, j, nb))
               h(i, j, nb) = (p(i, j, nb)/g2 + p(i, j, nb) + &
            & half*rho(i, j, nb)*(u(i, j, nb)**2 + v(i, j, nb)**2)) &
            & /rho(i, j, nb)
            end do
         end do
      end do

    end subroutine primitive_calc

  end module misc
