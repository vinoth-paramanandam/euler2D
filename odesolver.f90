module odesolver
  use constant
  use declaration
  use iso_fortran_env

  implicit none

  contains

    subroutine tvdrk2(rkstep)
      implicit none

      integer(i8), intent(in) :: rkstep
      integer(i8) :: nb, k
      integer(i16) :: i, j

      real(dp) :: dtmod

      if(rkstep .eq. 1) then
         do nb = 1, nblocks
            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  dtmod = dt_cell(i, j, nb)/a(i, j, nb)
                  do k = 1, 4
                     q_rk(1, k, i, j, nb) = dtmod*(-l_res(k, i, j, nb))
                     qn(k, i, j, nb) = qi(k, i, j, nb) + &
                          dtmod*(-l_res(k, i, j, nb))
                  end do
               end do
            end do
         end do
      else
         do nb = 1, nblocks
            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  dtmod = dt_cell(i, j, nb)/a(i, j, nb)
                  do k = 1, 4
                     q_rk(2, k, i, j, nb) = dtmod*(-l_res(k, i, j, nb))
                     qn(k, i, j, nb) = qi(k, i, j, nb) + &
                          half*q_rk(1, k, i, j, nb) + &
                          half*q_rk(2, k, i, j, nb)
                  end do
               end do
            end do
         end do
      end if
    end subroutine tvdrk2

  end module odesolver
