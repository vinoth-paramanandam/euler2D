module misc
  use constant
  use declaration
  use iso_fortran_env

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

    subroutine read_restart_file
      implicit none

      logical :: file_exist
      integer(i16) :: i, j
      integer(i8) :: nb

      inquire(file = restartfile, exist = file_exist)
      if (file_exist) then
         write(output_unit, *) &
              "Reading the restart fiole for loading the initial conditions"
         open(unit=23, file=restartfile, form='formatted')
         read(23, *) (((q(1, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         read(23, *) (((q(2, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         read(23, *) (((q(3, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         read(23, *) (((q(4, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         close(23)
      else
         write(error_unit, *) "Restart file not found"
         write(error_unit, *) " Aborting the program"
         stop
      end if

      call conservative2primitive
    end subroutine read_restart_file
  end module misc
