module initial
  use constant
  use declaration
  use iso_fortran_env

  implicit none

  contains

    subroutine read_initial_conds()
      implicit none

      integer(i16) :: i, j
      integer(i8) :: nb

      logical :: file_status

      inquire(file=initialfile, exist=file_status)
      if (file_status) then

         open(unit=19, file=initialfile)

         do nb = 1, nblocks
            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  read(19, *) rho(i, j, nb)
               end do
            end do

            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  read(19, *) u(i, j, nb)
               end do
            end do

            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  read(19, *) v(i, j, nb)
               end do
            end do

            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  read(19, *) p(i, j, nb)
               end do
            end do

         end do
         close(19)
      else
         write(error_unit, *) 'Initial conditions file is not found'
         write(error_unit, *) 'Aborting the program'
         stop
      end if

    end subroutine read_initial_conds

    subroutine nozzle_initial_conds
      implicit none

      integer(i8) :: nb
      integer(i16) :: i, j

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               rho(i, j, nb) = 1.2250d0
               u(i, j, nb) = zero
               v(i, j, nb) = zero
               p(i, j, nb) = 101325.0d0
            end do
         end do
      end do
    end subroutine nozzle_initial_conds
end module initial
