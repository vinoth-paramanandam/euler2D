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
               q(4, i, j, nb) = (cvval*t(i, j, nb) + half*&
                    (u(i, j, nb)**2 + v(i, j, nb)**2))*rho(i, j, nb)
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
                    (u(i, j, nb)**2 + v(i, j, nb)**2))*g2
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
      integer(i8) :: nb, k

      inquire(file = restartfile, exist = file_exist)
      if (file_exist) then
         write(output_unit, *) &
              "Reading the restart fiole for loading the initial conditions"
         open(unit=23, file=restartfile, form='formatted')
         read(23, *) ((((q(k, i, j, nb), k=1, 4), i = 1, nx(nb)), &
           j = 1, ny(nb)), nb =1, nblocks)
         !read(23, *) (((q(1, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         !read(23, *) (((q(2, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         !read(23, *) (((q(3, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         !nread(23, *) (((q(4, i, j, nb), i = 1, nx(nb)), j = 1, ny(nb)), nb = 1, nblocks)
         close(23)
      else
         write(error_unit, *) "Restart file not found"
         write(error_unit, *) " Aborting the program"
         stop
      end if

      call conservative2primitive
    end subroutine read_restart_file

    subroutine timestep
      implicit none

      integer(i16) :: i, j
      integer(i8) :: nb

      real(dp) :: qplus, asound, dt_cellvalue

      if (cfl) then
         dtmin = 1.0e8

         do nb = 1, nblocks
            do j = 1, ny(nb)
               do i = 1, nx(nb)
                  qplus = dsqrt(u(i, j, nb)**2 + v(i, j, nb)**2)
                  asound = c(i, j, nb)
                  dt_cellvalue = cfl_no*dl(i, j, nb)/qplus
                  dt_cell(i, j, nb) = dt_cellvalue
                  dtmin = dmin1(dtmin, dt_cellvalue)
               end do
            end do
         end do
         if (timeaccurate) then
            dt_cell = dtmin
         end if
      else
         dt_cell = dt
         dtmin = dt
      end if
    end subroutine timestep

    subroutine write_output
      integer(i8) :: nb
      integer(i16) :: i, j

      open(unit=unit_id, file=filename, form='formatted')
      nb = 1  ! need to find a fix for this
      write(unit_id, *) "VARIABLES = X, Y, DENSITY, U, V, P, T"
      write(unit_id, *) "ZONE I=", imax(nb), ", J=", jmax(nb), &
           "DATAPACKING=BLOCK, VARLOCATION=([3,4,5,6,7]=CELLCENTERED)"
      write(unit_id, *) "SOLUTIONTIME=", time

      do nb = 1, nblocks
         do j = 1, jmax(nb)
            do i = 1, imax(nb)
               write(unit_id, *) x(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         do j = 1, jmax(nb)
            do i = 1, imax(nb)
               write(unit_id, *) y(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               write(unit_id, *) rho(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               write(unit_id, *) u(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               write(unit_id, *) v(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               write(unit_id, *) p(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               write(unit_id, *) t(i, j, nb)
            end do
         end do
      end do

      close(unit_id)
    end subroutine write_output

    subroutine write_restart
      implicit none
      integer(i8) :: nb, k
      integer(i16) :: i, j

      open(unit=23, status='replace', file=restartfile, form='formatted')
      write(23, *) ((((q(k, i, j, nb), k=1, 4), i = 1, nx(nb)), &
           j = 1, ny(nb)), nb =1, nblocks)
      close(23)
    end subroutine write_restart
  end module misc
