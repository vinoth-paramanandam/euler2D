module grid
  use constant
  use declaration
  use iso_fortran_env

  implicit none

  contains

    subroutine meshreader
      implicit none

      integer(i8) :: ifile
      integer(i8) :: nb
      integer(i16) :: nxdomain, nydomain, nzdomain
      integer(i16) :: i, j, imxnb, jmxnb

      character(len=20) :: dummy

      logical :: file_status
      real(real64) :: x1, x2, x3, x4, y1, y2, y3, y4
      real(real64) :: dl1, dl2, dl3, dl4

      nxmax = 0
      nymax = 0
      nblocks = 0

      inquire(file=gridfile, exist=file_status)

      if (file_status) then
         open(unit=13, file=gridfile)
         do
            read(13, *, iostat=ifile) dummy, nxdomain, nydomain, nzdomain

            if (ifile < 0) then
               exit
            elseif (dummy(1:6) == 'domain') then
               nblocks = nblocks + 1
               nxmax = max(nxmax, nxdomain)
               nymax = max(nymax, nydomain)
            end if
         end do
         rewind(13)
      else
            write(error_unit, *) 'Grid file is not found.'
            write(error_unit, *) 'Aborting the program'
            stop
      end if

      ! Allocation of variables
      call mallocate_variables(nxmax, nymax, nblocks)

      ! Read the grid file to get the values
      do nb = 1, nblocks
         read(13, *) dummy, nxdomain, nydomain, nzdomain
         print *, 'Reading the Block ', dummy(8:)
         print *, 'Number of nodes in x direction ', nxdomain
         print *, 'Number of nodes in y direction ', nydomain

         imax(nb) = nxdomain
         jmax(nb) = nydomain
         nx(nb) = nxdomain - 1
         ny(nb) = nydomain - 1

         do j = 1, jmax(nb)
            do i = 1, imax(nb)
               read(13, *) x(i, j, nb), y(i, j, nb)
            end do
         end do
      end do

      do nb = 1, nblocks
         imxnb = imax(nb)
         do j = 1, jmax(nb)
            ! left boundary cells
            x(0, j, nb) = two*x(1, j, nb) - x(2, j, nb)
            x(-1, j, nb) = two*x(0, j, nb) - x(1, j, nb)
            x(-2, j, nb) = two*x(-1, j, nb) - x(0, j, nb)

            y(0, j, nb) = two*y(1, j, nb) - y(2, j, nb)
            y(-1, j, nb) = two*y(0, j, nb) - y(1, j, nb)
            y(-2, j, nb) = two*y(-1, j, nb) - y(0, j, nb)

            ! right boundary cells
            x(imxnb+1, j, nb) = two*x(imxnb, j, nb) - x(imxnb-1, j, nb)
            x(imxnb+2, j, nb) = two*x(imxnb+1, j, nb) - x(imxnb, j, nb)
            x(imxnb+3, j, nb) = two*x(imxnb+2, j, nb) - x(imxnb+1, j, nb)

            y(imxnb+1, j, nb) = two*y(imxnb, j, nb) - y(imxnb-1, j, nb)
            y(imxnb+2, j, nb) = two*y(imxnb+1, j, nb) - y(imxnb, j, nb)
            y(imxnb+3, j, nb) = two*y(imxnb+2, j, nb) - y(imxnb+1, j, nb)
        end do

        jmxnb = jmax(nb)
        do i = 1, imax(nb)
           x(i, 0, nb) = two*x(i, 1, nb) - x(i, 2, nb)
           x(i, -1, nb) = two*x(i, 0, nb) - x(i, 1, nb)
           x(i, -2, nb) = two*x(i, -1, nb) - x(i, 0, nb)

           y(i, 0, nb) = two*y(i, 1, nb) - y(i, 2, nb)
           y(i, -1, nb) = two*y(i, 0, nb) - y(i, 1, nb)
           y(i, -2, nb) = two*y(i, -1, nb) - y(i, 0, nb)

           x(i, jmxnb+1, nb) = two*x(i, jmxnb, nb) - x(i, jmxnb-1, nb)
           x(i, jmxnb+2, nb) = two*x(i, jmxnb+1, nb) - x(i, jmxnb, nb)
           x(i, jmxnb+3, nb) = two*x(i, jmxnb+2, nb) - x(i, jmxnb+1, nb)

           y(i, jmxnb+1, nb) = two*y(i, jmxnb, nb) - y(i, jmxnb-1, nb)
           y(i, jmxnb+2, nb) = two*y(i, jmxnb+1, nb) - y(i, jmxnb, nb)
           y(i, jmxnb+3, nb) = two*y(i, jmxnb+2, nb) - y(i, jmxnb+1, nb)
        end do
      end do

      do nb = 1, nblocks
         do j = 1, ny(nb)
            do i = 1, nx(nb)
               x1 = x(i, j, nb)
               x2 = x(i+1, j, nb)
               x3 = x(i+1, j+1, nb)
               x4 = x(i, j+1, nb)

               y1 = y(i, j, nb)
               y2 = y(i+1, j, nb)
               y3 = y(i+1, j+1, nb)
               y4 = y(i, j+1, nb)

               xc(i, j, nb) = quarter*(x1 + x2 + x3 + x4)
               yc(i, j, nb) = quarter*(y1 + y2 + y3 + y4)

               dl1 = dsqrt((x2 - x1)**2 + (y2 - y1)**2)
               dl2 = dsqrt((x3 - x2)**2 + (y3 - y2)**2)
               dl3 = dsqrt((x4 - x3)**2 + (y4 - y3)**2)
               dl4 = dsqrt((x1 - x4)**2 + (y1 - y4)**2)

               dl(i, j, nb) = dmin1(dl1, dl2, dl3, dl4)
               a(i, j, nb) = half*((x3 - x1)*(y4 - y2) - &
                    (x4 - x2)*(y3 - y1))
            end do
         end do
      end do
    end subroutine meshreader


    subroutine mallocate_variables(ixmax, jymax, nblock)
      implicit none

      integer(i16), intent(in):: ixmax, jymax
      integer(i8), intent(in):: nblock
      integer(i8) :: iostat
      integer(i16) :: istart, iend, jstart, jend
      integer(i8) :: nb

      ! staring and ending points of the grid points

      istart = 1 - gc
      jstart = 1 - gc
      iend = ixmax + gc
      jend = jymax + gc
      nb = nblock

      allocate(imax(nb), jmax(nb), nx(nb), ny(nb))

      allocate(x(istart:iend, jstart:jend, nb), y(istart:iend, jstart:jend, nb), STAT = iostat)
      if (iostat /= 0) STOP 'Error in allocating the points'

      allocate(bctype(istart:iend, jstart:jend, nb), STAT=iostat)
      if (iostat /= 0) stop 'Error in allocating boundary condition type'

      allocate(xc(istart:iend-1, jstart:jend-1, nb), yc(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the cell centers'

      allocate(a(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the areas'

      allocate(dl(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the areas'

      allocate(dt_cell(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the timestep'

      allocate(rho(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the density'

      allocate(u(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the x velocity'

      allocate(v(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the yvelocity'

      allocate(p(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the pressure'

      allocate(t(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the pressure'

      allocate(c(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the sound speed'

      allocate(h(istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the enthalpy'

      allocate(q(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the Conservative variables'

      ! allocate(uleft(4, istart:iend-1, jstart:jend-1, nb), uright(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      ! if (iostat /= 0) stop 'Error in allocating the weno reconstruction variables'

      ! allocate(delLeft(4, istart:iend-1, jstart:jend-1, nb), delRight(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      ! if (iostat /= 0) stop 'Error in allocating the weno reconstruction variables'

      allocate(l_res(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the residual variables'

      allocate(l_vres(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the residual variables'

      allocate(q_rk(4, 4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the timestepping variables'

      allocate(qn(4, istart:iend-1, jstart:jend-1, nb), qi(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the source variables'

      !initialisation of variables
      x = zero
      y = zero
      xc = zero
      yc = zero
      dl = zero
      a = zero
      q = zero
      q_rk = zero
      l_res = zero
      l_vres = zero
      ! delLeft = zero
      ! delRight = zero

      rho = zero
      u = zero
      v = zero
      p = zero
      c = zero
      h = zero
      t = zero

    end subroutine mallocate_variables

    subroutine deallocation
        implicit none

        deallocate(x, y, xc, yc, a, dl)
        deallocate(rho, u, v, p, t, c, h)
        deallocate(q, l_res, l_vres, q_rk, qn, qi)

    end subroutine deallocation
end module grid
