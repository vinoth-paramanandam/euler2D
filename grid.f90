module grid
  use constant
  use declaration
  use iso_fortran_env

  implicit none

  contains

    subroutine meshreader
      implicit none

      inquire(file=gridfile, exist=file_status)

      if (file_status) then
         do
            read(13, *, iostat=ifile) dummy, nxdomain, nydomain, nzdomain

            if (ifile < 0) then
               exit
            elseif (dummy(1:6) == 'domain') then
               nblocks = nblocks + 1
               nxmax = max(nxmax, nxdomain)
               nymax = max(nymax, nydomain)
            end if
            rewind(13)
         end do
      else
            write(error_unit, 'a') 'Grid file is not found.'
            write(error_unit, 'a') 'Aborting the program'
      end if

      ! Read the grid file to get the values
      do nb = 1, nblocks
         print *, 'Reading the Block ', dummy(8:)
         print *, 'Number of nodes in x direction ', nxdomain
         print *, 'Number of nodes in y direction ', nydomain

         imax(nb) = nxdomain
         jmax(nb) = nydomain
         nx(nb) = nxdomain - 1
         ny(nb) = nydomain - 1

         do j = 1, jmax(nb)
            do i = 1, imax(nb)
               read(3, *) x(i, j, nb), y(i, j, nb)
            end do
         end do
      end do

      close (13)

      call mallocate_variables(nxmax, nymax, nblocks)

    end subroutine meshreader


    subroutine mallocate_variables(ixmax, jymax, nblock)
      implicit none

      integer, intent(in):: ixmax, jymax, nblock
      integer :: iostat
      integer :: istart, iend, jstart, jend
      integer :: nb

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

      allocate(uleft(4, istart:iend-1, jstart:jend-1, nb), uright(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the weno reconstruction variables'

      allocate(delLeft(4, istart:iend-1, jstart:jend-1, nb), delRight(4, istart:iend-1, jstart:jend-1, nb), STAT = iostat)
      if (iostat /= 0) stop 'Error in allocating the weno reconstruction variables'

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
      delLeft = zero
      delRight = zero

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
        deallocate(q, uleft, uright, delLeft, delRight, l_res, l_vres, q_rk, qn, qi)

    end subroutine deallocation
end module grid
