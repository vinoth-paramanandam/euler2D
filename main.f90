program main
  use constant
  use declaration
  use read_input
  use grid
  use initial
  use misc
  use boundary
  use odesolver
  use flux
  use iso_fortran_env

  implicit none

  integer(i8) :: nb, k
  integer(i16) :: i, j

  real(dp) :: rss, sum2n, temp
  real(dp) :: diffro, diffmomx, diffmomy, diffenergy
  real(dp), dimension(4) :: sum1n

  logical :: maxcount = .false., l2norm = .false.

  call read_inputfile

  call meshreader

  if (restart_sim) then
     call read_restart_file
  else
     call nozzle_initial_conds
  end if

  call nozzleboundary

  open(unit=29, status='replace', file = outputfile, form='formatted', position='append')

  do
     qi = q

     call timestep

     time = time + dtmin

     do irkstep = 1, nsteps
        l_res = zero

        call isecondMUSCL

        call jsecondMUSCL

        call tvdrk2(irkstep)

        q = qn

        call conservative2primitive

        call nozzleboundary

     end do

     rss = 1.0d8

     diffro = zero
     diffmomx = zero
     diffmomy = zero
     diffenergy = zero

     if (mod(counter, print_val) == 0) then
        if (norm) then
           !$omp parallel do reduction(+:sum1n) &
           !$omp default(private) &
           !$omp shared(nblocks, nx, ny, q, qi)
           do nb = 1, nblocks
              do j = 1, ny(nb)
                 do i = 1, nx(nb)
                    do k = 1, 4
                       sum1n(k) = sum1n(k) + dabs(q(k, i, j, nb) - qi(k, i, j, nb))
                    end do
                 end do
              end do
           end do
           !$omp end parallel do

           do k = 1, 4
              sum1n(k) = dsqrt(sum1n(k))
           end do
           !$omp parallel do reduction(+:sum2n) &
           !$omp default(private) &
           !$omp shared(nblocks, nx, ny, q, qi)
           do nb = 1, nblocks
              do j = 1, ny(nb)
                 do i = 1, nx(nb)
                    temp = q(1, i, j, nb) - qi(1, i, j, nb)
                    sum2n = sum2n + temp*temp
                 end do
              end do
           end do
           !$omp end parallel do

           rss = dsqrt(sum2n)
          ! open(unit = 29, status = 'old', file = outputfile, form = 'formatted')
           write(29, '(1i8, 5es14.7)') counter, time, sum1n(1), &
                sum1n(2), sum1n(3), sum1n(4)
           write(output_unit, '(1i8, 3es14.5)') counter, time, sum1n(1), rss
        end if
     end if

     if (counter >= max_counter) maxcount = .true.
   !   if (rss <= 1.0d-8) l2norm = .true.

     if (maxcount .or. l2norm .or. ftime) then
        write(filename, '(a, i8.8, a)') "plot", counter, ".dat"
        unit_id = unit_id + 100

        call write_output
        call write_restart
        exit
     else
        if (mod(counter, print_iter) == 0) then
           write(filename, '(a, i8.8, a)') "plot", counter, ".dat"
           unit_id = unit_id + 100

           call write_output
           call write_restart
        end if
     end if

     counter = counter + 1
  end do

  call deallocation
  close(29)
end program main
