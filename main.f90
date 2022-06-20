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

  implicit none

  integer(i8) :: nb, k
  integer(i16) :: i, j

  real(dp) :: rss, sum2n, temp
  real(dp) :: diffro, diffmomx, diffmomy, diffenergy
  real(dp), dimension(4) :: sum1n

  logical :: maxcount, l2norm

  call read_inputfile

  call meshreader

  if (restart_sim) then
     call read_restart_file
  else
     call nozzle_initial_conds
  end if

  call nozzleboundary

  do
     qi = q

     call timestep

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
           do nb = 1, nblocks
              do j = 1, ny(nb)
                 do i = 1, nx(nb)
                    do k = 1, 4
                       sum1n(k) = sum1n(k) + dabs(q(k, i, j, nb) - qi(k, i, j, nb))
                    end do
                 end do
              end do
           end do

           do k = 1, 4
              sum1n(k) = dsqrt(sum1n(k))
           end do

           do nb = 1, nblocks
              do j = 1, ny(nb)
                 do i = 1, nx(nb)
                    temp = q(1, i, j, nb) - qi(1, i, j, k)
                    sum2n = sum2n + temp*temp
                 end do
              end do
           end do

           rss = dsqrt(sum2n)
           open(unit = 29, status = 'old', file = outputfile, form = 'formatted')
           write(29, '(1i16, 5es14.7)') counter, time, sum1n(1), &
                sum1n(2), sum1n(3), sum1n(4)
           write(output_unit, '(1i16, 3es14.5)') counter, time, sum1n, rss
        end if
     end if

     time = time + dtmin

     if (counter >= max_counter) maxcount = .true.
     if (rss <= 1.0d-8) l2norm = .true.

     if (maxcount .or. l2norm .or. ftime) then
        write(filename, '(a, i8.8, a)') "plot", counter, ".dat"
        unit_id = counter + 31

        call write_output
        call write_restart
        exit
     else
        if (mod(counter, print_iter) == 0) then
           write(filename, '(a, i8.8, a)') "plot", counter, ".dat"
           unit_id = counter + 31

           call write_output
           call write_restart
        end if
     end if

  end do

  call deallocation
end program main
