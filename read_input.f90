module read_input
  use iso_fortran_env
  use constant
  use declaration

  implicit none

  contains

  subroutine read_inputfile
    implicit none

    real(dp) :: isen_ratio
    logical :: file_status

    ! Creating namelists for the use in the input file
    namelist/pbl_params/cfl, ftime, timeaccurate, cfl_no, dt, &
         final_time, max_counter, print_iter, print_val, norm

    namelist/const_params/g1, amol, pr_lam, c1ref

    namelist/file_params/gridfile, restartfile, outputfile, &
         restart_sim, boundaryfile

    namelist/initial_params/mach_no, ismach, uinit, vinit, &
         pinit, tinit

    namelist/total_params/p0, t0

    namelist/specific_conds/problem_specific, npr

    namelist/rk_params/nsteps

    ! Reading the input file
    inquire(file=inputfile, exist=file_status)

    if(file_status) then
       open(unit = 11, file = 'input.in')

       read(11, nml = pbl_params)
       read(11, nml = const_params)
       read(11, nml = file_params)
       read(11, nml = initial_params)
       read(11, nml = specific_conds)
       read(11, nml = rk_params)
       close(11)
    else
       write(error_unit, '(a)') &
            'Input File is not found in the present directory'
       write(error_unit, '(a)') 'Aborting the program'
       stop
    end if

    ! Calculating various constants
    r = Runiv/amol
    g2 = g1 - one
    cpval = g1*r/g2
    cvval = r/g2
    c23 = two/three

    if (ismach) then
       isen_ratio = (one + g2*half*mach_no**2)
       tinit = 300.0d0/isen_ratio
       pinit = npr*101325.0d0/(isen_ratio**3.5d0)
       rinit = pinit/(r*tinit)
       uinit = mach_no*dsqrt(g1*pinit/rinit)
       vinit = zero
    end if

  end subroutine read_inputfile

end module read_input
