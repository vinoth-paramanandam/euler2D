!---------------------------------------------------------------------
! Indian Institute of Technology - Madras, PhD scholar
!---------------------------------------------------------------------
!
! MODULE:  Declaration
!
!> @author
!> Vinoth P}
!
! DESCRIPTION:
!>  Declaration of global variables for the program
!
! REVISION HISTORY:
! 03 June 2022 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!---------------------------------------------------------------------
module declaration
  use constant
  implicit none

  ! Problem initialisation parameters
  logical :: cfl, ftime
  real(dp) :: cfl_no, dt, final_time
  integer(i32) :: max_counter, print_iter, print_val

  ! Gas dynamics constants
  real(dp) :: g1, amol, pr_lam, c1ref
  real(dp) :: g2, r, c23, cpval, cvval

  ! Files for the problem
  character(len=50) :: gridfile, restartfile, outputfile, boundaryfile
  character(len=50) :: initialfile
  character(len=50) :: inputfile='input.in'
  logical :: restart_sim

  ! Counters initialised for the problem
  integer(i32) :: counter = 1
  integer(i8), parameter :: gc = 3

  ! variables for storing the grid
  integer(i16), allocatable, dimension(:) :: imax, jmax
  integer(i16), allocatable, dimension(:) :: nx, ny
  integer(i16) :: nxmax, nymax
  integer(i8) :: nblocks = 0

  real(dp), allocatable, dimension(:, :, :) :: x, y
  real(dp), allocatable, dimension(:, :, :) :: xc, yc
  real(dp), allocatable, dimension(:, :, :) :: a, dl

  ! Some initial conditions variable for the problem
  real(dp) :: mach_no, uinit, vinit, pinit, tinit, rinit, einit, ainit
  logical :: ismach

  ! Problem specific conditions
  real(dp) :: npr
  logical :: problem_specific

  ! fluid properties
  real(dp), allocatable, dimension(:, :, :) :: rho, u, v, p
  real(dp), allocatable, dimension(:, :, :) :: c, h, t

  ! Conservative varaibles
  real(dp), allocatable, dimension(:, :, :, :) :: q

  ! Runge Kutta Time Stepping
  real(dp), allocatable, dimension(:, :, :, :, :) :: q_rk
  real(dp), allocatable, dimension(:, :, :, :) :: l_res, l_vres

  ! Weno reconstruction variables
  ! real(dp), allocatable, dimension(:, :, :, :) :: uleft, uright

  ! File pointers and random filenames
  character(len=50) :: filename
  integer(i16) :: unit_id = 31

  ! Time variable storage
  real(dp), allocatable, dimension(:, :, :) :: dt_cell
  real(dp) :: time, dtmin
  logical :: timeaccurate

  ! Source term
  real(dp), allocatable, dimension(:, :, :, :) :: qn, qi

  ! RK time marching step
  integer(i8) :: irkstep, nsteps

  ! Wall pressure
  real(dp) :: pwall

  ! Boundary type and condition storage
  integer(i32), allocatable, dimension(:, :, :) :: bctype
  integer(i32), allocatable, dimension(:, :, :) :: connection
  integer(i32), allocatable, dimension(:, :, :) :: boundarycells

  ! Total conditions and norm variable
  real(dp) :: t0, p0
  logical :: norm
 end module
