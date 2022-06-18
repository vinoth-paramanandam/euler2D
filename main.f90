program main
  use constant
  use declaration
  use read_input
  use grid
  use initial
  use misc
  use boundary


  implicit none

  call read_inputfile

  call meshreader

  call nozzle_initial_conds

  call primitive_calc

  call nozzleboundary

  call primitive2conservative



end program main
