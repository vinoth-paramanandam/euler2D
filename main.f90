program main
  use constant
  use declaration
  use read_input
  use grid
  use initial
  use boundary


  implicit none

  call read_inputfile

  call meshreader

  call nozzle_initial_conds

end program main
