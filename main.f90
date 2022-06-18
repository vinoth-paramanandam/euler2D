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

  if (restart_sim) then
     call read_restart_file
  else
     call nozzle_initial_conds
  end if

  call nozzleboundary



end program main
