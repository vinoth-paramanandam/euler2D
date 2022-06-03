!---------------------------------------------------------------------
! Indian Institute of Technology-Madras, PhD scholar
!---------------------------------------------------------------------
!
! MODULE: Constant
!
!> @author
!> Vinoth P}
!
! DESCRIPTION:
!> Constants used to define the precision is provided in this module
!
! REVISION HISTORY:
! 03 June 2022 - Initial Version
!---------------------------------------------------------------------
module constant
  use iso_fortran_env
  implicit none

  ! Set the precision real variables for the program
  integer, parameter :: sp = real32
  integer, parameter :: dp = real64
  ! Set the precision integer variables for the program
  integer, parameter :: i8 = int8   ! -128 to 127
  integer, parameter :: i16 = int16 ! -32768 to 32767
  integer, parameter :: i32 = int32 ! -2147483648 to 2147483647
  integer, parameter :: i64 = int64 ! -9223372036854775808 to 9223372036854775807

  ! Some constants for the program
  real(dp), parameter :: pi = 3.14159265358979311599796346854d0
  real(dp), parameter :: Runiv = 8.3144621d0

  ! Some real numbers for the calculation
  real(dp), parameter :: one = 1.0d0
  real(dp), parameter :: two = 2.0d0
  real(dp), parameter :: three = 3.0d0
  real(dp), parameter :: four = 4.0d0
  real(dp), parameter :: five = 5.0d0
  real(dp), parameter :: six = 6.0d0
  real(dp), parameter :: seven = 7.0d0
  real(dp), parameter :: eight = 8.0d0
  real(dp), parameter :: nine = 9.0d0
  real(dp), parameter :: ten = 10.0d0
  real(dp), parameter :: eleven = 11.0d0
  real(dp), parameter :: twelve = 12.0d0
  real(dp), parameter :: thirteen = 13.0d0

  ! fraction
  real(dp), parameter :: half = 0.5d0
  real(dp), parameter :: quarter = 0.25d0
  real(dp), parameter :: threefourths = 0.75d0
  real(dp), parameter :: sixth = one/six

  ! Some small numbers
  real(dp), parameter :: e_tvd = 1.0d-8

 end module constant
