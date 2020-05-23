!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_io
  use iso_fortran_env, only: OUTPUT_UNIT
  use m_constants,     only: k_ml
  use m_useful,        only: error

  implicit none

  private
  type, public :: ioType
    character(len=k_ml) :: controlFile = "flow" !<  name of the input file from which the data is read
!>  name of the file in which the output is put
    character(len=k_ml) :: outputFile
!>   name of the file in which the debug data is written
    character(len=k_ml) :: debugFile
!>   name of the file in which the animation data is written
    character(len=k_ml) :: trajectoryFile = "HISTORY"
    character(len=k_ml) :: timesFile
    character(len=k_ml) :: configFile
    character(len=k_ml) :: fieldFile
!>  how much information to be Printed in the output file
    integer :: verbosity
!>  the debug level
    integer :: debug
!>  output on screen?
    logical :: stdout
!>  do we have trajectory
    logical :: istraj
!>  unit number for input file
    integer :: uinp
!>  unit number for output file
    integer :: uout
!>  unit number for debug file
    integer :: udeb = -1
!>  unit number for animation file
    integer :: utraj
    integer :: ufield
    integer :: utimeser
    integer :: uconfig
!> is it first time when is read?
    logical :: firstTime = .false.
  end type ioType

  public :: setIO
  public :: closeIO
contains

  subroutine setIO(io)
    type(ioType), intent(inout) :: io

    integer :: errno

    if (io%stdout) then
      io%uout = OUTPUT_UNIT
    else
      open (newunit=io%uout, file=trim(io%outputFile), status="unknown", action="write", iostat=errno)
      if (errno /= 0) then
        call error(__FILE__, &
                   __LINE__, &
                   -107, "Cannot open "//trim(io%outputFile)//" file")
      endif
    endif
    if (io%isTraj) then
      open (newunit=io%utraj, file=trim(io%trajectoryFile), status="unknown", action="write", iostat=errno)
      if (errno /= 0) then
        call error(__FILE__, &
                   __LINE__ &
                   , -108, "Cannot open "//trim(io%trajectoryFile)//" file")
      endif
    endif
      open (newunit=io%utimeser, file=trim(io%timesFile), status="unknown", action="write", iostat=errno)
      if (errno /= 0) then
        call error(__FILE__, &
                   __LINE__ &
                   , -109, "Cannot open "//trim(io%timesFile)//" file")
      endif

  end subroutine setIO

  subroutine closeIO(io)
    type(ioType), intent(inout) :: io

    if (io%isTraj) then
      close (io%utraj)
    endif
    close (io%utimeser)
  end subroutine closeIO
end module m_io
