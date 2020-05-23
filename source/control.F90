!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT

!> \brief deals with reading the input file(s)
!> \details it saves the gathered info in the right variables
!> \author Alin M. Elena (Daresbury Laboratory)
!> \date 11th of April 2016
!> \remarks
!
module m_control
  use iso_fortran_env, only: ERROR_UNIT
  use m_Constants,     only: k_ml,&
                             rp
  use m_Parser,        only: EndParse,&
                             ParseFile,&
                             PrintReport,&
                             getInteger,&
                             getLogical,&
                             getReal,&
                             getString
  use m_Useful,        only: error
  use m_io,            only: SetIO,&
                             ioType

  implicit none
  private

  type, public :: controlType
    real(rp) :: temperature
    real(rp) :: rc
    real(rp) :: lamda
    integer  :: steps
    logical  :: electrostatics
    real(rp) :: alpha
    real(rp) :: delta
    real(rp) :: elecPrecision
    real(rp) :: timestep
    real(rp) :: time = 0.0_rp
    integer  :: step = 0
    integer  :: freq
  end type controlType

  public :: readControl
  !
contains
  !

  subroutine lerror(ou, errno, filename)
    integer, intent(in)             :: ou, errno
    character(len=*), intent(in)    :: filename

    if (errno /= 0) then
      write (ou, *) "I can not open file ", trim(filename)
      write (ou, *) "Error no: ", errno
      stop
    end if

  end subroutine lerror

  subroutine readControl(io, control)
    type(ioType), intent(inout)      :: io
    type(controlType), intent(inout) :: control

    character(len=k_ml) :: dummy
    integer             :: errno = -1, id

    !
    open (newunit=io%uinp, file=trim(io%controlFile), status='old', action='read', iostat=errno)
    call lerror(ERROR_UNIT, errno, io%controlFile)
    ! parse file and in the same time
    ! makes a minimal check of corectness
    call ParseFile(io%uinp)

    io%debugFile = getString(0, "debug", "DEBUG", .false.)
    open (newunit=io%udeb, file=trim(io%debugFile), status='replace', action='write', iostat=errno)
    call lerror(ERROR_UNIT, errno, io%controlFile)
    id = io%udeb
    io%outputFile = getString(id, "output", "OUTPUT", .true.)
    io%trajectoryFile = getString(id, "trajname", "HISTORY", .true.)
    io%stdout = getLogical(id, "l_scr", .false.)
    io%timesFile = getString(id, "stats", "stats.times", .true.)
    io%isTraj = getLogical(id, "trajectory", .false.)
    call setIO(io)

    dummy = getString(id, "title", "no title", .true.)
    io%configFile = getString(id, "config", "CONFIG", .true.)
    io%fieldFile = getString(id, "field", "FIELD", .true.)
    control%temperature = getReal(id, "temperature", 298.0_rp)
    control%steps = getInteger(id, "steps", 42)
    control%rc = getReal(id, "rcut", 12.0_rp)
    control%lamda = getReal(id, "healing", 1.0_rp)
    control%delta = getReal(id, "skin", 3.5_rp)
    control%alpha = getReal(id, "alpha", 0.35_rp)
    control%timestep = getReal(id, "timestep", 0.001_rp)
    control%electrostatics = getLogical(id, "electrostatics", .true.)
    control%elecPrecision = getReal(id, "eprecision", 10e-6_rp)
    control%freq = getInteger(id, "frequency", 1)
    dummy = getString(id, "finish", "", .true.)

    call EndParse(id)
    call PrintReport(id)
    close (io%uinp)

  end subroutine readControl
end module m_control
