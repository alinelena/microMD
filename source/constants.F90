!    Copyright (c) 2006-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT

!> \brief defines constants that are needed in the code.
!> \author Alin M Elena
!> \date January 2007
module m_Constants
  implicit none
  private
!
  integer, parameter, public :: k_ml = 255 !< number of character per line
  integer, parameter, public :: k_mw = 42 !< number of character per word
  integer, parameter, public :: rp = kind(1.0d0)
  integer, parameter, public :: ip = 4
  integer, parameter, public :: elLen = 6

  ! this is not a computer constant but a problem constant
  real(rp), public :: engUnits = 1.0_rp
  real(rp), parameter, public :: pi = 4.0_rp * atan(1.0_rp)
  real(rp), parameter, public :: fpi = 4.0_rp * pi
  real(rp), parameter, public :: tpi = 2.0_rp * pi
  real(rp), parameter, public :: ifpi = 0.25_rp / pi
  real(rp), parameter, public :: kB = 8.31451115e-1_rp
end module m_Constants

