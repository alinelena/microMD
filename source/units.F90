!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_units
  use m_constants, only: engUnits,&
                         rp

  implicit none
  private
  ! unit of time      (to)    = 1.000000000 x 10**(-12) seconds   (pico-seconds)
  ! unit of length    (lo)    = 1.000000000 x 10**(-10) metres    (Angstroms)
  ! unit of mass      (mo)    = 1.660540200 x 10**(-27) kilograms (Daltons)
  ! unit of charge    (qo)    = 1.602177330 x 10**(-19) coulombs  (electrons)
  ! unit of energy    (eo)    = 1.660540200 x 10**(-23) joules    (10 J mol^-1)
  ! unit of pressure  (po)    = 1.660540200 x 10**(  7) pascals
  !                           = 1.638825760 x 10**(  2) atmospheres
  !
  ! multiply energy to get internal units.
  ! they are the same internals as DL_POLY
  public :: setEnergyUnits
contains

  subroutine setEnergyUnits(units)
    character(len=*), intent(in)    :: units

    select case (units)
    case ("kcal")
      engUnits = 418.4_rp
    case ("kj")
      engUnits = 100.0_rp
    case default
      engUnits = 1.0_rp
    end select
  end subroutine setEnergyUnits

end module m_units
