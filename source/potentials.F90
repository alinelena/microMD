!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_potentials
  use m_constants, only: k_mw,&
                         rp

  implicit none
  private
  integer, parameter, public :: LJ = 1
  integer, parameter, public :: LJAB = 2
  integer, parameter, public :: LJS = 3
  integer, parameter, public :: GLJ = 4

  type, public :: vdwType
    integer :: i = 0, j = 0, id = 0, np = 0
    character(len=k_mw) :: pot = ''
    character(len=k_mw) :: fpot = ''
    real(rp), allocatable :: param(:)
  end type vdwType
  public :: ljes, ljesS, ljABc, ljfrenkel
contains


  pure subroutine ljfrenkel(vdw,r2,energy,rdU)
    type(vdwType), intent(in)    :: vdw
    real(rp), intent(in)         :: r2
    real(rp), intent(out)        :: energy,rdU

    real(rp) :: ir, e,s2,rc2,rct,st

    rc2 = vdw%param(3)
    if (r2>rc2) then
      energy = 0.0_rp
      rdU = 0.0_rp
    else
      e = vdw%param(1)
      s2 = vdw%param(2)
      ir = 1.0_rp/r2
      st = s2*ir
      rct = rc2*ir
      energy = e*(st-1.0)*(rct-1.0)**2
      rdU = 4.0*e*rct*(rct-1.0)*(st-1.0)+2.0*e*(rct-1.0)**2*st
    end if
  end subroutine ljfrenkel
! watanabe, reinhardt, 1990
  pure real(rp) function Sr(r, rc, l) result(S)
    real(rp), intent(in)    :: r, rc, l

    real(rp) :: ar

    if (r < rc - l) then
      S = 1.0_rp
    elseif (r <= rc) then
      ar = (r - rc + l) / l
      S = 1.0_rp + ar * ar * (2.0_rp * ar - 3.0_rp)
    else
      S = 0.0_rp
    endif
  end function Sr

  pure subroutine ljABc(vdw, r2, U, rdU)
    type(vdwType), intent(in)    :: vdw
    real(rp), intent(in)         :: r2
    real(rp), intent(out)        :: U, rdU

    real(rp) :: a, b, c, s12, s6

    a = vdw%param(1)
    b = vdw%param(2)
    c = vdw%param(3)
    s6 = 1.0_rp / r2; s6 = s6 * s6 * s6
    s12 = s6 * s6

    U = a * s12 - b * c * s6
    rdU = -12.0_rp * (a * s12 - 0.5_rp * c * b * s6)

  end subroutine ljABc

  pure subroutine ljes(vdw, r2, U, rdU)
    type(vdwType), intent(in)    :: vdw
    real(rp), intent(in)         :: r2
    real(rp), intent(out)        :: U, rdU

    real(rp) :: c, e, s, s12, s6

    e = vdw%param(1)
    s = vdw%param(2)
    c = vdw%param(3)

    s6 = s * s / r2; s6 = s6 * s6 * s6
    s12 = s6 * s6
    U = 4.0_rp * e * (s12 - c * s6)
    rdU = 48.0_rp * e * (s12 - 0.5_rp * c * s6)

  end subroutine ljes

  pure subroutine ljesS(vdw, r2, rc, l, U)
    type(vdwType), intent(in)    :: vdw
    real(rp), intent(in)         :: r2, rc, l
    real(rp), intent(out)        :: U

    real(rp) :: c, e, s, s6

    e = vdw%param(1)
    s = vdw%param(2)
    c = vdw%param(3)

    s6 = s * s / r2; s6 = s6 * s6 * s6
    U = 4.0_rp * e * (s6 * s6 - c * s6) * Sr(sqrt(r2), rc, l)
  end subroutine ljesS
end module m_potentials
