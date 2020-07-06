!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_potentials
  use m_constants, only: k_mw,&
                         rp

  implicit none

  private
  integer, parameter, public :: LJ = 1
  integer, parameter, public :: LJ_AB = 2
  integer, parameter, public :: LJ_S = 3
  integer, parameter, public :: GLJ = 4
  integer, parameter, public :: LJ_C = 5

  character(len=*), parameter :: pot_labels(5)=[" Lennard Jones sig-eps",&
                                             "      Lennard Jones AB",&
                                             "  Lennard Jones smooth",&
                                             " Lennard Jones Frenkel",&
                                             "Lennard-Jones Cehesive"]

 type, abstract, public :: interaction
    integer :: j, i, id
  contains
    procedure(initialisation), deferred :: init
    procedure(potential), deferred :: pot
    procedure(lrc), deferred :: lrc
  end type

  abstract interface
    subroutine initialisation(t, p, i, j, id)
      import interaction, rp
      class(interaction), intent(inout) :: t
      real(rp), intent(in) :: p(:)
      integer, intent(in) :: i,j,id

    end subroutine initialisation

    subroutine potential(t, r2, U, rdU)
      import interaction, rp
      class(interaction), intent(inout) :: t
      real(rp), intent(in) :: r2
      real(rp), intent(out) :: U, rdU
    end subroutine potential

    subroutine lrc(t, rc, u, du)
      import interaction, rp
      class(interaction), intent(inout) :: t
      real(rp), intent(in) :: rc
      real(rp), intent(out) :: U, dU
    end subroutine lrc
  end interface

  type, public :: pot_holder
    class(interaction), allocatable :: h
  end type

  type, extends(interaction), public :: ljse
    real(rp) :: sig, eps
  contains
    procedure :: init => ljse_init
    procedure :: pot => ljse_i
    procedure :: lrc => ljse_lrc
  end type

  type, extends(interaction), public :: ljAB
    real(rp) :: a, b
  contains
    procedure :: init => ljab_init
    procedure :: pot => ljab_i
    procedure :: lrc => ljab_lrc
  end type

  type, extends(interaction), public :: ljsec
    real(rp) :: sig, eps, c
  contains
    procedure :: init => ljsec_init
    procedure :: pot => ljsec_i
    procedure :: lrc => ljsec_lrc
  end type

  type, extends(interaction), public :: ljf
    real(rp) :: sig, eps, rcut2
  contains
    procedure :: init => ljf_init
    procedure :: pot => ljf_i
    procedure :: lrc => ljf_lrc
  end type

  type, extends(interaction), public :: ljs
    real(rp) :: sig, eps, rc,l
  contains
    procedure :: init => ljs_init
    procedure :: pot => ljs_i
    procedure :: lrc => ljs_lrc
  end type

  contains

    subroutine ljse_init(t, p, i,j,id)
    class(ljse), intent(inout) :: t
    real(rp), intent(in)       :: p(:)
      integer, intent(in) :: i,j,id

    t%eps = p(1)
    t%sig = p(2)
    t%i = i
    t%j = j
    t%id = id
  end subroutine ljse_init

  pure subroutine ljse_i(t, r2, U, rdU)
    class(ljse), intent(inout) :: t
    real(rp), intent(in)       :: r2
    real(rp), intent(out)      :: U, rdU

    real(rp) :: s12, s6

    s6 = t%sig * t%sig / r2; s6 = s6 * s6 * s6
    s12 = s6 * s6
    U = 4.0_rp * t%eps * (s12 - s6)
    rdU = 48.0_rp * t%eps * (s12 - 0.5_rp * s6)

  end subroutine ljse_i

  pure subroutine ljse_lrc(t, rc, U, dU)
    class(ljse), intent(inout) :: t
    real(rp), intent(in)       :: rc
    real(rp), intent(out)      :: U, dU

    U = 4.0_rp*T%eps*(t%sig**12-3.0_rp*rc**6*t%sig**6)/(rc**9*9.0_rp)
    dU = 0.0_rp

  end subroutine ljse_lrc

  pure subroutine ljsec_init(t, p, i, j,id)
    class(ljsec), intent(inout) :: t
    real(rp), intent(in)      :: p(:)
    integer, intent(in)       :: i,j,id

    t%eps = p(1)
    t%sig = p(2)
    t%c = p(3)
    t%i = i
    t%j = j
    t%id = id
  end subroutine ljsec_init

  pure subroutine ljsec_i(t, r2, U, rdU)
    class(ljsec), intent(inout) :: t
    real(rp), intent(in)       :: r2
    real(rp), intent(out)      :: U, rdU

    real(rp) :: s12, s6

    s6 = t%sig * t%sig / r2; s6 = s6 * s6 * s6
    s12 = s6 * s6
    U = 4.0_rp * t%eps * (s12 - t%c*s6)
    rdU = 48.0_rp * t%eps * (s12 - 0.5_rp * t%c* s6)

  end subroutine ljsec_i

  pure subroutine ljsec_lrc(t, rc, U, dU)
    class(ljsec), intent(inout) :: t
    real(rp), intent(in)       :: rc
    real(rp), intent(out)      :: U, dU

    U = 4.0_rp*T%eps*(t%sig**12-3.0_rp*t%c*rc**6*t%sig**6)/(rc**9*9.0_rp)
    dU = 0.0_rp

  end subroutine ljsec_lrc





  pure subroutine ljf_init(t, p, i, j,id)
    class(ljf), intent(inout) :: t
    real(rp), intent(in)      :: p(:)
    integer, intent(in)       :: i,j,id

    real(rp) :: x

    t%sig = p(2)*p(2)
    t%rcut2 = p(3)*p(3)
    x = (t%rcut2/t%sig)
    t%eps = p(1) * 2.0_rp*x*(1.5_rp/(x-1.0_rp))**3
    t%i = i
    t%j = j
    t%id = id
  end subroutine ljf_init

  pure subroutine ljf_i(t, r2, U, rdU)
    class(ljf), intent(inout) :: t
    real(rp), intent(in)      :: r2
    real(rp), intent(out)     :: U, rdU

    real(rp) :: ir, rct, st

    if (r2 > t%rcut2) then
      U = 0.0_rp
      rdU = 0.0_rp
    else
      ir = 1.0_rp / r2
      st = t%sig * t%sig * ir
      rct = t%rcut2 * ir
      U = t%eps * (st - 1.0_rp) * (rct - 1.0_rp)**2
      rdU = 4.0_rp * t%eps * rct * (rct - 1.0_rp) * (st - 1.0_rp) + 2.0_rp * t%eps * (rct - 1.0_rp)**2 * st
    end if
  end subroutine ljf_i

  pure subroutine ljf_lrc(t, rc, U, dU)
    class(ljf), intent(inout) :: t
    real(rp), intent(in)       :: rc
    real(rp), intent(out)      :: U, dU

    U = 0.0_rp
    dU = 0.0_rp

  end subroutine ljf_lrc

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

    subroutine ljab_init(t, p, i,j,id)
    class(ljab), intent(inout) :: t
    real(rp), intent(in)       :: p(:)
      integer, intent(in) :: i,j,id

    t%a = p(1)
    t%b = p(2)
    t%i = i
    t%j = j
    t%id = id
  end subroutine ljab_init


  pure subroutine ljAB_i(t, r2, U, rdU)
    class(ljAB), intent(inout)    :: t
    real(rp), intent(in)         :: r2
    real(rp), intent(out)        :: U, rdU

    real(rp) ::  s12, s6

    s6 = 1.0_rp / r2; s6 = s6 * s6 * s6
    s12 = s6 * s6

    U = t%a * s12 - t%b * s6
    rdU = -12.0_rp * (t%a * s12 - 0.5_rp * t%b * s6)

  end subroutine ljAB_i

  pure subroutine ljAB_lrc(t, rc, U, dU)
    class(ljAB), intent(inout) :: t
    real(rp), intent(in)       :: rc
    real(rp), intent(out)      :: U, dU

    U = 0.0_rp
    dU = 0.0_rp

  end subroutine ljAB_lrc

    subroutine ljs_init(t, p, i,j,id)
    class(ljs), intent(inout) :: t
    real(rp), intent(in)       :: p(:)
      integer, intent(in) :: i,j,id

    t%eps = p(1)

    t%sig = p(2)
    t%rc = p(3)
    t%l  = p(4)
    t%i = i
    t%j = j
    t%id = id
  end subroutine ljs_init

  pure subroutine ljs_i(t, r2, U, rdU)
    class(ljs), intent(inout)    :: t
    real(rp), intent(in)         :: r2
    real(rp), intent(out)        :: U,rdU

    real(rp) :: s, s6


    s6 = t%sig * t%sig / r2; s6 = s6 * s6 * s6
    U = 4.0_rp * t%eps * (s6 * s6 - s6) * Sr(sqrt(r2), t%rc, t%l)
    rdU = 0.0_rp
  end subroutine ljs_i

  pure subroutine ljs_lrc(t, rc, U, dU)
    class(ljs), intent(inout) :: t
    real(rp), intent(in)       :: rc
    real(rp), intent(out)      :: U, dU

    U = 0.0_rp
    dU = 0.0_rp

  end subroutine ljs_lrc

end module m_potentials
