!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_forces
  use m_constants, only: ifpi,&
                         pi,&
                         rp,&
                         tpi
  use m_control,   only: controlType
  use m_particles, only: particlesType
  use m_useful,    only: getvdw,&
                         hs

  implicit none
  private
  public :: computeForces
  public :: long_range_correction
contains

  subroutine computeForces(ps, control)
    type(particlesType), intent(inout) :: ps
    type(controlType), intent(in)      :: control

    integer  :: i, is, j, js, k, l
    real(rp) :: eng, r(3), f(3), ir, r2, rdU, rij(3), si(3), U

    ps%fx = 0.0_rp
    ps%fy = 0.0_rp
    ps%fz = 0.0_rp

    !call ewaldForces(ps)
    do i = 1, ps%nGParticles
      is = ps%spec(i)
      r = [ps%x(i), ps%y(i), ps%z(i)]
      si = hs(ps%hi, r)
      eng = 0.0_rp
      f = 0.0_rp
      do l = 1, ps%neigh(i)%n
        j = ps%neigh(i)%list(l)
        rij = ps%mic(j, si)
        js = ps%spec(j)
        k = getVdw(is, js, ps%nSpecies)
        r2 = sum(rij * rij); ir = 1.0_rp / r2
        call ps%pots(k)%h%pot(r2, U, rdU)
        f = rij * rdU * ir
        eng = eng + U
        if (i < j) then
          ps%fx(j) = ps%fx(j) - f(1)
          ps%fy(j) = ps%fy(j) - f(2)
          ps%fz(j) = ps%fz(j) - f(3)
        endif
        ps%fx(i) = ps%fx(i) + f(1)
        ps%fy(i) = ps%fy(i) + f(2)
        ps%fz(i) = ps%fz(i) + f(3)
      enddo
      ps%engStress(i)%engPair = eng
    enddo
  end subroutine computeForces

  real(rp) function self_ewald(ps)
    type(particlesType), intent(inout) :: ps

    integer  :: i
    real(rp) :: q

    q = 0.0_rp
    do i = 1, ps%nGParticles
      q = q + ps%charge(i)**2
    end do
    ps%eeself = sqrt(ps%alpha / pi) * q
  end function self_ewald

  subroutine ewaldForcesRealSpace(ps)
    type(particlesType), intent(inout) :: ps

    integer  :: i, is, j, js, l
    real(rp) :: alpha, eng, f(3), ir, qa, qb, r, rij(3), sa, si(3)

    alpha = 0.3_rp
    sa = sqrt(alpha)
    do i = 1, ps%nGParticles
      is = ps%spec(i)
      qa = ps%charge(is)
      si = hs(ps%hi, [ps%x(i), ps%y(i), ps%z(i)])
      f = 0.0_rp
      eng = 0.0_rp
      do l = 1, ps%neigh(i)%n
        j = ps%neigh(i)%list(l)
        rij = ps%mic(j, si)
        js = ps%spec(j)
        qb = ps%charge(js)
        r = sqrt(sum(rij * rij)); ir = 1.0_rp / r
        eng = eng + qa * qb * ir * ifpi * erfc(sa * r)
      enddo
      ps%engStress(i)%engElec = eng * ifpi
    enddo
  end subroutine ewaldForcesRealSpace

  subroutine ewaldForcesKSpace(ps)
    type(particlesType), intent(inout) :: ps

  end subroutine ewaldForcesKSpace

  subroutine long_range_correction(ps, rc)
    type(particlesType), intent(inout) :: ps
    real(rp), intent(in)               :: rc

    integer  :: i, j, k
    real(rp) :: delrc, elrc

    ps%englrc = 0.0_rp
    do i = 1, ps%nSpecies
      do j = 1, i
        k = getVdw(i, j, ps%nSpecies)
        call ps%pots(k)%h%lrc(rc, elrc, delrc)
        if (i /= j) then
          elrc = elrc * 2.0_rp
        end if
        ps%englrc = ps%englrc + tpi * ps%isp(i) * ps%isp(j) / ps%vol * elrc
      end do
    end do
  end subroutine long_range_correction
end module m_forces
