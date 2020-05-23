!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT

module m_sampler
  use m_constants, only: rp
  use m_particles, only: particlesType

  implicit none

  private
  public :: nve_update_velocity
  public :: nve_update_position
contains

  subroutine nve_update_velocity(ps, dt)
    type(particlesType), intent(inout) :: ps
    real(rp), intent(in)               :: dt

    integer  :: i, sp
    real(rp) :: im

    do i = 1, ps%nGparticles
      sp = ps%spec(i)
      im = 1.0_rp / ps%mass(sp)
      ps%vx(i) = ps%vx(i) + 0.5_rp * ps%fx(i) * dt * im
      ps%vy(i) = ps%vy(i) + 0.5_rp * ps%fy(i) * dt * im
      ps%vz(i) = ps%vz(i) + 0.5_rp * ps%fz(i) * dt * im
    enddo
  end subroutine nve_update_velocity

  subroutine nve_update_position(ps, dt)
    type(particlesType), intent(inout) :: ps
    real(rp), intent(in)               :: dt

    integer :: i

    do i = 1, ps%nGparticles
      ps%x(i) = ps%x(i) + ps%vx(i) * dt
      ps%y(i) = ps%y(i) + ps%vy(i) * dt
      ps%z(i) = ps%z(i) + ps%vz(i) * dt
    enddo
  end subroutine nve_update_position

end module m_sampler
