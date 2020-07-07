!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_neighbours
  use m_constants, only: rp
  use m_particles, only: particlesType
  use m_useful,    only: hs

  implicit none
  private

  public :: setupNeighbours
  public :: buildNeighbours

contains

  subroutine setupNeighbours(ps)
    type(particlesType), intent(inout) :: ps

    allocate (ps%neigh(ps%nGParticles))

  end subroutine setupNeighbours

  subroutine buildNeighbours(ps, rc)
    type(particlesType), intent(inout) :: ps
    real(rp), intent(in)               :: rc

    integer  :: i, j
    real(rp) :: r2, rc2, rij(3), si(3),r(3)
    integer  :: mn(ps%nGParticles)

    rc2 = rc**2
    do i = 1, ps%nGParticles - 1
      r = [ps%x(i), ps%y(i), ps%z(i)]
      si = hs(ps%hi, r)
      do j = i + 1, ps%nGParticles
        rij = ps%mic(j, si)
        r2 = sum(rij * rij)
        if (r2 <= rc2) then
          ps%neigh(i)%n = ps%neigh(i)%n + 1
          mn(ps%neigh(i)%n) = j
        endif
      enddo
      if (.not. allocated(ps%neigh(i)%list)) allocate (ps%neigh(i)%list(ps%neigh(i)%n))
      if (ps%neigh(i)%n > size(ps%neigh(i)%list)) then
        deallocate (ps%neigh(i)%list)
        allocate (ps%neigh(i)%list(ps%neigh(i)%n))
      endif
      ps%neigh(i)%list = mn(1:ps%neigh(i)%n)
    enddo
  end subroutine buildNeighbours

end module m_neighbours
