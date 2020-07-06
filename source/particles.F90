!ke    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_particles
  use m_Constants,  only: elLen,&
                          ip,&
                          kB,&
                          k_ml,&
                          k_mw,&
                          rp
  use m_io,         only: ioType
  use m_potentials, only: pot_holder
  use m_useful,     only: gaussian,&
                          hs,&
                          isInSet

  implicit none
  private
  type, public :: neighType
    integer :: n = 0
    integer, allocatable :: list(:)
  end type neighType

  type, public :: engStressType
    real(rp)   :: engPair = 0.0_rp
    real(rp)   :: engElec = 0.0_rp
    real(rp)   :: stress(3, 3) = 0.0_rp
  end type engStressType

  type, public                        :: particlesType
    character(len=k_ml)               :: systemName
    character(len=k_ml)               :: fieldName
    integer                           :: nGParticles = 0
    integer(ip), allocatable          :: GIDs(:), spec(:)
    real(rp), allocatable             :: x(:), y(:), z(:)
    real(rp), allocatable             :: vx(:), vy(:), vz(:)
    real(rp), allocatable             :: fx(:), fy(:), fz(:)
    real(rp)                          :: h(3, 3) = 0.0_rp
    real(rp)                          :: hi(3, 3) = 0.0_rp
    integer                           :: pbc, nSpecies
    character(len=elLen), allocatable :: labels(:)
    character(len=elLen), allocatable :: species(:) ! unique set of species available
    real(rp), allocatable             :: mass(:) ! mass of species
    real(rp), allocatable             :: charge(:) ! charge of species
    real(rp), allocatable             :: colour(:) ! for nemd
    integer, allocatable              :: isp(:) ! count of each specie in the configuration
    type(neighType), allocatable      :: neigh(:) ! neighbours list
    type(engStressType), allocatable  :: engStress(:) ! energy and stress per particle
    type(pot_holder), allocatable      :: pots(:)
    real(rp)                          :: ke
    real(rp)                          :: eng
    real(rp)                          :: engPair
    real(rp)                          :: englrc
    real(rp)                          :: virlrc
    real(rp)                          :: engElec
    real(rp)                          :: temperature
    real(rp)                          :: eeself = 0.0_rp, alpha, nk(3)
    real(rp)                          :: stress(3, 3)
    real(rp)                          :: Vol = 0.0_rp
    real(rp)                          :: density = 0.0_rp
    real(rp)                          :: m = 0.0_rp
    real(rp)                          :: q = 0.0_rp
    real(rp)                          :: com(3)
    real(rp)                          :: cog(3)
    character(len=k_mw)               :: units
  contains
    private
    procedure, public                 :: summary => statsOnConfig
    procedure, public                 :: fieldSummary => statsOnField
    procedure, public                 :: getMass
    procedure, public                 :: getCharge
    procedure, public                 :: volume
    procedure, public                 :: energy
    procedure, public                 :: energyKinetic
    procedure, public                 :: centreOfMass
    procedure, public                 :: mic
    procedure, public                 :: writeTrajectory
    procedure, public                 :: init_velocities
    procedure, public                 :: scale_velocities
    final                             :: cleanup

  end type particlesType

contains

  subroutine cleanup(ps)
    type(particlesType) :: ps

    integer :: i

    if (allocated(ps%GIDs)) deallocate (ps%GIDs)
    if (allocated(ps%spec)) deallocate (ps%spec)
    if (allocated(ps%x)) deallocate (ps%x)
    if (allocated(ps%y)) deallocate (ps%y)
    if (allocated(ps%z)) deallocate (ps%z)
    if (allocated(ps%vx)) deallocate (ps%vx)
    if (allocated(ps%vy)) deallocate (ps%vy)
    if (allocated(ps%vz)) deallocate (ps%vz)
    if (allocated(ps%fx)) deallocate (ps%fx)
    if (allocated(ps%fy)) deallocate (ps%fy)
    if (allocated(ps%fz)) deallocate (ps%fz)
    if (allocated(ps%labels)) deallocate (ps%labels)
    if (allocated(ps%species)) deallocate (ps%species)
    if (allocated(ps%isp)) deallocate (ps%isp)
    if (allocated(ps%mass)) deallocate (ps%mass)
    if (allocated(ps%charge)) deallocate (ps%charge)
    if (allocated(ps%colour)) deallocate (ps%colour)
    if (allocated(ps%pots)) deallocate (ps%pots)
    do i = 1, size(ps%pots)
      if (allocated(ps%pots(i)%h)) deallocate (ps%pots(i)%h)
    enddo
    if (allocated(ps%pots)) deallocate (ps%pots)
    do i = 1, size(ps%neigh)
      if (allocated(ps%neigh(i)%list)) deallocate (ps%neigh(i)%list)
    enddo
    if (allocated(ps%neigh)) deallocate (ps%neigh)
    if (allocated(ps%engStress)) deallocate (ps%engStress)
  end subroutine cleanup

  subroutine statsOnConfig(ps, io)
    class(particlesType)        :: ps
    type(ioType), intent(in)    :: io

    integer :: i, k
    logical :: b

    write (io%uout, '(a,i0)') "Number of atoms in "//trim(io%configFile)//": ", ps%nGParticles
    call arrayToSet(ps)
    do i = 1, ps%nGParticles
      b = isInSet(ps%labels(i), ps%species, ps%nSpecies, k)
      ps%spec(i) = k
    enddo
    write (io%uout, '(a,i0)') "number of species: ", ps%nSpecies
    write (io%uout, '(a)', advance="no") "Composition: "
    do i = 1, ps%nSpecies
      write (io%uout, "(a,i0)", advance="no") trim(ps%species(i))//"_", ps%isp(i)
    enddo
    write (io%uout, *)
    write (io%uout, '(a,f16.8)') "Volume: ", ps%volume()
  end subroutine statsOnConfig

  function centreOfMass(ps)
    class(particlesType) :: ps
    real(rp)             :: centreOfMass(3)

    integer  :: i, sp
    real(rp) :: r(3)

    r = 0.0_rp
    do i = 1, ps%nGParticles
      sp = ps%spec(i)
      r = r + ps%mass(sp)*[ps%x(i), ps%y(i), ps%z(i)]
    enddo
    centreOfMass = r / ps%m
    ps%com = centreOfMass

  end function centreOfMass

  real(rp) function getMass(ps)
    class(particlesType) :: ps

    integer  :: i, sp
    real(rp) :: m

    m = 0.0_rp
    do i = 1, ps%nGParticles
      sp = ps%spec(i)
      m = m + ps%mass(sp)
    end do
    ps%m = m
    getMass = m
  end function getMass

  real(rp) function getCharge(ps)
    class(particlesType) :: ps

    integer  :: i, sp
    real(rp) :: c

    c = 0.0_rp
    do i = 1, ps%nGParticles
      sp = ps%spec(i)
      c = c + ps%charge(sp)
    end do
    ps%q = c
    getCharge = c
  end function getCharge

  function mic(ps, j, si)
    class(particlesType)    :: ps
    integer, intent(in)     :: j
    real(rp), intent(in)    :: si(3)
    real(rp), dimension(3)  :: mic

    real(rp) :: sij(3), sj(3)

    sj = hs(ps%hi, [ps%x(j), ps%y(j), ps%z(j)])
    sij = si - sj
    sij = sij - nint(sij)
    mic = hs(ps%h, sij)
  end function mic

  real(rp) function volume(ps)
    class(particlesType) :: ps

    volume = ps%h(1, 1) * (ps%h(2, 2) * ps%h(3, 3) - ps%h(2, 3) * ps%h(3, 2)) - &
             ps%h(1, 2) * (ps%h(2, 1) * ps%h(3, 3) - ps%h(2, 3) * ps%h(3, 1)) + &
             ps%h(1, 3) * (ps%h(2, 1) * ps%h(3, 2) - ps%h(2, 2) * ps%h(3, 1))
    ps%Vol = volume

  end function volume

  subroutine statsOnField(ps, io)
    class(particlesType)        :: ps
    type(ioType), intent(in)    :: io

    integer :: i

    write (io%uout, "(a8,3(a12))") "Specie|", "m|", "q|", "colour|"
    do i = 1, ps%nSpecies
      write (io%uout, "(a8,3(g11.4,1x))") trim(ps%species(i)), ps%mass(i), ps%charge(i), &
        ps%colour(i)
    enddo
    write (io%uout, '(a)') "Two body interactions:"
    write (io%uout, '(a8,a8,a5,a21,a36)') " A|", "B|", repeat(" ", 5), "Potential|", "Params|"
    do i = 1, size(ps%pots)
      write (io%uout, '(a8,a8,a5,a24,1x,4(g11.4,1x))') ps%species(ps%pots(i)%h%i), ps%species(ps%pots(i)%h%j)
    enddo
    write (io%uout, '(a,f16.8)') "Mass: ", ps%getMass()
    write (io%uout, '(a,f16.8)') "Density(ρ): ", ps%m / ps%Vol
    write (io%uout, '(a,f16.8)') "Charge: ", ps%getCharge()
  end subroutine statsOnField

  subroutine arrayToSet(ps)
    type(particlesType), intent(inout) :: ps

    character(len=elLen), allocatable :: stmp(:)
    integer                           :: i, ns, p
    integer, allocatable              :: istmp(:)

    allocate (stmp(ps%nGParticles))
    allocate (istmp(ps%nGParticles))
    ns = 0
    do i = 1, ps%nGParticles
      if (isInSet(ps%labels(i), stmp, ns, p)) then
        istmp(p) = istmp(p) + 1
      else
        ns = ns + 1
        stmp(ns) = ps%labels(i)
        istmp(ns) = 1
      end if
    end do

    ps%nSpecies = ns
    allocate (ps%species(ps%nSpecies))
    allocate (ps%isp(ps%nSpecies))
    allocate (ps%mass(ps%nSpecies))
    allocate (ps%charge(ps%nSpecies))
    allocate (ps%colour(ps%nSpecies))
    allocate (ps%pots(ps%nSpecies * (ps%nSpecies + 1) / 2))
    ps%species = stmp(1:ns)
    ps%isp = istmp(1:ns)
    deallocate (stmp, istmp)
  end subroutine arrayToSet

  subroutine energyKinetic(ps)
    class(particlesType), intent(inout) :: ps

    integer :: i

    ps%ke = 0.0_rp
    do i = 1, ps%nGParticles
      ps%ke = ps%ke + ps%mass(ps%spec(i)) * (ps%vx(i)**2 + ps%vy(i)**2 + ps%vz(i)**2)
    end do
    ps%temperature = ps%ke / (3.0_rp * kB * ps%nGParticles)
    ps%ke = 0.5_rp * ps%ke

  end subroutine energyKinetic

  subroutine energyPair(ps)
    type(particlesType), intent(inout) :: ps

    integer :: i

    ps%engPair = 0.0_rp
    do i = 1, ps%nGParticles
      ps%engPair = ps%engPair + ps%engStress(i)%engPair
    enddo
  end subroutine energyPair

  subroutine energyElectrostatics(ps)
    type(particlesType), intent(inout) :: ps

    integer :: i

    ps%engElec = 0.0_rp
    do i = 1, ps%nGParticles
      ps%engElec = ps%engElec + ps%engStress(i)%engElec
    enddo
    ps%engElec = ps%engElec + ps%eeself

  end subroutine energyElectrostatics

  subroutine energy(ps)
    class(particlesType), intent(inout) :: ps

    call energyPair(ps)
    call energyElectrostatics(ps)
    call ps%energyKinetic()
    ps%eng = ps%engPair + ps%engElec + ps%englrc + ps%ke

  end subroutine energy

  subroutine writeTrajectory(ps, io, nt, t, dt, level, isFirst)
    class(particlesType)             :: ps
    type(ioType), intent(in)         :: io
    integer, intent(in)              :: nt
    real(rp), intent(in)             :: t, dt
    integer, intent(in), optional    :: level
    logical, intent(in), optional    :: isFirst

    integer :: i, lev, sp
    logical :: header

    if (present(level)) then
      lev = level
    else
      lev = 2
    endif

    if (present(isFirst)) then
      header = isFirst
    else
      header = .false.
    endif

    if (header) then
      write (io%uTraj, *) trim(ps%systemName)
      write (io%uTraj, *) lev, "3", ps%nGParticles
    endif
    write (io%uTraj, *) "timestep", nt, ps%nGParticles, lev, " 3 ", dt, t
    write (io%uTraj, *) ps%h(1, :)
    write (io%uTraj, *) ps%h(2, :)
    write (io%uTraj, *) ps%h(3, :)
    do i = 1, ps%nGParticles
      sp = ps%spec(i)
      write (io%uTraj, *) ps%labels(i), i, ps%mass(sp), "0.0 ", "0.0"
      write (io%uTraj, *) ps%x(i), ps%y(i), ps%z(i)
      if (lev > 0) write (io%uTraj, *) ps%vx(i), ps%vy(i), ps%vz(i)
      if (lev > 1) write (io%uTraj, *) ps%fx(i), ps%fy(i), ps%fz(i)
    enddo
  end subroutine writeTrajectory

  subroutine scale_velocities(ps, T)
    class(particlesType) :: ps
    real(rp)             :: T

    integer  :: i, j
    real(rp) :: et(3)

    et = 0.0_rp
    do i = 1, ps%nGParticles

      j = ps%spec(i)
      et(1) = et(1) + ps%mass(j) * ps%vx(i)**2
      et(2) = et(2) + ps%mass(j) * ps%vy(i)**2
      et(3) = et(3) + ps%mass(j) * ps%vz(i)**2
    end do
    et = et / (kB * ps%nGParticles)
    et = sqrt(T / et)
    do i = 1, ps%nGParticles
      ps%vx(i) = ps%vx(i) * et(1)
      ps%vy(i) = ps%vy(i) * et(2)
      ps%vz(i) = ps%vz(i) * et(3)
    end do

  end subroutine scale_velocities

  subroutine init_velocities(ps, T)
    class(particlesType) :: ps
    real(rp)             :: T

    integer               :: i, j
    real(rp)              :: f, mt, vc(3)
    real(rp), allocatable :: x(:)

    allocate (x(3 * ps%nGParticles))

    vc = 0.0_rp
    mt = 0.0_rp
    call gaussian(x)
    do i = 1, ps%nGParticles
      j = ps%spec(i)
      f = sqrt(kb * T / ps%mass(j))
      ps%vx(i) = f * x(3 * (i - 1) + 1)
      ps%vy(i) = f * x(3 * (i - 1) + 2)
      ps%vz(i) = f * x(3 * (i - 1) + 3)
      vc = vc + ps%mass(j)*[ps%vx(i), ps%vy(i), ps%vz(i)]
      mt = mt + ps%mass(j)
    end do
    vc = vc / mt / ps%nGParticles
    do i = 1, ps%nGParticles
      ps%vx(i) = ps%vx(i) - vc(1)
      ps%vy(i) = ps%vy(i) - vc(2)
      ps%vz(i) = ps%vz(i) - vc(3)
    end do
    deallocate (x)
    call ps%scale_velocities(T)
  end subroutine init_velocities

end module m_particles
