!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT

!> \brief deals with reading the field file(s)
!> \details it saves the gathered info in the right variables
!> \author Alin M. Elena (Daresbury Laboratory)
!> \date 11th of April 2016
!> \remarks
!
module m_field
  use m_Constants,  only: elLen,&
                          engUnits,&
                          k_ml,&
                          k_mw,&
                          rp
  use m_Parser,     only: EndParse,&
                          GetBlock,&
                          GetLine,&
                          ParseFile,&
                          getString
  use m_Useful,     only: error,&
                          getVdw,&
                          isInSet
  use m_io,         only: ioType
  use m_particles,  only: particlesType
  use m_potentials, only: LJ,&
                          LJAB,&
                          LJS
  use m_units,      only: setEnergyUnits

  implicit none
  private
  public :: readField
  !
contains
  !
  !

  subroutine readPairs(io, ps)
    integer, intent(inout)             :: io
    type(particlesType), intent(inout) :: ps

    integer :: ub

    if (GetBlock(io, "pairs", ub)) then
      block
        integer :: i, j, k
        character(len=k_ml) :: line
        character(len=elLen) :: a, b
        character(len=k_mw) :: pot
        integer :: ierr
        real(rp) :: params(8)

        do while (GetLine(ub, line))
          params = 0.0_rp
          read (line, *, iostat=ierr) a, b, pot, params
          if (.not. isInSet(a, ps%species, ps%nSpecies, i)) then
            call error(__FILE__, __LINE__, -121, &
                       "error unknown specie "//trim(a))
          endif
          if (.not. isInSet(b, ps%species, ps%nSpecies, j)) then
            call error(__FILE__, __LINE__, -122, &
                       "error unknown specie "//trim(b))
          endif
          k = getVdw(i, j, ps%nSpecies)
          select case (trim (pot))
          case ("lj")
            ps%vdw(k)%i = i
            ps%vdw(k)%j = j
            ps%vdw(k)%id = LJ
            ps%vdw(k)%np = 3
            ps%vdw(k)%pot = trim(pot)
            ps%vdw(k)%fpot = 'Lennard-Jones ε,σ,c'
            allocate (ps%vdw(k)%param(3))
            ps%vdw(k)%param(1:3) = params(1:3)
            ps%vdw(k)%param(1) = ps%vdw(k)%param(1) * engUnits
          case ("ljs")
            ps%vdw(k)%i = i
            ps%vdw(k)%j = j
            ps%vdw(k)%id = LJS
            ps%vdw(k)%np = 3
            ps%vdw(k)%pot = trim(pot)
            ps%vdw(k)%fpot = 'Lennard-Jones ε,σ,c smoothed'
            allocate (ps%vdw(k)%param(3))
            ps%vdw(k)%param(1:3) = params(1:3)
            ps%vdw(k)%param(1) = ps%vdw(k)%param(1) * engUnits
          case ("ljAB")
            ps%vdw(k)%i = i
            ps%vdw(k)%j = j
            ps%vdw(k)%id = LJAB
            ps%vdw(k)%np = 3
            ps%vdw(k)%pot = trim(pot)
            ps%vdw(k)%fpot = 'Lennard-Jones ε,σ,c smoothed'
            allocate (ps%vdw(k)%param(3))
            ps%vdw(k)%param(1:3) = params(1:3)
            ps%vdw(k)%param(1) = ps%vdw(k)%param(1) * engUnits
            ps%vdw(k)%param(2) = ps%vdw(k)%param(2) * engUnits

          case default
            call error(__FILE__, __LINE__, -123, &
                       "unknown potential "//trim(pot))
          end select
        enddo
      end block
      close (ub, status="delete")
    else
      call error(__FILE__, __LINE__, -115, "Cannot read block pairs")
    endif

  end subroutine readPairs

  subroutine readSpecies(io, ps)
    integer, intent(inout)             :: io
    type(particlesType), intent(inout) :: ps

    integer :: ub

    if (GetBlock(io, "species", ub)) then
      block
        integer :: i, j
        character(len=elLen) :: el
        real(rp) :: m, c, col
        integer :: ierr
        logical :: specs(ps%nSpecies)

        specs = .true.
        do i = 1, ps%nSpecies
          read (ub, *, iostat=ierr) el, m, c, col
          if (ierr /= 0) then
            call error(__FILE__, __LINE__, -112, &
                       "Error reading data for specie "//trim(ps%species(i)))
          endif
          if (isInSet(el, ps%species, ps%nSpecies, j)) then
            if (specs(j)) then
              specs(j) = .false.
            else
              call error(__FILE__, __LINE__, -114, &
                         "Specie "//trim(el)//" already read")
            endif
            ps%mass(j) = m; ps%charge(j) = c
            ps%colour(j) = col
          else
            call error(__FILE__, __LINE__, -112, &
                       "Specie "//trim(el)//" not found!!!")
          endif
        enddo
      end block
      close (ub, status="delete")
    else
      call error(__FILE__, __LINE__, -111, "Cannot read block species")
    endif

  end subroutine readSpecies

  subroutine readField(io, ps)
    type(ioType), intent(inout)        :: io
    type(particlesType), intent(inout) :: ps

    character(len=k_ml) :: dummy
    integer             :: errno = -1, id

    open (newunit=io%ufield, file=trim(io%fieldFile), status='old', action='read', iostat=errno)
    if (errno /= 0) then
      call error(__FILE__, &
                 __LINE__, &
                 -110, "Cannot open "//trim(io%fieldFile)//" file")
    end if
    ! parse file and in the same time
    ! makes a minimal check of corectness
    call ParseFile(io%ufield)
    id = io%udeb
    dummy = getString(id, "title", "no title", .true.)
    ps%units = trim(getString(id, "units", "kcal/mol", .true.))
    call setEnergyUnits(ps%units)
    call readSpecies(id, ps)
    call readPairs(id, ps)

    dummy = getString(id, "finish", "", .true.)
    call EndParse(id)
  end subroutine readField
end module m_field
