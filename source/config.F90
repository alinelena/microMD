!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
module m_config
  use m_constants, only: k_ml,&
                         rp
  use m_io,        only: ioType
  use m_particles, only: particlesType
  use m_useful,    only: error,&
                         invertH

  implicit none
  private
  character(len=*), parameter :: mName = "io"
  public :: readConfig

contains

  subroutine readConfig(ps, io)
    type(particlesType), intent(inout) :: ps
    type(ioType), intent(in)           :: io

    character(len=*), parameter :: sName = "readConfig"

    character(len=k_ml) :: serr, title
    integer             :: fmps, i, ierr, imcon, nAtoms, nb, nl, ns, u

    nAtoms = -101
    ns = 0
    open (newunit=u, file=trim(io%configFile), action="read", iostat=ierr)
    if (ierr /= 0) then
      call error(__FILE__, __LINE__, -102, &
                 trim(mName)//"::"//trim(sName)//": error opening "//trim(io%configFile))
    endif
    read (u, '(a)') title
    ps%systemName = trim(title)
    read (u, *, iostat=ierr) fmps, imcon, nAtoms
    if (ierr /= 0 .and. nAtoms == -101) then
      nl = 1
      do
        read (u, *, iostat=ierr)
        if (is_iostat_end(ierr)) exit
        nl = nl + 1
      end do
      if (imcon /= 0) ns = 3
      select case (fmps)
      case (2)
        nb = 4
      case (1)
        nb = 3
      case default
        nb = 2
      end select
      nAtoms = (nl - ns) / nb
    else if (ierr /= 0) then
      call error(__FILE__, __LINE__, -101, &
                 trim(mName)//"::"//trim(sName)//": error reading line 2 in config file")
    endif
    rewind u
    read (u, *); read (u, *)
    if (imcon > 0) then
      read (u, *) ps%h(1, :)
      read (u, *) ps%h(2, :)
      read (u, *) ps%h(3, :)
    endif
    ps%nGParticles = nAtoms
    ps%pbc = imcon
    ps%hi = invertH(ps%h)
    allocate (ps%x(ps%nGParticles))
    allocate (ps%y(ps%nGParticles))
    allocate (ps%z(ps%nGParticles))
    allocate (ps%vx(ps%nGParticles))
    allocate (ps%vy(ps%nGParticles))
    allocate (ps%vz(ps%nGParticles))
    allocate (ps%fx(ps%nGParticles))
    allocate (ps%fy(ps%nGParticles))
    allocate (ps%fz(ps%nGParticles))
    allocate (ps%labels(ps%nGParticles))
    allocate (ps%GIDs(ps%nGParticles))
    allocate (ps%spec(ps%nGParticles))
    allocate (ps%engStress(ps%nGParticles))
    do i = 1, ps%nGParticles
      read (u, *, iostat=ierr) ps%labels(i), ps%GIDs(i)
      if (ierr /= 0) then
        write (serr, '(a,i0)') "error reading line 1 for atom ", i
        call error(__FILE__, __LINE__ - 3, -106, serr)
      endif
      read (u, *, iostat=ierr) ps%x(i), ps%y(i), ps%z(i)
      if (ierr /= 0) then
        write (serr, '(a,i0)') "error reading positions for atom ", i
        call error(__FILE__, __LINE__ - 3, -107, serr)
      endif
      if (fmps > 0) then
        read (u, *) ps%vx(i), ps%vy(i), ps%vz(i)
        if (ierr /= 0) then
          write (serr, '(a,i0)') "error reading velocities for atom ", i
          call error(__FILE__, __LINE__ - 3, -107, serr)
        endif
      else
        ps%vx(i) = 0.0_rp; ps%vy(i) = 0.0_rp; ps%vz(i) = 0.0_rp
      endif
      if (fmps > 1) then
        read (u, *) ps%fx(i), ps%fy(i), ps%fz(i)
        if (ierr /= 0) then
          write (serr, '(a,i0)') "error reading forces for atom ", i
          call error(__FILE__, __LINE__ - 3, -107, serr)
        endif
      else
        ps%fx(i) = 0.0_rp; ps%fy(i) = 0.0_rp; ps%fz(i) = 0.0_rp
      endif
    enddo
    close (u)
  end subroutine readConfig

end module m_config
