!    Copyright (c) 2016-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT
program microMD
  ! the name shall be μMD not microMD
  use m_Constants,  only: engUnits,&
                          k_ml,&
                          rp
  use m_Useful,     only: DateAndTime,&
                          compilerInfo,&
                          init_random
  use m_config,     only: readConfig
  use m_control,    only: controlType,&
                          readControl
  use m_field,      only: readField
  use m_forces,     only: computeForces,long_range_correction
  use m_io,         only: closeIO,&
                          ioType
  use m_neighbours, only: buildNeighbours,&
                          setupNeighbours
  use m_particles,  only: particlesType
  use m_sampler,    only: nve_update_position,&
                          nve_update_velocity
  use timer, only : timer_type,start_timer, stop_timer,init_timer_system,timer_report
  implicit none

  character(len=k_ml) :: dummy
  character(len=10) :: dt
  character(len=12) :: tm

  type(ioType)        :: io
  type(controlType)   :: control
  type(particlesType) :: particles
  type(timer_type) :: tmr

  call DateAndTime(dt, tm)
  write (dummy, '(a,a,a,a)') "Program µMD has started at ", dt, " ", tm
  if (command_argument_count() == 1) then
    call get_command_argument(1, io%controlFile)
  end if
  call readControl(io, control)
  call init_timer_system(tmr, io%uout,3)
  call compilerInfo(io%uout)
  write (io%uout, '(a)') trim(dummy)
  call start_timer(tmr, 'setup')
  call init_random(control%seed)
  call readConfig(particles, io)
  call particles%summary(io)
  call readField(io, particles)
  call particles%fieldSummary(io)
  call setupNeighbours(particles)
  call stop_timer(tmr, 'setup')

  write (io%uout, '(a,a)') "energy units: ", trim(particles%units)
  call long_range_correction(particles,control%rc)
  write (io%uout, '(a,es16.8)') "energy lrc: ",particles%englrc/engUnits
  block
    integer :: t
    character(len=k_ml) ::  fmte,fmts
    call start_timer(tmr, 'neighbours')
    call buildNeighbours(particles, control%rc)
    call stop_timer(tmr, 'neighbours')

    call start_timer(tmr, 'short range')
    call computeForces(particles, control)
    call stop_timer(tmr, 'short range')

    call particles%init_velocities(control%temperature)
    call particles%energy()
    call start_timer(tmr, 'ouput')
    if (io%isTraj) then
      call particles%writeTrajectory(io, 0, 0.0_rp, 0.0_rp, isFirst=.true., level=2)
    endif
    fmte = '('//"i8,1x,5(es13.6,1x)"//')'
    fmts = '('//"a8,1x,5(a13,1x)"//')'
    write (io%utimeser, fmt=trim(fmts)) "Timestep","Time", "TotEng", "T", "Epair","E_tail"
    write (io%utimeser, fmt=trim(fmte)) 0,control%time, particles%eng / engUnits, particles%temperature, &
      particles%engPair / engUnits, particles%englrc / engUnits
    call stop_timer(tmr, 'ouput')
    do t = 1, control%steps
      control%time = control%time + control%timestep
      control%step = control%step + 1
    call start_timer(tmr, 'neighbours')
    call buildNeighbours(particles, control%rc)
    call stop_timer(tmr, 'neighbours')
      call nve_update_velocity(particles, control%timestep)
      call nve_update_position(particles, control%timestep)
      call start_timer(tmr, 'short range')
      call computeForces(particles, control)
      call stop_timer(tmr, 'short range')
      call nve_update_velocity(particles, control%timestep)
      call start_timer(tmr, 'ouput')
      if ( (t==1) .or.(mod(t,control%freq) == 0) ) then
        call particles%energy()
        write (io%utimeser, fmt=trim(fmte)) t,control%time, particles%eng / engUnits, particles%temperature, &
          particles%engPair / engUnits, particles%englrc / engUnits
      end if
      if (io%isTraj) then
        call particles%writeTrajectory(io, control%step, control%time, control%timestep, isFirst=.false., level=2)
      endif
      call stop_timer(tmr, 'ouput')
    enddo
    call timer_report(tmr)
  end block
  call closeIO(io)
  call DateAndTime(dt, tm)
  write (io%uout, '(a,a,a,a)') "Program µMD has ended at ", dt, " ", tm

end program microMD
