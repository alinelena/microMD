!    Copyright (c) 2006-2016 Alin Marin Elena <alinm.elena@gmail.com>
!    The MIT License http://opensource.org/licenses/MIT

!> \brief deals with all subroutines related with the parsing
!> \author Alin M. Elena
!> \date 14th of January 2006
!> \remark 17th of July 2007 removed the save attribute from variables
!> \remark July 2016 extended and adapted
!> \internal suspected memory leak
module m_Parser
  use iso_fortran_env, only: ERROR_UNIT
  use m_Constants,     only: k_ml,&
                             k_mw,&
                             rp
  use m_Useful,        only: cstr,&
                             error

  implicit none
  private
  ! equalize lengths by space padding at the end
  character(len=*), parameter :: allowedDups(3) = ['io   ', 'no   ', 'print']
  public :: ParseFile, EndParse
  public :: GetLogical, GetString, GetReal, GetInteger, GetLine, GetBlock
  public :: PrintReport
!
!> data structure of a node in the parsing tree
  type, private :: names
    character(len=k_mw)  :: nam = '' !< name of the entity read in file
    character(len=k_ml)  :: value = '' !< the value associated with the node
    logical               :: processed = .false.
    integer :: ut = 0 !< unit number where write the info to be parsed (write(ut,*))
    integer :: lines = 0
    type(names), pointer :: next => null() !< next element in tree
  end type names
!
!> data structure of info about the parsed files
  type, private :: infoParseType
    integer :: comments = 0 !< number of comments parsed
    integer :: empty = 0 !< number of empty lines
    integer :: tokens = 0 !< number of tokens
    integer :: lines = 0 !< number of lines
    integer :: blocks = 0 !< number of blocks
    integer :: blocklines = 0 !< number of lines in blocks
    character(len=256) :: filename
  end type infoParseType
  type(infoParseType) :: ireport
!
  type(names), pointer :: tnames, currentt !< blocks and tokens lists for the files plus the current ones
  type(names), pointer :: bnames, currentb !< blocks and tokens lists for the files plus the current ones
!
contains
!>  \brief opens the input and error files
!>  \details allocates the lists for blocks and tokens
!>  starts the parsing
!> \author Alin M. Elena (Belfast)
!> \date 14th-15th of January 2006
!> \param ioLoc is type(ioType) (see m_Types::ioType)
!> \remark 20th of January 2007 added the parameter \em ioLoc
  subroutine ParseFile(uinp)
    integer, intent(in)    :: uinp

    ireport%empty = 0
    ireport%comments = 0
    ireport%tokens = 0
    ireport%lines = 0
    ireport%blocks = 0
    ireport%blocklines = 0
!
!allocate the lists for name of the blocks and tokens
    allocate (tnames)
    allocate (bnames)

    call parse(uinp)
    inquire (unit=uinp, name=ireport%filename)
  end subroutine ParseFile
!
  subroutine PrintReport(iu)
    integer, intent(in)    :: iu

    write (iu, '(a42)') repeat('=', 42)
    write (iu, '(a)') "Parsing report for file: "//trim(ireport%filename)
    write (iu, '(a)') "1. Token labels"
    call PrintName(tnames, iu)
    write (iu, '(a)') "2. Block labels"
    call PrintName(bnames, iu)
    write (iu, '(a)') "3. General info"
    write (iu, '(a,i5)') "No of lines read", ireport%lines
    write (iu, '(a,i5)') "No of comment lines", ireport%comments
    write (iu, '(a,i5)') "No of empty lines", ireport%empty
    write (iu, '(a,i5)') "No of tokens", ireport%tokens
    write (iu, '(a,i5)', advance='no') "No of blocks", ireport%blocks
    write (iu, '(a,i5,a)') " containing ", ireport%blocklines, " lines"
    write (iu, '(a42)') repeat('=', 42)
  end subroutine PrintReport
!
!> \brief    recursive subroutine which effectively parses the file
!> \details    associated with unit=nounit
!>          if an include statement is found the routine is called again to parse the new file
!> \author Alin M. Elena (Belfast)
!> \date 14th-15th of January 2006
!> \param nounit integer, represents the unit number of a file which was previous opened and now
!> is parsed
!> \param  ioLoc type(ioType) (see ioType)
!> \param  ireport type(infoParseType),  keeps the info about the parsed files (see m_Types::infoParseType type)
!> \remark 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added parameter \em ioLoc
!
  subroutine parse(nounit)
    integer, intent(in)    :: nounit

    character(len=*), parameter :: myName = "parse"

    character(len=k_ml) :: blockfile, line, lineaux, nam, storeline
    character(len=k_mw) :: blockname, endblockname, tokename
    integer             :: blines, errno, lineno, nt
    logical             :: bool

!

!
    lineno = 0
    inquire (unit=nounit, name=nam)
    do while (GetLine(nounit, line))
      lineaux = adjustl(line)
      lineno = lineno + 1
      if (lineaux(1:1) == "!" .or. lineaux(1:1) == "#" .or. lineaux(1:2) == "//") then
        ireport%comments = ireport%comments + 1
! just count the comments nothing else
      else if (len(trim(line)) == 0) then
        ireport%empty = ireport%empty + 1
      else if (cstr(lineaux(1:5), "block")) then
        ! defines the behaviour of a block - endblock statement
        ireport%blocks = ireport%blocks + 1
        read (lineaux(6:k_ml), *, iostat=errno) blockname
        if (errno /= 0) then
          call ParseErr("invalid name after block statement", myname, nam, lineno)
        endif
        bool = .true.
        blines = 0
        blockfile = trim(blockname)//".blk"
        !  create a scratch file
        open (newunit=nt, file=trim(blockfile), status="replace", action="write")
        do while (bool)
          bool = GetLine(nounit, line)
          lineno = lineno + 1
          lineaux = adjustl(line)
          ! take some action if we have a endblock
          if (cstr(lineaux(1:8), "endblock")) then
            ! check the name
            read (lineaux(9:k_ml), *, iostat=errno) endblockname
            if (errno /= 0) then
              call ParseErr("invalid name after endblock statement", myname, nam, lineno)
            endif
            ! closing the wrong block
            if (.not. cstr(trim(endblockname), trim(blockname))) then
              call ParseErr("closing wrong block, expected endblock "//&
               & trim(blockname)//" found endblock "//trim(endblockname), myname, nam, lineno)
            endif
            exit
          end if
          ! parse the content of the block line by line and get rid of commented and empty lines
          if (lineaux(1:1) == "!" .or. lineaux(1:1) == "#" .or. lineaux(1:2) == "//") then
            ireport%comments = ireport%comments + 1
          else if (len(trim(line)) == 0) then
            ireport%empty = ireport%empty + 1
          else
            ! count only not obviuos without value lines
            ! get rid of any inline comment
            storeline = ParseLine(lineaux)
            write (nt, '(a)') trim(storeline)
            blines = blines + 1
          end if
        end do
        ! check if we have reached end-of-file
        if (.not. bool) call ParseErr("missing endblock "//trim(blockname), myname, nam, -1)
        ! check empty blocks and give an warning, empty blocks are ignored
        if (blines == 0) then
          call ParseWar("detected empty(probably only comments and empty lines) block "//trim(blockname), myname, nam, lineno)
          ! empty files we delete them immediately
        end if
        ! check the uniquness of the block
        ireport%blocklines = blines + ireport%blocklines
        !
        close (nt)
        ! check for the block in the existent list
        if (ireport%blocks == 1) then
          currentb => bnames
          bnames%nam = trim(blockname)
          bnames%lines = blines
          bnames%value = trim(blockfile)
        else
          if (FindName(trim(blockname), bnames)) call ParseErr("found block "//trim(blockname)//" duplicated", myname, nam, &
            & lineno)
          call AddName(trim(blockname), currentb, blines)
          currentb%value = trim(blockfile)
        end if
      else if (cstr(lineaux(1:8), "endblock")) then
        ! endblock without block
        call ParseErr("endblock without block", myname, nam, lineno)
      else
        ireport%tokens = ireport%tokens + 1
        call getToken(lineaux, tokename, storeline)
! check for the token in the existent list and add it if is new
        if (ireport%tokens == 1) then
          currentt => tnames
          tnames%nam = "title"
          if (cstr(trim(tokename), "title")) then
            tnames%value = trim(storeline)
          else
            tnames%value = trim(lineaux)
          end if
        else
          if (FindName(trim(tokename), tnames)) then
            call ParseErr("found token "//trim(tokename)//" duplicated", myname, nam, lineno)
          end if
          call AddName(trim(tokename), currentt)
          currentt%value = trim(storeline)
        end if
        if (cstr(trim(tokename), "finish")) exit
      end if
    end do
    ireport%lines = ireport%lines + lineno
  end subroutine parse

  subroutine getToken(lineaux, tokename, storeline)
    character(len=*), intent(in)    :: lineaux
    character(len=*), intent(out)   :: tokename, storeline

    character(len=k_ml) :: lineaux2
    character(len=k_mw) :: t2
    integer             :: errno

    read (lineaux, *, iostat=errno) tokename
    lineaux2 = ParseLine(lineaux)
    if (isInList(tokename, allowedDups)) then
      read (lineaux2(len(trim(tokename)) + 1:k_ml), *, iostat=errno) t2
      tokename = trim(tokename)//" "//trim(t2)
    end if
    read (lineaux2(len(trim(tokename)) + 1:k_ml), *, iostat=errno) storeline
    if (errno /= 0) storeline = ''
  end subroutine getToken

  logical function isInList(token, list)
    character(len=*), intent(in)    :: token, list(:)

    integer :: i

    isInList = .false.
    do i = 1, size(list)
      if (cstr(trim(token), trim(list(i)))) then
        isInList = .true.
        exit
      end if
    end do
  end function isInList

!
!> \brief logical function reads a line from a unit
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param uno integer, unit number from where to read the line
!> \param line character(len=*), the line that was read
!> \return .true. if successfull in reading the line, .false. otherwise
!> \remarks
!
  logical function GetLine(uno, line)
    integer, intent(in)                :: uno
    character(len=k_ml), intent(out)   :: line

    integer :: errno

!
    inquire (unit=uno, iostat=errno)
    GetLine = .false.
    if (errno /= 0) then
      write (*, *) "Unexpected error opening the input file(s)"
      stop
    end if
    read (uno, fmt='(a)', iostat=errno) line
    if (errno == 0) GetLine = .true.
  end function GetLine
!
!> \brief adds a new node in the list
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param word character(len=*), value for field nam of the node
!> \param current pointer to the last node in the list
!> \param lines integer, optional value for field lines
!
  subroutine AddName(word, current, lines)
    character(len=*), intent(in)     :: word
    type(names), pointer             :: current
    integer, intent(in), optional    :: lines

    type(names), pointer :: node

!
    allocate (node)
    node%nam = trim(word)
    current%next => node
    if (present(lines)) node%lines = lines
    current => node
!
  end subroutine AddName
!
!> \brief logical function finds a field in a list starting at root
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param word character(len=*), value of the field nam to be searched for
!> \param root type(names), pointer, starting point in search
!> \param loc tppe(names), pointer, optional returns the location where the info was found.
  logical function FindName(word, root, loc)
    character(len=*), intent(in)    :: word
    type(names), pointer            :: root
    type(names), optional, pointer  :: loc

    type(names), pointer :: current

!
    current => root
    FindName = .false.
    do while (associated(current))
      if (cstr(trim(word), trim(current%nam))) then
        FindName = .true.
        if (present(loc)) loc => current
        exit
      end if
      current => current%next
    end do
  end function FindName
!
!> \brief Prints at a specified unit the nam field  from a list
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param root type(names), pointer  starting node in the Printing list
!> \param ioLoc type(ioType) conatins the unit to Print (see m_Types::ioType)
!> \remarks
  subroutine PrintName(root, ounit)
    type(names), pointer   :: root
    integer, intent(in)    :: ounit

    type(names), pointer :: current

!
    current => root
    do while (associated(current))
      write (ounit, '(a)') trim(current%nam)
      current => current%next
    end do
  end subroutine PrintName
!
!> \brief recursive function that deallocates all the nodes starting with root
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param root type(names), pointer the starting node
!> \remarks
  recursive subroutine DeleteList(root)
    type(names), pointer :: root

    type(names), pointer :: current

!
    current => root%next
    if (associated(current)) then
      call DeleteList(current)
    else
      if (root%ut /= 0) close (root%ut)
      deallocate (root)
    end if
  end subroutine DeleteList
!
!> \brief parses a line removing comments
!> \details everything after ! \# or // is considered a comment
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param line character(len=*) the line to be parsed
  character(k_ml) function ParseLine(line)
    character(len=*), intent(in)    :: line

    integer :: p(1:3), pos

!
    p(1) = scan(line, "!")
    p(2) = scan(line, "#")
    p(3) = index(line, "//")
    where (p == 0) p = k_ml + 1
    pos = minval(p) - 1
    if (pos /= k_ml) then
      ParseLine = line(1:pos)
    else
      ParseLine = line
    end if
  end function ParseLine
!> \brief Prints an error message and aborts the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param message the error message
!> \param routine the subprogram that generated the error
!> \param filename raeding this \em filename the error occured
!> \param lineno the line number that generated the error
!> \param ioLoc type(ioType) (see m_Types::ioType)
  subroutine ParseErr(message, routine, filename, lineno)
    character(len=*), intent(in)    :: message, routine, filename
    integer, intent(in)             :: lineno

!
    write (ERROR_UNIT, '(a,a)') "Error: ", trim(message)
    write (ERROR_UNIT, '(a,a)') "Routine: ", trim(routine)
    write (ERROR_UNIT, '(a,a)') "File: ", trim(filename)
    write (ERROR_UNIT, '(a,i7)') "Line number: ", lineno
    write (ERROR_UNIT, '(a)') "User stop"
    write (ERROR_UNIT, '(a)') "User stop, but probably this is not what you want"
    stop
  end subroutine ParseErr
!
!> \brief Prints a warning message occured in the parsing process
!> \details
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param message the error message
!> \param routine the subprogram that generated the warning
!> \param filename raeding this \em filename the warning occured
!> \param lineno the line number that generated the warning
!> \param ioLoc type(ioType) (see m_Types::ioType)
!
  subroutine ParseWar(message, routine, filename, lineno)
    character(len=*), intent(in)    :: message, routine, filename
    integer, intent(in)             :: lineno

!
    write (ERROR_UNIT, '(a,a)') "Warning: ", trim(message)
    write (ERROR_UNIT, '(a,a)') "Routine: ", trim(routine)
    write (ERROR_UNIT, '(a,a)') "File: ", trim(filename)
    write (ERROR_UNIT, '(a,i7)') "Line number: ", lineno
    write (ERROR_UNIT, '(a)') "User warning"
  end subroutine ParseWar
!
!> \brief  returns a logical value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to .true. if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  logical function GetLogical(io, label, dflt)
    integer, intent(in)              :: io
    character(len=*), intent(in)     :: label
    logical, intent(in), optional    :: dflt

    character(len=10)    :: def
    type(names), pointer :: found

!default value is set to .true. if no dflt parameter is present
!
    def = ''
    if (FindName(trim(label), tnames, found)) then
      if ((cstr(trim(found%value), "yes")) .or. (cstr(trim(found%value), "true")) .or. (cstr(trim(found%value), ".true.")) .or. &
     & (cstr(trim(found%value), "t")) .or. (cstr(trim(found%value), "y")) .or. (cstr(trim(found%value), "1"))) then
        GetLogical = .true.
      else if ((cstr(trim(found%value), "no")) .or. (cstr(trim(found%value), "false")) .or. (cstr(trim(found%value), ".false.")) &
     & .or. (cstr(trim(found%value), "f")) .or. (cstr(trim(found%value), "n")) .or. (cstr(trim(found%value), "0"))) then
        GetLogical = .false.
      elseif (len_trim(found%value) == 0) then
        GetLogical = .true.
      else
        call error(__FILE__, __LINE__, -105, &
                   "Unsuported value: "//trim(label)//" "//trim(found%value))
      end if
      found%processed = .true.
    else
      if (present(dflt)) then
        GetLogical = dflt
      else
        GetLogical = .true.
      end if
      def = "! default"
    end if
    write (io, '(a,2x,l1,a)') trim(label), GetLogical, trim(def)
  end function GetLogical
!
!> \brief  returns a string value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to empty string if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \param Print optional logical determines if the token is Printed or not. if the parameter is missing is assuk_med .false.
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  character(len=k_ml) function GetString(io, label, dflt, print)
    integer, intent(in)                       :: io
    character(len=*), intent(in)              :: label
    character(len=*), intent(in), optional    :: dflt
    logical, intent(in), optional             :: print

    character(len=10)    :: def
    type(names), pointer :: found

! default value if no dflt parameter is present is empty string
! so be sure that you provide a default parameter
    def = ''
    if (FindName(trim(label), tnames, found)) then
      GetString = trim(found%value)
      found%processed = .true.
    else
      if (present(dflt)) then
        GetString = trim(dflt)
      else
        GetString = ''
      end if
      def = '! default'
    end if
    if (present(print)) then
      if (print) then
        write (io, '(a,2x,a,1x,a)') trim(label), trim(GetString), trim(def)
      end if
    end if
  end function GetString
!
!> \brief  returns a real(kind=pr) value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to 0.0_pr if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  real(rp) function GetReal(io, label, dflt)
    integer, intent(inout)            :: io
    character(len=*), intent(in)      :: label
    real(rp), intent(in), optional    :: dflt

    character(len=10)    :: def
    character(len=k_mw)  :: aux
    integer              :: errno
    type(names), pointer :: found

! default value if no dflt parameter is present is 0.0_pr
    def = ''
    if (FindName(trim(label), tnames, found)) then
      aux = trim(found%value)
      read (aux, *, iostat=errno) GetReal
      if (errno /= 0) call error(__FILE__, __LINE__, -104, &
                                 "wrong value "//trim(found%value)//" suplied for label "//trim(label))
      found%processed = .true.
    else
      if (present(dflt)) then
        GetReal = dflt
      else
        GetReal = 0.0_rp
      end if
      def = '! default'
    end if
    write (io, '(a,2x,g32.16,1x,a)') trim(label), GetReal, trim(def)
  end function GetReal
!
!> \brief  returns an integer value found in the input file(s) associated with the token \em label
!> \details  if no token is found the value is set to dflt
!> default value is set to 0 if no dflt parameter is present
!> output is put in the units indicated by \em ioLoc
!> \author Alin M Elena
!> \date 14th of January 2006
!> \param label character(len=*), the token to search for
!> \param dflt optional default value for the token \em label
!> \param ioLoc type(ioType) (see m_Types::ioType)
!> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  integer function GetInteger(io, label, dflt)
    integer, intent(inout)           :: io
    character(len=*), intent(in)     :: label
    integer, intent(in), optional    :: dflt

    character(len=10)    :: def
    character(len=k_mw)  :: aux
    integer              :: errno
    type(names), pointer :: found

! default value if no dflt parameter is present is 0
!
    def = ''
    if (FindName(trim(label), tnames, found)) then
      aux = trim(found%value)
      read (aux, *, iostat=errno) GetInteger
      if (errno /= 0) call error(__FILE__, __LINE__, -103, &
                                 "wrong value "//trim(found%value)//" suplied for label "//trim(label))
      found%processed = .true.
    else
      if (present(dflt)) then
        GetInteger = dflt
      else
        GetInteger = 0
      end if
      def = '! default'
    end if
    write (io, '(a,2x,i0,a)') trim(label), GetInteger, trim(def)
  end function GetInteger

  !> \brief  returns a logical value to indicate if it found in the input file(s) a valid block associated with the token \em label
  !> \details  if no token is found the value returned is .false.
  !> a unit from where the content of the block cand be read is returned via nt output is put in the units indicated by ioLoc
  !> \author Alin M Elena
  !> \date 14th of January 2006
  !> \param label character(len=*), the token to search for
  !> \param nt unit number from where to read the block data
  !> \param ioLoc type(ioType) (see m_Types::ioType)
  !> \remarks 20th of January, 2007, by Alin M Elena (Queen's University Belfast), added \em ioLoc parameter
  logical function GetBlock(io, label, nt)
    integer, intent(in)             :: io
    character(len=*), intent(in)    :: label
    integer, intent(out)            :: nt

    character(len=k_mw)  :: line
    type(names), pointer :: found

    !
    if (FindName(trim(label), bnames, found)) then
      GetBlock = .true.
      line = trim(found%value)
      open (newunit=nt, action="read", status="old", file=trim(line))
      found%ut = nt
      found%processed = .true.
      write (io, '(a,2x,l1)') trim(label), GetBlock
    else
      write (io, '(a)') "Block "//trim(label)//" was not found"
      GetBlock = .false.
    end if
  end function GetBlock
  !
  !> \brief ends the parsing process and cuts the \em tree (deallocates the memory used during the parsing process)

  integer function countRead(root)

    type(names), pointer :: root

    integer              :: c
    type(names), pointer :: current

!
    current => root
    c = 0
    do while (associated(current))
      if (current%processed) c = c + 1
      current => current%next
    end do
    countRead = c
  end function countRead

  subroutine printUnprocessed(root)
    type(names), pointer :: root

    type(names), pointer :: current

!
    write (ERROR_UNIT, '(a)') "Unprocessed tokens: "
    current => root
    do while (associated(current))
      if (.not. current%processed) then
        write (ERROR_UNIT, '(a)') trim(current%nam)
      end if
      current => current%next
    end do
  end subroutine printUnprocessed

!> \author Alin M Elena
!> \date 14th of January 2006
!> \warning after the call to this function none of the get_* function will work
!> all the input has to be done before the call to it.
  subroutine EndParse(io)
    integer, intent(inout) :: io

! to be called only when there is nothing to be read

    if (countRead(tnames) /= ireport%tokens) then
      write (ERROR_UNIT, '(a)') "***Unprocessed tokens found: "
      call printUnprocessed(tnames)
      call error(__FILE__, __LINE__, -102, &
                 "Unprocessed tokens exist.!!!")
    else
      write (io, '(a)') "***All found tokens processed"
    end if
    call DeleteList(tnames)
    ! check for unprocessed blocks and go mad if any found
    if (countRead(bnames) /= ireport%blocks) then
      write (ERROR_UNIT, '(a)') "***Unprocessed blocks found: "
      call printUnprocessed(bnames)
      call error(__FILE__, __LINE__, -110, &
                 "Unprocessed blocks exist.!!!")
    else
      write (io, '(a)') "***All found blocks processed"
    end if
    call DeleteList(bnames)
  end subroutine EndParse
!
end module m_Parser
