!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2020  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Simple command line argument parser
module dftbp_argparser
  use dftbp_charmanip, only : toLower
  implicit none
  private

  public :: TArgParser, init


  !> Command line argument
  type :: TArgument

    !> Has the argument been processed
    logical :: unused

    !> Corresponds to an existing file
    logical :: isFile

    !> Looks like a command line flag
    logical :: isFlag

    !> Raw value of the argument
    character(len=:), allocatable :: raw

  end type TArgument


  !> Command line argument parser
  type :: TArgParser
    private

    !> All command line arguments
    type(TArgument), allocatable :: arg(:)

    !> Start of command line flags
    integer :: flagStart

    !> End of command line flags
    integer :: flagEnd

  contains

    !> Basic search of argument list
    procedure, private :: findArg

    !> Basic search of argument list, also markes argument as read
    procedure, private :: getArg

    !> Find a given flag
    procedure :: getFlag

    !> Get file name
    procedure :: countFiles

    !> Get file name
    procedure :: getFile

  end type TArgParser


  !> Initialize command line argument parser
  interface init
    module procedure :: initArgParser
    module procedure :: initArgument
  end interface init


contains


  !> Initialized command line argument
  subroutine initArgument(self, iArg)

    !> Command argument
    type(TArgument), intent(out) :: self

    !> Position of argument
    integer, intent(in) :: iArg

    integer :: length

    self%unused = .true.

    call get_command_argument(iArg, length=length)
    allocate(character(len=length) :: self%raw)
    call get_command_argument(iArg, self%raw)

    inquire(file=self%raw, exist=self%isFile)
    self%isFlag = index(self%raw, '-') == 1

  end subroutine initArgument


  !> Initialize command line argument parser
  subroutine initArgParser(self)

    !> Instance of the parser
    type(TArgParser), intent(out) :: self

    integer :: nArg, iArg

    nArg = command_argument_count()
    allocate(self%arg(nArg))

    do iArg = 1, nArg
      call init(self%arg(iArg), iArg)
    end do

    self%flagStart = 1
    iArg = self%getArg('--')
    if (iArg > 0) then
      self%flagEnd = iArg
    else
      self%flagEnd = nArg
    end if

  end subroutine initArgParser


  !> Basic search of argument list
  pure function findArg(self, arg, iStart, iEnd) result(position)

    !> Instance of the parser
    class(TArgParser), intent(in) :: self

    !> Argument
    character(len=*), intent(in) :: arg

    !> Start position
    integer, intent(in), optional :: iStart

    !> End position
    integer, intent(in), optional :: iEnd

    !> Position of argument
    integer :: position

    integer :: iArg, jStart, jEnd

    if (present(iStart)) then
      jStart = max(1, iStart)
    else
      jStart = 1
    end if

    if (present(iEnd)) then
      jEnd = min(iEnd, size(self%arg))
    else
      jEnd = size(self%arg)
    end if

    position = 0
    do iArg = jStart, jEnd
      if (self%arg(iArg)%unused) then
        if (toLower(arg) == toLower(self%arg(iArg)%raw)) then
          position = iArg
          exit
        end if
      end if
    end do

  end function findArg


  !> Basic search of argument list, also markes argument as read
  function getArg(self, arg, iStart, iEnd) result(position)

    !> Instance of the parser
    class(TArgParser), intent(inout) :: self

    !> Argument
    character(len=*), intent(in) :: arg

    !> Start position
    integer, intent(in), optional :: iStart

    !> End position
    integer, intent(in), optional :: iEnd

    !> Position of argument
    integer :: position

    position = self%findArg(arg, iStart, iEnd)
    if (position > 0) then
      self%arg(position)%unused = .false.
    end if

  end function getArg


  !> Search for a command line flag
  function getFlag(self, flag) result(hasFlag)

    !> Instance of the parser
    class(TArgParser), intent(inout) :: self

    !> Flag
    character(len=*), intent(in) :: flag

    !> Flag is present
    logical :: hasFlag

    hasFlag = self%getArg(flag, self%flagStart, self%flagEnd) > 0

  end function getFlag


  !> Count number of unprocessed files
  function countFiles(self) result(nFiles)

    !> Instance of the parser
    class(TArgParser), intent(in) :: self

    !> Number of unprocessed files
    integer :: nFiles

    nFiles = count(self%arg%unused .and. self%arg%isFile)

  end function countFiles


  !> Return file name
  subroutine getFile(self, file)

    !> Instance of the parser
    class(TArgParser), intent(inout) :: self

    character(len=:), allocatable, intent(out) :: file

    integer :: iArg

    do iArg = 1, size(self%arg)
      if (self%arg(iArg)%unused .and. self%arg(iArg)%isFile) then
        file = self%arg(iArg)%raw
        self%arg(iArg)%unused = .false.
        exit
      end if
    end do

  end subroutine getFile


end module dftbp_argparser
