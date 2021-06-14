!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Utility containers to store and handle integrals
module dftbp_integrals
  use dftbp_accuracy, only : dp
  implicit none
  private

  public :: TIntegralCont, TIntegralCont_init

  !> Container to store Hamiltonian and overlap related integrals
  type :: TIntegralCont

    !> Non-SCC part of the Hamiltonian matrix in sparse storage
    real(dp), allocatable :: hcore(:)

    !> Sparse Hamiltonian matrix
    real(dp), allocatable :: hamiltonian(:, :)

    !> Imaginary part of the sparse Hamiltonian matrix
    real(dp), allocatable :: iHamiltonian(:, :)

    !> Sparse Overlap matrix
    real(dp), allocatable :: overlap(:)

    !> Sparse dipole moment integrals, moment operator on bra
    real(dp), allocatable :: dipoleBra(:, :)

    !> Sparse dipole moment integrals, moment operator on ket
    real(dp), allocatable :: dipoleKet(:, :)

    !> Sparse quadrupole moment integrals, moment operator on bra
    real(dp), allocatable :: quadrupoleBra(:, :)

    !> Sparse quadrupole moment integrals, moment operator on ket
    real(dp), allocatable :: quadrupoleKet(:, :)
  end type TIntegralCont


contains


  !> Constructor for a new integral container
  subroutine TIntegralCont_init(this)

    !> Instance of the integral container
    type(TIntegralCont), intent(out) :: this
  end subroutine TIntegralCont_init

end module dftbp_integrals
