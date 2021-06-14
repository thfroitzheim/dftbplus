!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2021  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!
#:include "common.fypp"
#:include "error.fypp"

!> Proxy module for interfacing with the tblite library.
!>
!> Todo: Conventions for spherical harmonics are different in tblite and DFTB+.
!>       This requires shuffling them correctly when returned from the library.
!>
!> ang.  | DFTB+
!> ----- | ------------------------------
!> 0 (s) | 0
!> 1 (p) | -1, 0, 1
!> 2 (d) | -2, -1, 0, 1, 2
!> 3 (f) | -3, -2, -1, 0, 1, 2, 3
!> 4 (g) | -4, -3, -2, -1, 0, 1, 2, 3, 4
module dftbp_tblite
  use dftbp_accuracy, only : dp
  use dftbp_blasroutines, only : gemv
  use dftbp_charges, only : getSummedCharges
  use dftbp_commontypes, only : TOrbitals
  use dftbp_environment, only : TEnvironment
  use dftbp_message, only : error
  use dftbp_periodic, only : TNeighbourList
  use dftbp_schedule, only : distributeRangeInChunks, assembleChunks
  use dftbp_simplealgebra, only : determinant33
#:if WITH_TBLITE
  use mctc_env, only : error_type
  use mctc_io, only : structure_type, new
  use mctc_io_symbols, only : symbol_length
  use tblite_basis_type, only : get_cutoff, basis_type
  use tblite_context_type, only : context_type
  use tblite_coulomb_cache, only : coulomb_cache
  use tblite_cutoff, only : get_lattice_points
  use tblite_disp_cache, only : dispersion_cache
  use tblite_integral_multipole, only : multipole_cgto, multipole_grad_cgto
  use tblite_integral_overlap, only : overlap_cgto, overlap_grad_cgto, maxl, msao
  use tblite_scf_potential, only : potential_type, new_potential
  use tblite_scf_info, only : scf_info, atom_resolved, shell_resolved, orbital_resolved, &
    & not_used
  use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
  use tblite_xtb_calculator, only : xtb_calculator
  use tblite_xtb_gfn2, only : new_gfn2_calculator
  use tblite_xtb_gfn1, only : new_gfn1_calculator
  use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_occupation, &
      & get_hamiltonian_gradient, tb_hamiltonian
  use tblite_xtb_ipea1, only : new_ipea1_calculator
  use tblite_xtb_singlepoint, only : xtb_singlepoint
  use tblite_version, only : get_tblite_version
#:endif
  implicit none
  private

  public :: TTBLite, TTBLiteInput, TTBLite_init, writeTBLiteInfo
  public :: gfn2xtb, gfn1xtb, ipea1xtb


  !> Type to implement an enumerated method selector
  type :: TMethodSelector
    private
    integer :: id
  end type TMethodSelector

  !> Selector for GFN2-xTB Hamiltonian
  type(TMethodSelector), parameter :: gfn2xtb = TMethodSelector(2)

  !> Selector for GFN1-xTB Hamiltonian
  type(TMethodSelector), parameter :: gfn1xtb = TMethodSelector(1)

  !> Selector for IPEA1-xTB Hamiltonian
  type(TMethodSelector), parameter :: ipea1xtb = TMethodSelector(11)


  !> Input for the library
  type :: TTBLiteInput

    !> Selected method
    type(TMethodSelector) :: method

  end type TTBLiteInput


  !> Library interface handler
  type :: TTBLite
  #:if WITH_TBLITE
    !> Molecular structure data
    type(structure_type) :: mol

    !> Calculation context
    type(context_type) :: ctx

    !> Parametrisation data
    type(xtb_calculator) :: calc

    !> Wavefunction data
    type(wavefunction_type) :: wfn

    !> Density-dependent potential data
    type(potential_type) :: pot

    !> Reuseable data for Coulombic interactions
    type(coulomb_cache) :: cache

    !> Reuseable data for Dispersion interactions
    type(dispersion_cache) :: dcache
  #:endif

    !> Mapping between species and identifiers
    integer, allocatable :: sp2id(:)

    !> Coordination number
    real(dp), allocatable :: cn(:)

    !> Derivative of the coordination number w.r.t. atomic displacements
    real(dp), allocatable :: dcndr(:, :, :)

    !> Derivative of the coordination number w.r.t. strain deformations
    real(dp), allocatable :: dcndL(:, :, :)

    !> Diagonal elements of the Hamiltonian
    real(dp), allocatable :: selfenergy(:)

    !> Derivatives of the diagonal elements w.r.t. the coordination number
    real(dp), allocatable :: dsedcn(:)

    !> Repulsion energy
    real(dp) :: erep

    !> Non-self consistent dispersion energy
    real(dp) :: edisp

    !> Self-consistent dispersion energy
    real(dp) :: escd

    !> Electrostatic energy
    real(dp) :: ees

    !> Contributions to the gradient
    real(dp), allocatable :: gradient(:, :)

    !> Contributions to the virial
    real(dp) :: sigma(3, 3)

  contains

    !> Update internal copy of coordinates
    procedure :: updateCoords

    !> Update internal copy of lattice vectors
    procedure :: updateLatVecs

    !> Get real space cutoff
    procedure :: getRCutoff

    !> Get energy contributions
    procedure :: getEnergies

    !> Get force contributions
    procedure :: addGradients

    !> Get stress tensor contributions
    procedure :: getStress

    !> Updates with changed populations for the instance.
    procedure :: updatePopulations

    !> Returns shifts per atom
    procedure :: getShifts

    !> Get orbital information
    procedure :: getOrbitalInfo

    !> Get reference occupation
    procedure :: getReferenceN0

    !> Returns the equivalence to get the correct mixing of charge dependent contributions
    procedure :: getOrbitalEquiv

    !> Construct Hamiltonian and overlap related integrals
    procedure :: buildSH0

    !> Evaluate shift related derivatives from Hamiltonian and overlap related integrals
    procedure :: buildDerivativeShift

  end type TTBLite


contains


  !> Constructor for the library interface
  subroutine TTBLite_init(this, input, nAtom, species0, speciesNames, coords0, latVecs)

    !> Instance of the library interface
    type(TTBLite), intent(out) :: this

    !> Input to construct the library interface from
    type(TTBLiteInput), intent(in) :: input

    !> Nr. of atoms in the system
    integer, intent(in) :: nAtom

    !> Species of every atom in the unit cell
    integer, intent(in) :: species0(:)

    !> Atomic coordinates in the unit cell
    real(dp), intent(in) :: coords0(:,:)

    !> Symbols of the species
    character(len=*), intent(in) :: speciesNames(:)

    !> Lattice vectors, if the system is periodic
    real(dp), intent(in), optional :: latVecs(:,:)

  #:if WITH_TBLITE
    type(scf_info) :: info
    character(len=symbol_length), allocatable :: symbol(:)

    symbol = speciesNames(species0)
    call new(this%mol, symbol, coords0, lattice=latVecs)

    select case(input%method%id)
    case default
      call error("Unknown method selector")
    case(gfn2xtb%id)
      call new_gfn2_calculator(this%calc, this%mol)
    case(gfn1xtb%id)
      call new_gfn1_calculator(this%calc, this%mol)
    case(ipea1xtb%id)
      call new_ipea1_calculator(this%calc, this%mol)
    end select

    info = this%calc%variable_info()
    if (info%charge > shell_resolved) then
      call error("Library interface does not support orbital-resolved charge communication")
    end if
    if (info%dipole > atom_resolved) then
      call error("Library interface does not support shell/orbital-resolved dipole moment communication")
    end if
    if (info%quadrupole > atom_resolved) then
      call error("Library interface does not support shell/orbital-resolved quadrupole moment communication")
    end if

    call new_wavefunction(this%wfn, this%mol%nat, this%calc%bas%nsh, this%calc%bas%nao, &
        & 0.0_dp)

    call new_potential(this%pot, this%mol, this%calc%bas)

    if (allocated(this%calc%ncoord)) then
      allocate(this%cn(this%mol%nat))
      allocate(this%dcndr(3, this%mol%nat, this%mol%nat), this%dcndL(3, 3, this%mol%nat))
    end if

    allocate(this%selfenergy(this%calc%bas%nsh), this%dsedcn(this%calc%bas%nsh))

    allocate(this%gradient(3, this%mol%nat))

    allocate(this%sp2id(maxval(species0)))
    call getSpeciesIdentifierMap(this%sp2id, species0, this%mol%id)
  #:else
    call not_implemented_error
  #:endif
  end subroutine TTBLite_init


  subroutine getSpeciesIdentifierMap(sp2id, species, id)

    !> Mapping from species to identifiers
    integer, intent(out) :: sp2id(:)

    !> Element species used in DFTB+
    integer, intent(in) :: species(:)

    !> Element identifiers used in tblite
    integer, intent(in) :: id(:)

    integer :: nSpecies, iAt, iSp, iId
    logical, allocatable :: done(:)

    nSpecies= maxval(species)
    allocate(done(nSpecies))
    done(:) = .false.
    do iAt = 1, size(species)
      iId = id(iAt)
      iSp = species(iAt)
      if (done(iSp)) cycle
      sp2id(iSp) = iId
      done(iSp) = .true.
    end do
  end subroutine getSpeciesIdentifierMap


  !> Write information about library setup
  subroutine writeTBLiteInfo(unit, this)

    !> Formatted unit for output
    integer, intent(in) :: unit

    !> Data structure
    type(TTBLite), intent(in) :: this

    character(len=:), allocatable :: version_string

  #:if WITH_TBLITE
    call get_tblite_version(string=version_string)
    write(unit, '(a, ":", t30, a)') "tblite library version", version_string
  #:else
    call not_implemented_error
  #:endif
  end subroutine writeTBLiteInfo


  !> Update internal stored coordinates
  subroutine updateCoords(this, env, neighList, img2CentCell, coords, species0)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> List of neighbours to atoms
    type(TNeighbourList), intent(in) :: neighList

    !> Image to central cell atom index
    integer, intent(in) :: img2CentCell(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Central cell chemical species
    integer, intent(in) :: species0(:)

  #:if WITH_TBLITE
    real(dp) :: cutoff
    real(dp), allocatable :: lattr(:, :)

    this%mol%xyz(:, :) = coords(:, :this%mol%nat)
    this%erep = 0.0_dp
    this%edisp = 0.0_dp
    this%gradient(:, :) = 0.0_dp
    this%sigma(:, :) = 0.0_dp

    if (allocated(this%calc%repulsion)) then
      cutoff = 25.0_dp
      call get_lattice_points(this%mol%periodic, this%mol%lattice, cutoff, lattr)
      call this%calc%repulsion%get_engrad(this%mol, lattr, cutoff, this%erep, &
          & this%gradient, this%sigma)
    end if

    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%update(this%mol, this%dcache)
      call this%calc%dispersion%get_engrad(this%mol, this%dcache, this%edisp, &
          & this%gradient, this%sigma)
    end if

    call new_potential(this%pot, this%mol, this%calc%bas)
    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%update(this%mol, this%cache)
    end if

    if (allocated(this%calc%ncoord)) then
      call this%calc%ncoord%get_cn(this%mol, this%cn, this%dcndr, this%dcndL)
    end if

    call get_selfenergy(this%calc%h0, this%mol%id, this%calc%bas%ish_at, &
        & this%calc%bas%nsh_id, cn=this%cn, selfenergy=this%selfenergy, dsedcn=this%dsedcn)
  #:else
    call not_implemented_error
  #:endif
  end subroutine updateCoords


  !> Update internal copy of lattice vectors
  subroutine updateLatVecs(this, latVecs)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Lattice vectors
    real(dp), intent(in) :: latVecs(:,:)

  #:if WITH_TBLITE
    this%mol%lattice(:, :) = latVecs
  #:else
    call not_implemented_error
  #:endif
  end subroutine updateLatVecs


  !> Get energy contributions
  subroutine getEnergies(this, energies)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Energy contributions for each atom
    real(dp), intent(out) :: energies(:)

  #:if WITH_TBLITE
    energies(:) = (this%erep + this%edisp + this%escd + this%ees) / size(energies)
  #:else
    call not_implemented_error
  #:endif
  end subroutine getEnergies


  !> Get force contributions
  subroutine addGradients(this, env, neighList, species, coords, img2CentCell, gradients)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Specie for each atom.
    integer, intent(in) :: species(:)

    !> Coordinate of each atom.
    real(dp), intent(in) :: coords(:,:)

    !> Mapping of atoms to cetnral cell.
    integer, intent(in) :: img2CentCell(:)

    !> Gradient contributions for each atom
    real(dp), intent(inout) :: gradients(:,:)

  #:if WITH_TBLITE
    gradients(:, :) = gradients + this%gradient
  #:else
    call not_implemented_error
  #:endif
  end subroutine addGradients


  !> Get stress tensor contributions
  subroutine getStress(this, stress)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Stress tensor contributions
    real(dp), intent(out) :: stress(:,:)

  #:if WITH_TBLITE
    stress(:, :) = this%sigma / abs(determinant33(this%mol%lattice))
  #:else
    call not_implemented_error
  #:endif
  end subroutine getStress


  !> Distance cut off for dispersion interactions
  function getRCutoff(this) result(cutoff)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Resulting cutoff
    real(dp) :: cutoff

  #:if WITH_TBLITE
    cutoff = get_cutoff(this%calc%bas)
  #:else
    call not_implemented_error
  #:endif
  end function getRCutoff


  !> Updates with changed populations for the instance.
  !>
  !> This routine will be called for both potential and energy calculations.
  !> In case of the energy calculations the orbital charges qq are the actual
  !> output charges, while in case of the potential calculation the orbital
  !> charges in qq are an incomplete output from the mixer and are only accurate
  !> up to the requested charges from the variable_info routine.
  !>
  !> This routine gets orbital populations as input, rather than charges, therefore
  !> we are working with the opposite sign convention w.r.t. the library.
  subroutine updatePopulations(this, env, species, neighList, img2CentCell, orb, &
      & qq, q0, dipAtom, quadAtom)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Species, shape: [nAtom]
    integer, intent(in) :: species(:)

    !> Neighbour list.
    type(TNeighbourList), intent(in) :: neighList

    !> Orbital populations.
    real(dp), intent(in) :: qq(:,:,:)

    !> Reference orbital populations.
    real(dp), intent(in) :: q0(:,:,:)

    !> Dipole populations for each atom
    real(dp), intent(in), optional :: dipAtom(:,:)

    !> Quadrupole populations for each atom
    real(dp), intent(in), optional :: quadAtom(:,:)

    !> Mapping on atoms in central cell.
    integer, intent(in) :: img2CentCell(:)

    !> Orbital information
    type(TOrbitals), intent(in) :: orb

  #:if WITH_TBLITE
    integer :: iAt, iSh, ii
    real(dp), allocatable :: dQAtom(:), dQShell(:, :)

    call this%pot%reset
    this%escd = 0.0_dp
    this%ees = 0.0_dp

    allocate(dQAtom(this%mol%nat), dQShell(orb%mShell, this%mol%nat))
    call getSummedCharges(species, orb, qq, q0, dQAtom=dQAtom, dQShell=dQShell)

    this%wfn%qat(:) = -dQAtom
    do iAt = 1, size(dQShell, 2)
      ii = this%calc%bas%ish_at(iAt)
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        this%wfn%qsh(ii+iSh) = -dQShell(iSh, iAt)
      end do
    end do

    if (present(dipAtom)) then
      this%wfn%dpat(:, :) = -dipAtom
    end if

    if (present(quadAtom)) then
      this%wfn%qpat(:, :) = -quadAtom
    end if

    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%get_potential(this%mol, this%cache, this%wfn, this%pot)
    end if
    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%get_potential(this%mol, this%dcache, this%wfn, this%pot)
    end if

    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%get_energy(this%mol, this%cache, this%wfn, this%ees)
    end if
    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%get_energy(this%mol, this%dcache, this%wfn, this%escd)
    end if
  #:else
    call not_implemented_error
  #:endif
  end subroutine updatePopulations


  !> Returns atom-resolved and shell-resolved shifts from library.
  !>
  !> Potential shifts in tblite are calculated as derivatives w.r.t. the partial
  !> charges / multipole moments, while in DFTB+ they are calculated as derivative
  !> w.r.t. the population. This means we have to flip the sign here.
  subroutine getShifts(this, shiftPerAtom, shiftPerShell, shiftPerDipole, shiftPerQuadrupole)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Shift per atom
    real(dp), intent(out) :: shiftPerAtom(:)

    !> Shift per shell
    real(dp), intent(out) :: shiftPerShell(:,:)

    !> Shift per dipole moment
    real(dp), intent(out) :: shiftPerDipole(:,:)

    !> Shift per quadrupole moment
    real(dp), intent(out) :: shiftPerQuadrupole(:,:)

  #:if WITH_TBLITE
    integer :: iAt, iSh, ii

    shiftPerAtom(:) = -this%pot%vat

    shiftPerShell(:,:) = 0.0_dp
    do iAt = 1, size(shiftPerShell, 2)
      ii = this%calc%bas%ish_at(iAt)
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        shiftPerShell(iSh, iAt) = -this%pot%vsh(ii+iSh)
      end do
    end do

    shiftPerDipole(:, :) = -this%pot%vdp
    shiftPerQuadrupole(:, :) = -this%pot%vqp
  #:else
    call not_implemented_error
  #:endif
  end subroutine getShifts


  !> Create orbital information. The orbital information is already generated in
  !> this%calc%bas, but might be incomplete w.r.t. the information required here.
  !>
  !> Both this%mol%id and species are identifiers for unique groups of elements,
  !> since we cannot guarantee that the species in DFTB+ are identical to the ones
  !> found in the library we will reconstruct the species information from the atoms.
  subroutine getOrbitalInfo(this, species0, orb)

    !> Data structure
    class(TTBLite), intent(in) :: this

    !> Species of each atom, shape: [nAtom]
    integer, intent(in) :: species0(:)

    !> Orbital information
    type(TOrbitals), intent(out) :: orb

  #:if WITH_TBLITE
    integer :: nShell, mShell, nSpecies, iAt, iSp, iId, iSh, ind, ii
    logical, allocatable :: done(:)

    nSpecies= maxval(species0)
    mShell = maxval(this%calc%bas%nsh_id)
    allocate(done(nSpecies))
    allocate(orb%angShell(mShell, nSpecies))
    done(:) = .false.
    do iAt = 1, size(species0)
      iId = this%mol%id(iAt)
      iSp = species0(iAt)
      if (done(iSp)) cycle
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        orb%angShell(iSh, iSp) = this%calc%bas%cgto(iSh, iId)%ang
      end do
      done(iSp) = .true.
    end do

    allocate(orb%nShell(nSpecies))
    allocate(orb%nOrbSpecies(nSpecies))
    allocate(orb%nOrbAtom(size(species0)))

    done(:) = .false.
    do iAt = 1, size(species0)
      iSp = species0(iAt)
      ii = this%calc%bas%ish_at(iAt)
      nShell = this%calc%bas%nsh_at(iAt)
      orb%nOrbAtom(iAt) = sum(this%calc%bas%nao_sh(ii+1:ii+nShell))
      if (done(iSp)) cycle
      orb%nOrbSpecies(iSp) = orb%nOrbAtom(iAt)
      orb%nShell(iSp) = nShell
      done(iSp) = .true.
    end do

    orb%mShell = mShell
    orb%mOrb = maxval(orb%nOrbSpecies)
    orb%nOrb = sum(orb%nOrbAtom)

    allocate(orb%iShellOrb(orb%mOrb, nSpecies))
    allocate(orb%posShell(orb%mShell+1, nSpecies))
    do iSp = 1, nSpecies
      ind = 1
      do iSh = 1, orb%nShell(iSp)
        orb%posShell(iSh, iSp) = ind
        orb%iShellOrb(ind:ind+2*orb%angShell(iSh, iSp), iSp) = iSh
        ind = ind + 2 * orb%angShell(iSh, iSp) + 1
      end do
      orb%posShell(orb%nShell(iSp)+1, iSp) = ind
    end do
  #:else
    call not_implemented_error
  #:endif
  end subroutine getOrbitalInfo


  !> Get reference occupation numbers.
  subroutine getReferenceN0(this, species0, referenceN0)

    !> Data structure
    class(TTBLite), intent(in) :: this

    !> Species of each atom, shape: [nAtom]
    integer, intent(in) :: species0(:)

    !> Reference occupation numbers
    real(dp), intent(out) :: referenceN0(:, :)

  #:if WITH_TBLITE
    integer :: iAt, iSp, iId, iSh
    logical, allocatable :: done(:)

    referenceN0(:,:) = 0.0_dp
    allocate(done(maxval(species0)))
    done(:) = .false.
    do iAt = 1, size(species0)
      iId = this%mol%id(iAt)
      iSp = species0(iAt)
      if (done(iSp)) cycle
      do iSh = 1, this%calc%bas%nsh_at(iAt)
        referenceN0(iSh, iSp) = this%calc%h0%refocc(iSh, iId)
      end do
      done(iSp) = .true.
    end do
  #:else
    call not_implemented_error
  #:endif
  end subroutine getReferenceN0


  !> Returns the equivalence to get the correct mixing of charge dependent contributions
  subroutine getOrbitalEquiv(this, equiv, requireDipole, requireQuadrupole, orb, species)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> The equivalence vector on return
    integer, intent(out) :: equiv(:,:,:)

    !> Require mixing for cumulative atomic dipole moments
    logical, intent(out) :: requireDipole

    !> Require mixing for cumulative atomic quadrupole moments
    logical, intent(out) :: requireQuadrupole

    !> Information about the orbitals and their angular momenta
    type(TOrbitals), intent(in) :: orb

    !> Species of each atom
    integer, intent(in) :: species(:)

  #:if WITH_TBLITE
    type(scf_info) :: info
    integer :: nAtom, iCount, iSpin, nSpin, iAt, iSp, ii

    nAtom = size(equiv, dim=2)
    nSpin = size(equiv, dim=3)

    equiv(:,:,:) = 0

    info = this%calc%variable_info()
    select case(info%charge)
    case(atom_resolved)
      iCount = 0
      do iSpin = 1, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          iCount = iCount + 1
          do ii = 1, orb%nOrbSpecies(iSp)
            equiv(ii, iAt, iSpin) = iCount
          end do
        end do
      end do
    case(shell_resolved)
      iCount = 0
      do iSpin = 1, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          do ii = 1, orb%nOrbSpecies(iSp)
            equiv(ii, iAt, iSpin) = iCount + orb%iShellOrb(ii, iSp)
          end do
          iCount = iCount + orb%nShell(iSp)
        end do
      end do
    case(orbital_resolved)
      iCount = 0
      do iSpin = 1, nSpin
        do iAt = 1, nAtom
          iSp = species(iAt)
          do ii = 1, orb%nOrbSpecies(iSp)
            iCount = iCount + 1
            equiv(ii, iAt, iSpin) = iCount
          end do
        end do
      end do
    end select

    requireDipole = info%dipole == atom_resolved
    requireQuadrupole = info%quadrupole == atom_resolved
  #:else
    call not_implemented_error
  #:endif
  end subroutine getOrbitalEquiv


  !> Build atomic block sparse compressed Hamiltonian and overlap related integrals
  subroutine buildSH0(this, env, species, coords, nNeighbour, iNeighbours, img2centCell, &
      & iPair, orb, hamiltonian, overlap, dpintBra, dpintKet, qpintBra, qpintKet)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: hamiltonian(:)

    !> Returns the non-self-consistent Hamiltonian
    real(dp), intent(out) :: overlap(:)

    !> Dipole moment integral matrix
    real(dp), intent(inout) :: dpintBra(:, :), dpintKet(:, :)

    !> Quadrupole moment integral matrix
    real(dp), intent(inout) :: qpintBra(:, :), qpintKet(:, :)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Mapping of images back to atoms in the central cell
    integer, intent(in) :: img2centCell(:)

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms
    integer, intent(in) :: iPair(0:,:)

    !> Information about the orbitals in the system
    type(TOrbitals), intent(in) :: orb

  #:if WITH_TBLITE
    integer :: nAtom, iAtFirst, iAtLast

    nAtom = size(nNeighbour)
    hamiltonian(:) = 0.0_dp
    overlap(:) = 0.0_dp
    dpintBra(:, :) = 0.0_dp
    dpintKet(:, :) = 0.0_dp
    qpintBra(:, :) = 0.0_dp
    qpintKet(:, :) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    call buildDiagonalBlocks(iAtFirst, iAtLast, this%mol%id, iPair, orb%nOrbAtom, &
        & this%calc%bas, this%selfenergy, hamiltonian, overlap, dpintBra, dpintKet, &
        & qpintBra, qpintKet)

    call buildDiatomicBlocks(iAtFirst, iAtLast, this%mol%id, coords, &
        & nNeighbour, iNeighbours, img2centCell, iPair, orb%nOrbAtom, &
        & this%calc%bas, this%calc%h0, this%selfenergy, hamiltonian, overlap, &
        & dpintBra, dpintKet, qpintBra, qpintKet)

    call assembleChunks(env, hamiltonian)
    call assembleChunks(env, overlap)
  #:else
    call not_implemented_error
  #:endif
  end subroutine buildSH0


#:if WITH_TBLITE
  !> Put the on-site energies into the Hamiltonian,
  !> and <lm|l'm'> = delta_l,l' * delta_m,m' for the overlap
  !>
  !> This routine requires access to the internals of the tblite library,
  !> breaking the abstraction layer of the library interface.
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiagonalBlocks(iAtFirst, iAtLast, species, iPair, nOrbAtom, &
      & bas, selfenergy, hamiltonian, overlap, dpintBra, dpintKet, qpintBra, qpintKet)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Diagonal elememts of the Hamiltonian
    real(dp), intent(in) :: selfenergy(:)

    !> Overlap integral matrix
    real(dp), intent(inout) :: overlap(:)

    !> Dipole moment integral matrix
    real(dp), intent(inout) :: dpintBra(:, :), dpintKet(:, :)

    !> Quadrupole moment integral matrix
    real(dp), intent(inout) :: qpintBra(:, :), qpintKet(:, :)

    !> Effective Hamiltonian
    real(dp), intent(inout) :: hamiltonian(:)

    integer :: iAt, iZp, ind, iOrb, ii, jj, nBlk, is, io, ish, jsh, iao, jao, nao
    integer :: ij, isparse
    real(dp), parameter :: r2 = 0.0_dp, vec(3) = 0.0_dp
    real(dp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)

    allocate(stmp(msao(bas%maxl)**2))
    allocate(dtmp(3, msao(bas%maxl)**2), qtmp(6, msao(bas%maxl)**2))

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(iAtFirst, iAtLast, species, bas, iPair, nOrbAtom, selfenergy, &
    !$omp& overlap, hamiltonian, dpintBra, dpintKet, qpintBra, qpintKet) &
    !$omp private(iAt, iZp, is, io, ind, ish, jsh, ii, jj, iao, jao, nBlk, nao, &
    !$omp& stmp, dtmp, qtmp, ij, isparse)
    do iAt = iAtFirst, iAtLast
      iZp = species(iAt)
      is = bas%ish_at(iAt)
      io = bas%iao_sh(is+1)
      ind = iPair(0, iAt)
      nBlk = nOrbAtom(iAt)
      do iSh = 1, bas%nsh_id(iZp)
        ii = bas%iao_sh(is+iSh) - io
        do iao = 1, msao(bas%cgto(iSh, iZp)%ang)
          isparse = ind + ii+iao + nBlk*(ii+iao-1)

          overlap(isparse) = 1.0_dp
          hamiltonian(isparse) = selfenergy(is+iSh)
        end do
      end do

      do iSh = 1, bas%nsh_id(iZp)
        ii = bas%iao_sh(is+iSh) - io
        do jSh = 1, bas%nsh_id(iZp)
          jj = bas%iao_sh(is+jSh) - io
          call multipole_cgto(bas%cgto(jSh, iZp), bas%cgto(iSh, iZp), &
              & r2, vec, bas%intcut, stmp, dtmp, qtmp)

          nao = msao(bas%cgto(jSh, iZp)%ang)
          do iao = 1, msao(bas%cgto(iSh, iZp)%ang)
            do jao = 1, nao
              ij = jao + nao*(iao-1)
              isparse = ind + jj+jao + nBlk*(ii+iao-1)

              dpintBra(:, isparse) = dtmp(:, ij)
              dpintKet(:, isparse) = dtmp(:, ij)

              qpintBra(:, isparse) = qtmp(:, ij)
              qpintKet(:, isparse) = qtmp(:, ij)
            end do
          end do
        end do
      end do
    end do
  end subroutine buildDiagonalBlocks


  !> The Hamiltonian is saved in an atomic block sparse compressed format.
  !> We calculate a shell pair as a contiguous blocks and spread it to the
  !> contiguous atomic block.
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiatomicBlocks(iAtFirst, iAtLast, species, coords, &
      & nNeighbour, iNeighbours, img2centCell, iPair, nOrbAtom, bas, h0, selfenergy, &
      & hamiltonian, overlap, dpintBra, dpintKet, qpintBra, qpintKet)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Mapping of images back to atoms in the central cell
    integer, intent(in) :: img2centCell(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Hamiltonian interaction data
    type(tb_hamiltonian), intent(in) :: h0

    !> Diagonal elememts of the Hamiltonian
    real(dp), intent(in) :: selfenergy(:)

    !> Overlap integral matrix
    real(dp), intent(inout) :: overlap(:)

    !> Dipole moment integral matrix
    real(dp), intent(inout) :: dpintBra(:, :), dpintKet(:, :)

    !> Quadrupole moment integral matrix
    real(dp), intent(inout) :: qpintBra(:, :), qpintKet(:, :)

    !> Effective Hamiltonian
    real(dp), intent(inout) :: hamiltonian(:)

    integer :: iAt, jAt, iZp, jZp, iNeigh, img, ind, io, jo, ij, isparse
    integer :: iSh, jSh, is, js, ii, jj, iao, jao, nao, nBlk
    real(dp) :: rr, r2, vec(3), hij, shpoly, dtmpj(3), qtmpj(6), sval
    real(dp), allocatable :: stmp(:)
    real(dp), allocatable :: dtmpi(:, :), qtmpi(:, :)

    allocate(stmp(msao(bas%maxl)**2))
    allocate(dtmpi(3, msao(bas%maxl)**2), qtmpi(6, msao(bas%maxl)**2))

    !$omp parallel do schedule(runtime) default(none) &
    !$omp shared(iatfirst, iatlast, species, bas, overlap, hamiltonian, h0, selfenergy, &
    !$omp& nNeighbour, ineighbours, img2centCell, ipair, nOrbAtom, coords, &
    !$omp& dpintBra, dpintKet, qpintBra, qpintKet) &
    !$omp private(iAt, jAt, iZp, jZp, is, js, iSh, jSh, ii, jj, iao, jao, nao, nBlk, ij, &
    !$omp& io, jo, iNeigh, img, ind, r2, vec, stmp, dtmpi, dtmpj, qtmpi, qtmpj, hij, &
    !$omp& shpoly, rr, sval, isparse)
    do iAt = iAtFirst, iAtLast
      iZp = species(iAt)
      is = bas%ish_at(iAt)
      io = bas%iao_sh(is+1)
      do iNeigh = 1, nNeighbour(iAt)
        img = iNeighbours(iNeigh, iAt)
        jAt = img2centCell(img)
        jZp = species(jAt)
        js = bas%ish_at(jAt)
        ind = iPair(iNeigh, iAt)
        jo = bas%iao_sh(js+1)
        nBlk = nOrbAtom(jAt)

        vec(:) = coords(:, iAt) - coords(:, img)
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        rr = sqrt(sqrt(r2) / (h0%rad(jZp) + h0%rad(iZp)))
        do ish = 1, bas%nsh_id(iZp)
          ii = bas%iao_sh(is+iSh) - io
          do jsh = 1, bas%nsh_id(jZp)
            jj = bas%iao_sh(js+jSh) - jo
            ! Moment operator is placed on ket (iAt) atom
            call multipole_cgto(bas%cgto(jSh, jZp), bas%cgto(iSh, iZp), &
                & r2, vec, bas%intcut, stmp, dtmpi, qtmpi)

            shpoly = (1.0_dp + h0%shpoly(iSh, iZp)*rr) &
                & * (1.0_dp + h0%shpoly(jSh, jZp)*rr)

            hij = 0.5_dp * (selfenergy(is+iSh) + selfenergy(js+jSh)) &
                & * h0%hscale(jSh, iSh, jZp, iZp) * shpoly

            nao = msao(bas%cgto(jSh, jZp)%ang)
            do iao = 1, msao(bas%cgto(iSh, iZp)%ang)
              do jao = 1, nao
                ij = jao + nao*(iao-1)
                isparse = ind + jj+jao + nBlk*(ii+iao-1)
                sval = stmp(ij)
                call shiftOperator(qtmpj, dtmpj, qtmpi(:, ij), dtmpi(:, ij), sval, vec)

                overlap(isparse) = stmp(ij)
                hamiltonian(isparse) = stmp(ij) * hij

                dpintKet(:, isparse) = dtmpi(:, ij)
                dpintBra(:, isparse) = dtmpj

                qpintKet(:, isparse) = qtmpi(:, ij)
                qpintBra(:, isparse) = qtmpj
              end do
            end do

          end do
        end do

      end do
    end do
  end subroutine buildDiatomicBlocks


  !> Perform the shift of the moment operator from the Ket (iAt) to the Bra (jAt) atom.
  pure subroutine shiftOperator(qj, dj, qi, di, s, vec)
    !> Quadrupole moment integral with operator on Bra (jAt) atom
    real(dp), intent(out) :: qj(:)
    !> Dipole moment integral with operator on Bra (jAt) atom
    real(dp), intent(out) :: dj(:)
    !> Quadrupole moment integral with operator on Ket (iAt) atom
    real(dp), intent(in) :: qi(:)
    !> Dipole moment integral with operator on Ket (iAt) atom
    real(dp), intent(in) :: di(:)
    !> Overlap integral
    real(dp), intent(in) :: s
    !> Displacement vector between i and j
    real(dp), intent(in) :: vec(:)

    real(dp) :: qtr

    ! Shift dipole operator from ket (iAt) to bra (jAt) atom
    dj(:) = di - vec*s
    ! Quadrupole integrals on ket (iAt) are already traceless,
    ! therefore we create only the shift term for the operator first
    qj([1, 3, 6]) = -2*vec*di + vec**2*s
    qj([2, 4, 5]) = &
      & - vec([1, 1, 2])*di([2, 3, 3]) &
      & - vec([2, 3, 3])*di([1, 1, 2]) &
      & + vec([1, 1, 2])*vec([2, 3, 3])*s
    ! Collect trace from the shift terms
    qtr = 0.5_dp * (qj(1) + qj(3) + qj(6))
    ! Rescale and remove the trace from the shift term, before creating
    ! the full quadrupole integral on the bra (jAt) atom
    qj(1) = qi(1) + 1.5_dp * qj(1) - qtr
    qj(2) = qi(2) + 1.5_dp * qj(2)
    qj(3) = qi(3) + 1.5_dp * qj(3) - qtr
    qj(4) = qi(4) + 1.5_dp * qj(4)
    qj(5) = qi(5) + 1.5_dp * qj(5)
    qj(6) = qi(6) + 1.5_dp * qj(6) - qtr
  end subroutine shiftOperator
#:endif


  !> Evaluate derivatives of potential shifts from Hamiltonian and overlap related
  !> integrals.
  subroutine buildDerivativeShift(this, env, DM, EDM, coords, species, &
      & nNeighbour, iNeighbours, img2CentCell, iPair, orb, shift)

    !> Data structure
    class(TTBLite), intent(inout) :: this

    !> Computational environment settings
    type(TEnvironment), intent(in) :: env

    !> density matrix in packed format
    real(dp), intent(in) :: DM(:,:)

    !> energy-weighted density matrix in packed format
    real(dp), intent(in) :: EDM(:)

    !> list of all atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> list of all atomic species
    integer, intent(in) :: species(:)

    !> number of neighbours of each atom
    integer, intent(in) :: nNeighbour(:)

    !> neighbour list for atoms
    integer, intent(in) :: iNeighbours(0:,:)

    !> indexing array for periodic image atoms
    integer, intent(in) :: img2CentCell(:)

    !> indexing array for the Hamiltonian
    integer, intent(in) :: iPair(0:,:)

    !> Information about the shells and orbitals in the system.
    type(TOrbitals), intent(in) :: orb

    !> block shift from the potential
    real(dp), intent(in) :: shift(:,:,:,:)

  #:if WITH_TBLITE
    integer :: nAtom, iAtFirst, iAtLast
    real(dp), allocatable :: gradient(:, :), sigma(:, :), dEdcn(:)

    if (allocated(this%calc%coulomb)) then
      call this%calc%coulomb%get_gradient(this%mol, this%cache, this%wfn, this%gradient, &
          & this%sigma)
    end if

    if (allocated(this%calc%dispersion)) then
      call this%calc%dispersion%get_gradient(this%mol, this%dcache, this%wfn, this%gradient, &
          & this%sigma)
    end if

    nAtom = size(nNeighbour)
    allocate(gradient(3, nAtom), sigma(3, 3), dEdcn(nAtom))
    gradient(:, :) = 0.0_dp
    sigma(:, :) = 0.0_dp
    dEdcn(:) = 0.0_dp

    call distributeRangeInChunks(env, 1, nAtom, iAtFirst, iAtLast)

    call buildDiagonalDerivs(iAtFirst, iAtLast, this%mol%id, iPair, orb%nOrbAtom, &
        & this%calc%bas, this%dsedcn, DM, dEdcn)
    call buildDiatomicDerivs(iAtFirst, iAtLast, this%mol%id, coords, nNeighbour, &
        & iNeighbours, img2centCell, iPair, orb%nOrbAtom, this%calc%bas, this%calc%h0, &
        & this%selfenergy, this%dsedcn, DM, EDM, shift, dEdcn, gradient, sigma)

    call assembleChunks(env, dEdcn)
    call assembleChunks(env, gradient)
    call assembleChunks(env, sigma)

    call gemv(gradient, this%dcndr, dEdcn, beta=1.0_dp)
    call gemv(sigma, this%dcndL, dEdcn, beta=1.0_dp)

    this%gradient(:, :) = this%gradient + gradient
    this%sigma(:, :) = this%sigma + sigma
  #:else
    call not_implemented_error
  #:endif
  end subroutine buildDerivativeShift


#:if WITH_TBLITE
  !> Calculate derivatives for all diagonal elements.
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiagonalDerivs(iAtFirst, iAtLast, species, iPair, nOrbAtom, &
      & bas, dsedcn, pmat, dEdcn)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Derivative of diagonal elememts of the Hamiltonian w.r.t. coordination number
    real(dp), intent(in) :: dsedcn(:)

    !> Density matrix
    real(dp), intent(in) :: pmat(:, :)

    !> Derivative of energy w.r.t. coordination number
    real(dp), intent(inout) :: dEdcn(:)

    integer :: iat, izp, is, ish, ii, iao, io, ind, nBlk
    real(dp) :: dhdcni, dcni

    !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn) &
    !$omp shared(iAtFirst, iAtLast, species, bas, dsedcn, pmat, iPair, nOrbAtom) &
    !$omp private(iat, izp, is, ish, ii, iao, io, ind, nBlk, dcni, dhdcni)
    do iat = iAtFirst, iAtLast
      izp = species(iat)
      is = bas%ish_at(iat)
      io = bas%iao_sh(is+1)
      ind = iPair(0, iAt)
      nBlk = nOrbAtom(iAt)
      do ish = 1, bas%nsh_id(izp)
        ii = bas%iao_sh(is+ish) - io
        dhdcni = dsedcn(is+ish)
        dcni = 0.0_dp
        do iao = 1, msao(bas%cgto(ish, izp)%ang)
          dcni = dcni + dhdcni * pmat(ind + ii+iao + nBlk*(ii+iao-1), 1)
        end do
        dEdcn(iat) = dEdcn(iat) + dcni
      end do
    end do
  end subroutine buildDiagonalDerivs


  !> Calculate derivatives for all diatomic pairs.
  !>
  !> Todo: Handle actual block potential shifts
  !>
  !> Candidate for (partial) upstreaming in tblite library.
  subroutine buildDiatomicDerivs(iAtFirst, iAtLast, species, coords, &
      & nNeighbour, iNeighbours, img2centCell, iPair, nOrbAtom, bas, h0, &
      & selfenergy, dsedcn, pmat, xmat, shift, dEdcn, gradient, sigma)

    !> Atom range for this processor to evaluate
    integer, intent(in) :: iAtFirst, iAtLast

    !> Chemical species of each atom
    integer, intent(in) :: species(:)

    !> Atomic coordinates
    real(dp), intent(in) :: coords(:,:)

    !> Number of surrounding neighbours for each atom
    integer, intent(in) :: nNeighbour(:)

    !> List of surrounding neighbours for each atom
    integer, intent(in) :: iNeighbours(0:,:)

    !> Shift vector, where the interaction between two atoms starts in the sparse format.
    integer, intent(in) :: iPair(0:,:)

    !> Size of the block in spare format for each atom
    integer, intent(in) :: nOrbAtom(:)

    !> Mapping of images back to atoms in the central cell
    integer, intent(in) :: img2centCell(:)

    !> Basis set information
    type(basis_type), intent(in) :: bas

    !> Hamiltonian interaction data
    type(tb_hamiltonian), intent(in) :: h0

    !> Diagonal elememts of the Hamiltonian
    real(dp), intent(in) :: selfenergy(:)

    !> Derivative of diagonal elememts of the Hamiltonian w.r.t. coordination number
    real(dp), intent(in) :: dsedcn(:)

    !> Density matrix in packed form with spin information
    real(dp), intent(in) :: pmat(:, :)

    !> Energy weighted density matrix in packed form
    real(dp), intent(in) :: xmat(:)

    !> Block potential shift, shape: [nOrb, nOrb, nAtom, nSpin]
    real(dp), intent(in) :: shift(:, :, :, :)

    !> Derivative of energy w.r.t. coordination number
    real(dp), intent(inout) :: dEdcn(:)

    !> Energy derivative w.r.t. atomic displacements
    real(dp), intent(inout) :: gradient(:, :)

    !> Energy derivative w.r.t. strain deformations
    real(dp), intent(inout) :: sigma(:, :)

    integer :: iat, jat, izp, jzp, itr, io, jo, nblk, iNeigh, img, ind
    integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
    real(dp) :: rr, r2, vec(3), cutoff2, hij, shpoly, dshpoly, dG(3), hscale
    real(dp) :: sval, dcni, dcnj, dhdcni, dhdcnj, hpij, pij
    real(dp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)
    real(dp), allocatable :: dstmp(:, :), ddtmpi(:, :, :), dqtmpi(:, :, :)
    real(dp), allocatable :: ddtmpj(:, :, :), dqtmpj(:, :, :)

    allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & dtmp(3, msao(bas%maxl)**2), ddtmpi(3, 3, msao(bas%maxl)**2), &
      & qtmp(6, msao(bas%maxl)**2), dqtmpi(3, 6, msao(bas%maxl)**2), &
      & ddtmpj(3, 3, msao(bas%maxl)**2), dqtmpj(3, 6, msao(bas%maxl)**2))

    !$omp parallel do schedule(runtime) default(none) reduction(+:dEdcn, gradient, sigma) &
    !$omp shared(iAtFirst, iAtLast, species, bas, h0, selfenergy, dsedcn, coords, &
    !$omp& nNeighbour, iNeighbours, img2centCell, iPair, nOrbAtom, pmat, xmat, shift) &
    !$omp private(iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, ind, &
    !$omp& iNeigh, io, jo, nblk, img, r2, vec, stmp, dtmp, qtmp, dstmp, ddtmpi, ddtmpj, &
    !$omp& dqtmpi, dqtmpj, hij, shpoly, dshpoly, dG, dcni, dcnj, dhdcni, dhdcnj, hpij, rr, &
    !$omp& sval, hscale, pij)
    do iat = iAtFirst, iAtLast
      izp = species(iat)
      is = bas%ish_at(iat)
      io = bas%iao_sh(is+1)
      do iNeigh = 1, nNeighbour(iAt)
        img = iNeighbours(iNeigh, iAt)
        jAt = img2centCell(img)
        jZp = species(jAt)
        js = bas%ish_at(jAt)
        ind = iPair(iNeigh, iAt)
        jo = bas%iao_sh(js+1)
        nBlk = nOrbAtom(jAt)

        vec(:) = coords(:, iAt) - coords(:, img)
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        rr = sqrt(sqrt(r2) / (h0%rad(jzp) + h0%rad(izp)))
        do ish = 1, bas%nsh_id(izp)
          ii = bas%iao_sh(is+ish) - io
          do jsh = 1, bas%nsh_id(jzp)
            jj = bas%iao_sh(js+jsh) - jo
            !call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
            !    & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddtmpj, dqtmpj, &
            !    & ddtmpi, dqtmpi)
            call overlap_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                & r2, vec, bas%intcut, stmp, dstmp)

            shpoly = (1.0_dp + h0%shpoly(ish, izp)*rr) &
              & * (1.0_dp + h0%shpoly(jsh, jzp)*rr)
            dshpoly = ((1.0_dp + h0%shpoly(ish, izp)*rr)*h0%shpoly(jsh, jzp)*rr &
              & + (1.0_dp + h0%shpoly(jsh, jzp)*rr)*h0%shpoly(ish, izp)*rr) &
              & * 0.5_dp / r2

            hscale = h0%hscale(jsh, ish, jzp, izp)
            hij = 0.5_dp * (selfenergy(is+ish) + selfenergy(js+jsh)) * hscale
            dhdcni = dsedcn(is+ish) * shpoly * hscale
            dhdcnj = dsedcn(js+jsh) * shpoly * hscale

            dG(:) = 0.0_dp
            dcni = 0.0_dp
            dcnj = 0.0_dp
            nao = msao(bas%cgto(jsh, jzp)%ang)
            do iao = 1, msao(bas%cgto(ish, izp)%ang)
              do jao = 1, nao

                ij = jao + nao*(iao-1)
                pij = pmat(ind + jj+jao + nBlk*(ii+iao-1), 1)
                hpij = pij * hij * shpoly
                sval = 2*hpij - 2*xmat(ind + jj+jao + nBlk*(ii+iao-1)) &
                  & + pij * (shift(jj+jao, jj+jao, jat, 1) + shift(ii+iao, ii+iao, iat, 1))

                dG(:) = dG + sval * dstmp(:, ij) &
                  & + 2*hpij*stmp(ij) * dshpoly / shpoly * vec
                  !& + pij * matmul(ddtmpi(:, :, ij), pot%vdp(:, iat)) &
                  !& + pij * matmul(ddtmpj(:, :, ij), pot%vdp(:, jat)) &
                  !& + pij * matmul(dqtmpi(:, :, ij), pot%vqp(:, iat)) &
                  !& + pij * matmul(dqtmpj(:, :, ij), pot%vqp(:, jat))

                dcni = dcni + dhdcni * pij * stmp(ij)
                dcnj = dcnj + dhdcnj * pij * stmp(ij)
              end do
            end do
            dEdcn(iat) = dEdcn(iat) + dcni
            dEdcn(jat) = dEdcn(jat) + dcnj
            gradient(:, iat) = gradient(:, iat) + dG
            gradient(:, jat) = gradient(:, jat) - dG
            sigma(:, :) = sigma + 0.5_dp * (spread(vec, 1, 3) * spread(dG, 2, 3) &
              & + spread(dG, 1, 3) * spread(vec, 2, 3))

          end do
        end do

      end do
    end do
  end subroutine buildDiatomicDerivs
#:endif


#:if not WITH_TBLITE
  subroutine not_implemented_error

    call error("DFTB+ compiled without support for tblite library")
  end subroutine not_implemented_error
#:endif


end module dftbp_tblite
