module MicroStrucFeatures

  use GrainIdentification, only: Grain
  use ReferenceLattice

  implicit none
!***********************************************************************!

  Type StructFeatures

    type(Grain) :: grains
    type(RefLattice) :: reference

  end type StructFeatures

  type(StructFeatures) :: features

end module MicroStrucFeatures