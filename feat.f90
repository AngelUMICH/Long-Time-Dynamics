module MicroStrucFeatures

  use GrainIdentification, only: Grain
  use ReferenceLattice
  use Clusters

  implicit none
!***********************************************************************!

  Type StructFeatures

    type(Grain) :: grains
    type(RefLattice) :: reference
    type(DefectClusters) :: defects

  end type StructFeatures

!***********************************************************************!

  type(StructFeatures) :: features

end module MicroStrucFeatures