module MicroStrucFeatures

  use GrainIdentification, only: Grain

  implicit none
!***********************************************************************!

  Type StructFeatures

    type(Grain) :: grains

  end type StructFeatures

  type(StructFeatures) :: features

end module MicroStrucFeatures