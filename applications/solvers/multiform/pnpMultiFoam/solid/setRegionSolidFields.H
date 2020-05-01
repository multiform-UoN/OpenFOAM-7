fvMesh& mesh = solidRegions[i];

//Set current fields to solid region i
volScalarField& V = VSolid[i];
volScalarField& T = TSolid[i];
surfaceScalarField& phiNP = phiNPSolid[i];
volVectorField& U = USolid[i];
volScalarField& rho = rhoSolid[i];

//Creates and initialises the relative face-flux field phi.
#include "createPhi.H"

//Set current concentrations to solid region i
PtrList<volScalarField>& cPtrL = cPtrLSolid[i];

//Set species dictionary to solid region i
PtrList<dictionary>& species = speciesSolid[i];

//Set effective dielectric to solid region i
dimensionedScalar& epsilon_ = epsilon_Solid[i];

//Set electro-neutrality and diffusion dominant flags to fluid region i
const bool&	electroNeutrality = electroNeutralitySolid[i];
const bool& diffusionDominant = diffusionDominantSolid[i];


//Set reference cell and value to fluid region i
const label& VRefCell = VRefCellSolid[i];
const scalar& VRefValue = VRefValueSolid[i];

//Find # of loop controls set for solid region i
//solidNoLoopControl& pimple = pimples.solid(i);
