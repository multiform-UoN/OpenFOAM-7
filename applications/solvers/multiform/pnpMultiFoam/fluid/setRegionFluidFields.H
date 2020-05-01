//Set current mesh to fluid region i
fvMesh& mesh = fluidRegions[i];

//Set current fields to fluid region i
volScalarField& V = VFluid[i];
volScalarField& T = TFluid[i];
surfaceScalarField& phiNP = phiNPFluid[i];
volVectorField& U = UFluid[i];
volScalarField& rho = rhoFluid[i];

//Creates and initialises the relative face-flux field phi.
#include "createPhi.H"

//Set current concentrations to fluid region i
PtrList<volScalarField>& cPtrL = cPtrLFluid[i];

//Set species dictionary to fluid region i
PtrList<dictionary>& species = speciesFluid[i];

//Set effective dielectric to fluid region i
dimensionedScalar& epsilon_ = epsilon_Fluid[i];

//Set electro-neutrality and diffusion dominant flags to fluid region i
const bool&	electroNeutrality = electroNeutralityFluid[i];
const bool& diffusionDominant = diffusionDominantFluid[i];

//Set reference cell and value to fluid region i
const label& VRefCell = VRefCellFluid[i];
const scalar& VRefValue = VRefValueFluid[i];


//Find # of loop controls set for fluid region i
pimpleNoLoopControl& pimple = pimples.pimple(i);
