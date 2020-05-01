#include "createRDeltaT.H"

//- V is the electric potential
Info<< "Reading field V\n" << endl;
volScalarField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//- Read temperature field
Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating Nernst-Planck face flux field phiNP\n" << endl;

surfaceScalarField phiNP
(
    IOobject
    (
        "phiNP",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -fvc::flux(e_*fvc::grad(V)/(k_*T))
);



// Initialize reference value and cell
label VRefCell = 0;
scalar VRefValue = 0.0;

//Update RefCell and RefValue if need to
    setRefCell(
        V,
        mesh.solutionDict().subDict("couplingControls"),
        VRefCell,
        VRefValue
    );


//- U is the velocity field
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//Creates and initialises the relative face-flux field phi.
#include "createPhi.H"

//- Create charge density field
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(
        "rho",
        dimensionSet(0,-3,1,0,0,1,0),
        scalar(0)
    )
);


Info<< "Reading electrokineticProperties\n" << endl;

IOdictionary electrokineticProperties
(
    IOobject
    (
        "electrokineticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const bool electroNeutrality(
    electrokineticProperties.lookupOrDefault<bool>("electroNeutrality",false) //readBool
);

Info<< "Electro-neutrality flag set to "<< electroNeutrality<<endl;

//- Read fluid relative dielectric
const dimensionedScalar epsilonR_
(
    electrokineticProperties.lookup("epsilonR")
);

//- Calculate effective dielectric
const dimensionedScalar epsilon_( epsilonR_ * epsilonZero_);

//- Read chemical species
PtrList<dictionary> species;

PtrList<entry> specEntries
(
    electrokineticProperties.lookup("species")
);

species.setSize(specEntries.size());

//- Create concentration fields
PtrList<volScalarField> cPtrL(species.size());

forAll(species,sp)
{

    cPtrL.set(
        sp,
        new volScalarField
        (
            IOobject
            (
                specEntries[sp].keyword(),
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )

    );

    species.set(sp,new dictionary(specEntries[sp].dict()));

}

//Initialize list of Courant numbers for all species to 0
List<scalar> speciesCoNum(specEntries.size(),0.0);
List<scalar> speciesMeanCoNum(specEntries.size(),0.0);

//Initialize max Co # of all species to 0
scalar CoNum = 0.0;
