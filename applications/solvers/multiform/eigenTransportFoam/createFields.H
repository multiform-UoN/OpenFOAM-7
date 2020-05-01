//- Get specific controls for power method

const dictionary& pwrDict(pwrctrl.dict());

scalar eigenTol(readScalar(pwrDict.lookup("eigenTol")));
Switch sameEigen(pwrDict.lookup("sameEigenvalue"));

Info<< "Reading field psi\n" << endl;

volScalarField psi
(
    IOobject
    (
        "psi",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field psiAdj\n" << endl;

volScalarField psiAdj
(
    IOobject
    (
        "psiAdj",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


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

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity Gamma\n" << endl;

dimensionedScalar Gamma
(
    transportProperties.lookup("Gamma")
);

//- Eigenvalue
dimensionedScalar lambda("lambda", dimless/dimTime, scalar(1));
dimensionedScalar lambdaAdj("lambda", dimless/dimTime, scalar(1));
dimensionedScalar vol("vol",dimVolume, gSum(mesh.V()));

//- Create adjoint flux
surfaceScalarField phiAdj("phiAdj",-phi);



