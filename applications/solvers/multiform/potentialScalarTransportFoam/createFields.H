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


Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    transportProperties.lookup("DT")
);

#include "createPhi.H"

#include "createFvOptions.H"

Info<< "Reading potential properties\n" << endl;

dimensionedScalar refDistance
(
    transportProperties.lookup("potRefD")
);

dimensionedScalar minDistance
(
    transportProperties.lookup("potMinD")
);

dimensionedScalar refSpeed
(
    transportProperties.lookup("potRefU")
);

const volScalarField wd(wallDist(mesh).y()/refDistance);

volScalarField pot
(
    IOobject
    (
        "pot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimless
);

pot = pow(wd+minDistance/refDistance,-1);

const volVectorField Up("Up",fvc::grad(pot*refDistance)*refSpeed);
Up.write();

surfaceScalarField phip
(
    IOobject
    (
        "phip",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::flux(Up)
);
