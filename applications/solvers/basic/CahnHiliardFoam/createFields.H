Info<< "Reading field alpha (phase field) \n" << endl;

volScalarField alpha
(
    IOobject
    (
        "alpha",
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


Info<< "Reading mobility M\n" << endl;

dimensionedScalar M
(
    transportProperties.lookup("M")
);

Info<< "Reading diffuse interface length epsilon\n" << endl;

dimensionedScalar epsilon
(
    transportProperties.lookup("epsilon")
);

#include "createFvOptions.H"


Info<< "Creating field mu (chemical potential) \n" << endl;

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
    scalar(100) * alpha * (alpha*alpha-scalar(1) )
);

volScalarField pot_imp
(
    IOobject
    (
        "pot",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    scalar(100) * (scalar(3)*alpha*alpha-scalar(1) )
);

volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
