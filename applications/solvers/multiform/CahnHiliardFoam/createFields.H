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


Info<< "Reading free energy constant W0\n" << endl;

dimensionedScalar W0
(
    transportProperties.lookup("W0")
);

Info<< "Reading diffuse interface length epsilon\n" << endl;

dimensionedScalar epsilon
(
    transportProperties.lookup("epsilon")
);

#include "createFvOptions.H"


Info<< "\nCreating field mu (chemical potential)" << endl;

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
    W0 * alpha * ( sqr(alpha)-scalar(1.) )
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
    W0 * ( scalar(3.)*sqr(alpha)-scalar(1.) )
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
    //pot - epsilon*epsilon * fvc::laplacian(alpha)
);

const dimensionedScalar pi("pi",dimless/dimLength,M_PI);
