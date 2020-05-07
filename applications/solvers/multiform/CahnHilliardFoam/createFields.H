//- Read transport properties

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

dictionary freeEnergy(transportProperties.subDict("freeEnergy"));

Info<< "Reading free energy constants \n" << endl;

const dimensionedScalar fea
(
    freeEnergy.lookup("a")
);

const dimensionedScalar feb
(
    freeEnergy.lookup("b")
);

Info<< "Reading diffuse interface length epsilon\n" << endl;

const dimensionedScalar epsilon
(
    transportProperties.lookup("epsilon")
);

const dimensionedScalar epsSq(pow(epsilon,2));

//- Reading fields

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

//- Create dfdc
volScalarField fPrime(fea*pow(alpha,3.) - feb*alpha);

Info<< "Creating field mu (chemical potential) \n" << endl;

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

mu = fPrime - epsSq*fvc::laplacian(alpha);
mu.correctBoundaryConditions();
