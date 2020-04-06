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

dimensionedScalar epsSq(pow(epsilon,2));

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

volScalarField potImp
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

const dimensionedScalar vol( "vol", dimVol, gSum( mesh.V() ) );
const dimensionedScalar one( "one", dimless/dimTime, scalar(1) );

Info<< "Reading relaxation coefficients\n" << endl;

scalar thetaMu
(
    readScalar(pimple.dict().lookup("muRelax"))
);

scalar thetaAlpha
(
    readScalar(pimple.dict().lookup("alphaRelax"))
);

bool debugCH
(
    pimple.dict().lookupOrDefault<bool>("debugCH",false)
);
