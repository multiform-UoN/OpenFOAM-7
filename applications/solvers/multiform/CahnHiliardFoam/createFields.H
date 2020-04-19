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

const dimensionedScalar M
(
    transportProperties.lookup("M")
);

Info<< "Reading free energy constant W0\n" << endl;

const dimensionedScalar W0
(
    transportProperties.lookup("W0")
);

Info<< "Reading diffuse interface length epsilon\n" << endl;

const dimensionedScalar epsilon
(
    transportProperties.lookup("epsilon")
);

const dimensionedScalar epsSq(pow(epsilon,2));

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

const scalar thetaMu
(
    readScalar(pimple.dict().lookup("muRelax"))
);

const scalar thetaAlpha
(
    readScalar(pimple.dict().lookup("alphaRelax"))
);

const scalar thetaMu0
(
    (scalar(1)-thetaMu)/thetaMu
);

const scalar thetaAlpha0
(
    (scalar(1)-thetaAlpha)/thetaAlpha
);

const dimensionedScalar lambdaMu
(
    "lambdaMu", dimless/dimTime, readScalar(pimple.dict().lookup("lambdaMu"))
);

const dimensionedScalar lambdaAlpha
(
    "lambdaAlpha", dimless/(dimLength*dimLength), readScalar(pimple.dict().lookup("lambdaAlpha"))
);

const scalar alphaRelaxTime
(
  readScalar(pimple.dict().lookup("alphaRelaxTime"))
);

const scalar muRelaxTime
(
  readScalar(pimple.dict().lookup("muRelaxTime"))
);

const dimensionedScalar pseudoTimeStepAlpha
(
    "pseudoTimeStepAlpha",
    dimless/(dimLength*dimLength),
    (scalar(1)-alphaRelaxTime)
    /
    (runTime.deltaT().value()*alphaRelaxTime)
);

const dimensionedScalar pseudoTimeStepMu
(
    (scalar(1)-muRelaxTime)
    /
    (runTime.deltaT()*muRelaxTime)
);

const bool debugCH
(
    pimple.dict().lookupOrDefault<bool>("debugCH",false)
);

const dimensionedScalar thetaAlphaEps(scalar(1)/(epsSq*thetaAlpha));

dimensionedScalar meanpot(fvc::domainIntegrate(pot)/vol);
