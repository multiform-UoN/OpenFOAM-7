const dimensionedScalar vol("vol",pow(dimLength,3),gSum(mesh.V()));

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

volTensorField DT
(
    IOobject
    (
        "DT",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


volScalarField source
(
    IOobject
    (
        "source",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

source -= fvc::domainIntegrate(source)/vol; // force zero integral


Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


Info<< "Reading average value\n" << endl;

dimensionedScalar avgT
(
    "avgT",T.dimensions(),readScalar(transportProperties.lookup("average"))
);


dimensionedScalar rate
(
    "lagrangeCoeff",dimless/dimTime/T.dimensions(),readScalar(transportProperties.lookup("lagrangeCoeff"))
);


// volScalarField coeff
// (
//     IOobject
//     (
//         "coeff",
//         runTime.timeName(),
//         mesh,
//         IOobject::READ_IF_PRESENT,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     rate*(avgT)
// );
