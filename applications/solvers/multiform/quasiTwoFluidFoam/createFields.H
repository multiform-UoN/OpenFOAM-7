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

    const dictionary rheology(transportProperties.subDict("rheology"));

    dimensionedScalar Dc
    (
        transportProperties.lookup("Dc")
    );

    dimensionedScalar nu0
    (
        rheology.lookup("nu")
    );

    dimensionedScalar cm
    (
        rheology.lookup("cm")
    );

    dimensionedScalar Ca
    (
        rheology.lookup("Ca")
    );

    dimensionedScalar Cb
    (
        rheology.lookup("Cb")
    );

    dimensionedScalar Cc
    (
        rheology.lookup("Cc")
    );

    // Stokes number
    dimensionedScalar St
    (
        transportProperties.lookup("St")
    );

    dimensionedTensor I
    (
    "I",
    dimensionSet(0, 0, 0, 0, 0, 0, 0),
    tensor(1,0,0,0,1,0,0,0,1)
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
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

    Info<< "Reading field c\n" << endl;
    volScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        nu0 * pow( (1.0- (c/cm)), -2.0)
    );

    #include "createPhi.H"

    volVectorField Uslip
    (
        IOobject
        (
            "Uslip",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

//    volVectorField Uslip("Uslip",U*0);
    surfaceScalarField phislip("phislip",phi*0);


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    #include "createMRF.H"       //added
    #include "createFvOptions.H" //added

    const volScalarField invSt("invSt",nu/nu/St);

/////////////////////////////////////////////////////////////////
