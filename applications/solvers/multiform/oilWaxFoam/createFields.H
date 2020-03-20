
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

    // concentration diffusivity
    dimensionedScalar Dc
    (
        transportProperties.lookup("Dc")
    );

    // wax diffusivity
    dimensionedScalar Dpsi
    (
        transportProperties.lookup("Dpsi")
    );

    // temperature diffusivity
    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    // reference viscosity
    dimensionedScalar nu0
    (
        transportProperties.lookup("nu")
    );

    // critical wax concentration
    dimensionedScalar psim
    (
        transportProperties.lookup("psim")
    );

    // relaxation time oversaturation
     dimensionedScalar tau1
    (
        transportProperties.lookup("tau1")
    );

    // relaxation time undersaturation
     dimensionedScalar tau2
    (
        transportProperties.lookup("tau2")
    );

    // Stokes number
    dimensionedScalar St
    (
        transportProperties.lookup("St")
    );

    // reference temperature
    dimensionedScalar Trif
    (
        transportProperties.lookup("Trif")
    );

    // reference concentration
    dimensionedScalar crif
    (
        transportProperties.lookup("crif")
    );

    // concentration per unit temperature
    dimensionedScalar cT
    (
        transportProperties.lookup("cT")
    );

    // scaling factor for source terms
    dimensionedScalar psifact
    (
        transportProperties.lookup("psifact")
    );

    // viscosity destructuring term
    dimensionedScalar alfa
    (
        transportProperties.lookup("alfa")
    );

    // viscosity structuring/hardening/ageing coefficient
    dimensionedScalar teta
    (
        transportProperties.lookup("teta")
    );

    // Identity tensor
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
    volTensorField Gr(fvc::grad(U));

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

    Info<< "Reading field lamda\n" << endl;
    volScalarField lamda
    (
        IOobject
        (
            "lamda",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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
        nu0 * (scalar(1)+lamda*psi)
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

    //volVectorField Uslip("Uslip",U*0);
    surfaceScalarField phislip("phislip",phi*0);


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PIMPLE"), pRefCell, pRefValue);
    mesh.setFluxRequired(p.name());

    #include "createMRF.H"       //added
    #include "createFvOptions.H" //added

    volScalarField invSt("invSt",nu/nu/St);

/////////////////////////////////////////////////////////////////
