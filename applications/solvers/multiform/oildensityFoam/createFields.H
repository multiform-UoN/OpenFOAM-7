
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
    // reference viscosity
    dimensionedScalar nu0
    (
        transportProperties.lookup("nu")
    );
    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")
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
    // reference temperature viscosity
    dimensionedScalar Ti
    (
        transportProperties.lookup("Ti")
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
    //Molecular mass
    dimensionedScalar MM
    (
        transportProperties.lookup("MM")
    );
    //constant for solve density
    dimensionedScalar as
    (
        transportProperties.lookup("as")
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
    //specific heat
    dimensionedScalar cp
    (
        transportProperties.lookup("cp")
    );
    //thermal conductivity
    dimensionedScalar k
    (
        transportProperties.lookup("k")
    );
    // Identity tensor
    dimensionedTensor I
      (
      "I",
      dimensionSet(0, 0, 0, 0, 0, 0, 0),
      tensor(1,0,0,0,1,0,0,0,1)
      );

    //surfaceScalarField psifact("psifact",ni/MM);
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

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
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
        turbulence->nu() * (scalar(1)+lamda*psi)
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho0-as*(T - Ti)
    );

//
    #include "createMRF.H"       //added
    #include "createFvOptions.H" //added

    volScalarField invSt("invSt",nu/nu/St);

/////////////////////////////////////////////////////////////////
