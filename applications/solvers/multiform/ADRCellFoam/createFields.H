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

    volScalarField Tadj
    (
        IOobject
        (
            "Tadj",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        T
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
            IOobject::NO_WRITE
        ),
        mesh
    );

    volVectorField Uplus
    (
        IOobject
        (
            "Uplus",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        U
    );

    Info<< "Reading field Xi\n" << endl;

    volVectorField Xi
    (
        IOobject
        (
            "Xi",
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

    dimensionedScalar R
    (
        transportProperties.lookup("R")
    );


    const dimensionedScalar totalVol("totalVol", dimVolume, gSum(mesh.V()));

    const dimensionedScalar unitLength("unitLength",dimLength,1);

    const dimensionedVector Uav = fvc::domainIntegrate(U) / totalVol;

    const dimensionedTensor I
    (
        "I",
        dimensionSet(0, 0, 0, 0, 0, 0, 0),
        tensor(1, 0, 0, 0, 1, 0, 0, 0, 1)
    );

    //- Read solution controls
    dictionary spectralCellControls
    (
        mesh.solutionDict().subDict("spectralCellControls")
    );

    scalar powerIterationTolerance
    (
        readScalar(spectralCellControls.lookup("powerIterationTolerance"))
    );

    label maxPowerIter
    (
        readLabel(spectralCellControls.lookup("powerIterationMaxIter"))
    );

    label nNonOrthCorrCell
    (
        readLabel(spectralCellControls.lookup("nNonOrthogonalCorrectors"))
    );

    dictionary dispersionControls
    (
        mesh.solutionDict().subDict("dispersionControls")
    );

    label nNonOrthCorr
    (
        readLabel(dispersionControls.lookup("nNonOrthogonalCorrectors"))
    );

    scalar dispersionTolerance
    (
        readScalar(dispersionControls.lookup("dispersionTolerance"))
    );

    label maxDispersionIter
    (
        readLabel(dispersionControls.lookup("dispersionMaxIter"))
    );
