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

    // Linearised isotherm
    scalar dSeqdC
    (
        transportProperties.lookupOrDefault<scalar>("dSeqdC",scalar(1))
    );

    // Source term
    const dimensionedScalar A
    (
        transportProperties.lookup("A")
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
        spectralCellControls.lookupOrDefault<scalar>("powerIterationTolerance",1e-10)
    );

    label nNonOrthCorrCell
    (
        spectralCellControls.lookupOrDefault<label>("nNonOrthogonalCorrectors",2)
    );

    dimensionedScalar lambda1
    (
      "eigenvalue",
      dimTime,
      scalar(1)/spectralCellControls.lookupOrDefault<scalar>("lambda0",1)
    );

    volScalarField T0(T);

    label powerIter(0);

    dimensionedScalar residualT("residualT",dimless,scalar(0.));

    dimensionedScalar normT0("normT",T.dimensions(),scalar(0.));

    dimensionedScalar lambdaRes("lambdaRes",dimless,scalar(0.));

    label patchI = mesh.boundaryMesh().findPatchID("solid");

    if (!isA<dynamicRetentionFvPatchScalarField>(T.boundaryField()[patchI]))
    {
      FatalError << "Not a dynamicRetention patch!!" << endl;
    }

    dynamicRetentionFvPatchScalarField& fluxPatch
    (
        refCast<dynamicRetentionFvPatchScalarField>
        (
            T.boundaryFieldRef()[patchI]
        )
    );

    scalar normS(fluxPatch.Savg());
