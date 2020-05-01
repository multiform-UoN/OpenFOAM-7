{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pisoDict = coupling.dict();

    scalar maxCo
    (
        pisoDict.lookupOrDefault<scalar>("maxCo", 0.8)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pisoDict.lookupOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pisoDict.lookupOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pisoDict.lookupOrDefault<scalar>("maxDeltaT", great)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    // Set the reciprocal time-step from the local Courant number
    volScalarField& C_ini(cPtrL[0]);
    dimensionedScalar D_ini(species[0].lookup("D"));
    dimensionedScalar Z_ini(species[0].lookup("Z"));


    //Set reciprocal time-step from the local Courant number of first species
    rDeltaT.ref() = max
    (
        1/dimensionedScalar(dimTime, maxDeltaT),
        fvc::surfaceSum(mag(phi + phiNP*Z_ini*D_ini))()()
       /((2*maxCo)*mesh.V())
    );

    //Go through all other species and check to see if rDeltaT needs updating
    forAll(species,sp){

        //Gather pointers to specific ion species sp
        volScalarField&   C(cPtrL[sp]);
        dimensionedScalar D(species[sp].lookup("D"));
        dimensionedScalar Z(species[sp].lookup("Z"));

        rDeltaT.ref() = max
        (
            rDeltaT.ref(),
            fvc::surfaceSum(mag(phi + phiNP*Z*D))()()
            /((2*maxCo)*mesh.V())
        );
    }

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }
}
