    const fvMesh& mesh = fluidRegions[i];

    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];
    volScalarField& T = fieldsFluid[i];

    fv::options& fvOptions = fluidFvOptions[i];

 
    pimpleNoLoopControl& pimple = pimples.pimple(i);


