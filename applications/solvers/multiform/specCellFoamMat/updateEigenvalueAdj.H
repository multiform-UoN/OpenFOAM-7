//- The new Eigenvalues are calculated using the Rayleight quotient
//  such that lambda corresponds to the spectral radius of the ADE operator
dimensionedScalar lambdaAdjRayleigh
(
    one
    //lambdaAdj
    *
    fvc::domainIntegrate(psiAdj*psiAdj.prevIter())
    /
    fvc::domainIntegrate(psiAdj*psiAdj)
);

lambdaAdj = shift*lambdaAdj + theta*lambdaAdjRayleigh;
lambda = lambdaAdj;

//- Calculate new residual
lambdaRes = mag(lambdaOld - lambda.value());

Info<< "Eigenvalue : " << lambda.value() << endl;
Info<< "Residual on eigenvalue : " << lambdaRes  << endl;
