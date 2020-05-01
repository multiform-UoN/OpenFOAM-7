//- The new Eigenvalues are calculated using the Rayleight quotient
//  such that lambda corresponds to the spectral radius of the ADE operator
dimensionedScalar lambdaRayleigh
(
    one
    //lambda
    *
    fvc::domainIntegrate(psi*psi.prevIter())
    /
    fvc::domainIntegrate(psi*psi)
);

lambda = shift*lambda + theta*lambdaRayleigh;
lambdaAdj = lambda;

//- Calculate new residual
lambdaRes = mag(lambdaOld - lambda.value());

Info<< "Eigenvalue : " << lambda.value() << endl;
Info<< "Residual on eigenvalue : " << lambdaRes  << endl;
