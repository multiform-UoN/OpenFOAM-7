{

    X -= fvc::domainIntegrate(X)/vol;
    X.correctBoundaryConditions();
}
