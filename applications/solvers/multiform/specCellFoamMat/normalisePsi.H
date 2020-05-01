{

    psi = mag(psi) *
    (
        Foam::sqrt
        (
          fvc::domainIntegrate(psi.prevIter()*psi.prevIter())
          /
          fvc::domainIntegrate(psi*psi)
        )
    );

    psi.correctBoundaryConditions();

    volScalarField psiFieldRes
    (
        magSqr( psi - psi.prevIter() )
    );

    psiRes = Foam::sqrt
                (
                  fvc::domainIntegrate(psiFieldRes)
                ).value();
}

Info<<"Residual on psi : " << psiRes << endl;
