{

    psiAdj = mag(psiAdj) *
    (
        Foam::sqrt
        (
          fvc::domainIntegrate(psiAdj.prevIter()*psiAdj.prevIter())
          /
          fvc::domainIntegrate(psiAdj*psiAdj)
        )
    );

    psiAdj.correctBoundaryConditions();

    volScalarField psiAdjFieldRes
    (
        magSqr( psiAdj - psiAdj.prevIter() )
    );

    psiAdjRes = Foam::sqrt
                (
                  fvc::domainIntegrate(psiAdjFieldRes)
                ).value();

}

Info<<"Residual on psiAdj : " << psiAdjRes << endl;
