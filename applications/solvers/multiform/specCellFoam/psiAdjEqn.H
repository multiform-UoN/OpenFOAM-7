
    //- Equation for the adjoint eigenfunction psi
    fvScalarMatrix psiAdjEqn
    (
        fvm::div(phiAdj,psiAdj)
      - fvm::laplacian(D, psiAdj)
    );

    psiAdjEqn.relax();
    solve(psiAdjEqn == lambdaAdj*psiAdj.prevIter());
