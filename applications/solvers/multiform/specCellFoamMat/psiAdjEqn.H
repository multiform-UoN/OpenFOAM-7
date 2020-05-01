
//- Equation for the adjoint eigenfunction psi
fvScalarMatrix psiAdjEqn
(
  fvm::div(phiAdj,psiAdj)
  -
  fvm::laplacian(D, psiAdj)
  ==
  fvm::Sp(shift*thetaShift*lambdaAdj+fvc::div(phiAdj),psiAdj)
  +
  shift*(scalar(1)-thetaShift)*lambdaAdj*psiAdj
);

//psiAdjEqn.setReference(0,scalar(1));

psiAdjEqn.relax();
// solve(psiAdjEqn);
// solve(psiAdjEqn == lambdaAdj*psiAdj.prevIter());
solve(psiAdjEqn == one*psiAdj.prevIter());

// Info << "Norm psi "
// << (sqrt(fvc::domainIntegrate(magSqr(psiAdj)))).value()
// << endl;
// Info << "Norm residual "
// << (sqrt(fvc::domainIntegrate(magSqr(psiAdjEqn.A()*psi - psiAdjEqn.H())))).value()
// << endl;
// Info << "Norm residual "
// << (sqrt(fvc::domainIntegrate(magSqr(fvc::div((U+Up)*psiAdj)+fvc::laplacian(D*psiAdj)+fvc::div(phiAdj)*psiAdj)))).value()
// << endl;
