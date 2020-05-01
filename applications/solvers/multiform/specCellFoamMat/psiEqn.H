
//- Equation for the eigenfunction psi
fvScalarMatrix psiEqn
(
  fvm::div(phi,psi)
  -
  fvm::laplacian(D, psi)
  ==
  fvm::Sp(shift*thetaShift*lambda,psi)
  +
  shift*(scalar(1)-thetaShift)*lambda*psi
);

//psiEqn.setReference(0,scalar(1));

psiEqn.relax();
// solve(psiEqn);
// solve(psiEqn == lambda*psi.prevIter());
solve(psiEqn == one*psi.prevIter());

// Info << "Norm psi "
// << (sqrt(fvc::domainIntegrate(magSqr(psi)))).value()
// << endl;
// Info << "Norm residual "
// << (sqrt(fvc::domainIntegrate(magSqr(psiEqn.A()*psi - psiEqn.H())))).value()
// << endl;
// Info << "Norm residual "
// << (sqrt(fvc::domainIntegrate(magSqr(fvc::div((U+Up)*psi)-fvc::laplacian(D*psi))))).value()
// << endl;
