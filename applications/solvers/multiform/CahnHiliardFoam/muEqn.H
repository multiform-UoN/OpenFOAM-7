// -- DEBUG
if (debugCH)
{
  volScalarField laplmu
  (
    "laplmu",
    (thetaMu-scalar(1))*M*fvc::laplacian(mu)
  );
  laplmu.write();
  volScalarField ddtAlpha
  (
    "ddtAlpha",
    thetaMu*(alpha-alpha.prevIter())/runTime.deltaT()
  );
  ddtAlpha.write();
  mu.write();
  alpha.write();
  runTime++;
}
// -- DEBUG

while (pimple.correctNonOrthogonal())
{

  const dimensionedScalar muCorr
  (
    fvc::domainIntegrate(mu-pot)/vol
  );

  fvScalarMatrix muEqn
  (
    fvc::ddt(alpha)/thetaMu
    +
    fvm::Sp(pseudoTimeStepMu,mu)
    -
    pseudoTimeStepMu*mu.prevIter()
    -
    fvm::laplacian(M, mu)
    ==
    thetaMu0*M*fvc::laplacian(mu)//.prevIter()) // use previous time?
    // This term is for solvability
    -
    fvc::Sp(lambdaMu*mag(muCorr),mu)
    +
    lambdaMu*(mag(muCorr)-muCorr)*mu
  );

  // muEqn.mySetReference(0,scalar(0),true); // no longer needed
  muEqn.relax();
  muEqn.solve();

  // // mu defined up to a constant
  // // current implementation of setReference is not appropriate
  // // for the moment just rescaling to have zero mean
  // mu -= fvc::domainIntegrate(mu-pot)/vol;
  // mu.correctBoundaryConditions();
  // // no longer needed

}
