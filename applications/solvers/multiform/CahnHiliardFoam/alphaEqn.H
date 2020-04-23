// -- DEBUG
if (debugCH)
{
  volScalarField laplalpha
  (
    "laplalpha",
    epsSq*(thetaAlpha-scalar(1))*fvc::laplacian(alpha)
  );
  laplalpha.write();

  pot.write();
  volScalarField mupot("mupot",mu-pot);
  mupot.write();
}
// -- DEBUG

//TODO: try potential formulation substituting c into the equation for mu and using linear part of df/dc
//      into the equation for c.

while (pimple.correctNonOrthogonal())
{

  // Shall we update the non-linearity here?
  // #include "updatePot.H"

  const dimensionedScalar massCor
  (
    fvc::domainIntegrate(alpha-alpha.prevIter())/vol
  );

  fvScalarMatrix alphaEqn
  (
    fvm::Sp(pseudoTimeStepAlpha,alpha)
    -
    pseudoTimeStepAlpha*alpha.prevIter()
    -
    fvm::laplacian(alpha)
    ==
    thetaAlpha0*fvc::laplacian(alpha)// .prevIter()) // use previous time?
    -
    thetaAlphaEps*(pot - mu) // .prevIter() // use previous time?
    // -- Linearisation of pot
    -
    fvm::Sp(thetaAlphaEps*potImp,alpha)
    +
    thetaAlphaEps*potImp*alpha.prevIter()
    // -- This term is to conserve mass
    -
    fvc::Sp(lambdaAlpha*mag(massCor),alpha)
    +
    lambdaAlpha*(mag(massCor)-massCor)*alpha
  );

  alphaEqn.relax();
  alphaEqn.solve();

}
