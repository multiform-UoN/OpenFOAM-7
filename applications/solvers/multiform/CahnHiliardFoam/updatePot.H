// Double-well potential
pot = W0 * alpha * ( sqr(alpha)-scalar(1.) ) ;
potImp = W0 * ( scalar(3.)*sqr(alpha)-scalar(1.) );

// solvability condition
// pot -= fvc::domainIntegrate(pot-mu)/vol;
