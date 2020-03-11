//Viscosity Evolution Equation
const volScalarField SlImp
  ("SlImp", alfa*gamma);


  while (pimple.correctNonOrthogonal())
  {
   fvScalarMatrix lamdaEqn
   (
        fvm::ddt(lamda)
        + fvm::div(phipsi,lamda)
        - fvm::laplacian(Dpsi,lamda)
        ==
        - fvm::Sp(psi*SlImp,lamda)
        + psi*pow(teta,-1)
    );
   lamdaEqn.relax();
   lamdaEqn.solve();

  }
