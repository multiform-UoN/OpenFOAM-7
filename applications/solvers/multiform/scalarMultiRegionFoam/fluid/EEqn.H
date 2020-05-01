{

#include "CourantNo.H"

while (pimple.correctNonOrthogonal())
  {

    fvScalarMatrix EEqn
    (
        fvm::ddt(T) + fvm::div(phi, T)
      - fvm::laplacian(DT, T)
     ==
       fvOptions(T)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(T);

  }
}
