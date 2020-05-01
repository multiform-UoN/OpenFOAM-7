{
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix hEqn
        (
            fvm::ddt(T)
          - fvm::laplacian(DTs, T)
            ==
            fvOptions(T)
        );

        hEqn.relax();

        fvOptions.constrain(hEqn);

        hEqn.solve();

        fvOptions.correct(T);
    }
}

