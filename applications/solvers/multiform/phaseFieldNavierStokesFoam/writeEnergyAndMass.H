{
    scalar vol(gSum(mesh.V()));
    Info<< "Avg Energy : "
        <<(
            fvc::domainIntegrate
            (
                sqr(sqr(2.0*alpha-1.0)-scalar(1.))/scalar(4.)
                + epsSq*magSqr(fvc::grad(alpha))*scalar(2.)
            )/vol
        ).value()
        << endl;

    Info<< "Avg mass : "
        << (fvc::domainIntegrate(alpha)/vol).value() << endl;
}
