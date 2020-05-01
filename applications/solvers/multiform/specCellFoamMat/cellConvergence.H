{
    currCellIter++;

    volTensorField gradX(fvc::grad(X));

    //- Calculate effective diffusion coefficient
    dimensionedTensor Deff
    (
        "Deff",
        fvc::domainIntegrate
        (
            (
                 beta*D*(tensor::I + gradX)
                 + ( (Ueff-Ustar)*X )
             )
        ) / vol
    );

    //- Calculate alternative diffusion coefficient
    dimensionedTensor DeffII
    (
        "DeffII",
        fvc::domainIntegrate
        (
            beta*D
            * (tensor::I + gradX)
            &
            T(tensor::I + gradX)
        )  / vol
    );


    if(deltaDRes)
    {
        XRes = mag(tr(Deff).value() -tr(DeffII).value())
            /(mag(tr(Deff).value())+small);
    }
    else
    {
        XRes = fvc::domainIntegrate(magSqr(X - X.prevIter())).value()
            /(fvc::domainIntegrate(magSqr(X)).value());
    }

    cellConverged =
        (
            XRes < corrTol
        )
        ||
        (
            currCellIter > maxCellProbIter
        );

    Info<<"Residual : " << XRes << " Converged : " << cellConverged
        <<endl;

}
