//- Calculate effective diffusion coefficient
dimensionedTensor Deff
(
    "Deff",
    fvc::domainIntegrate
    (
        (
             beta*D*(tensor::I + fvc::grad(X))
             + ( (beta*Ueff-Ustar)*X )
         )
    ) / vol
);

//- Calculate first order corrector
dimensionedVector CI
(
    "CI",
    fvc::domainIntegrate(beta*X) / vol
);

//- Calculate alternative diffusion coefficient
dimensionedTensor DeffII
(
    "DeffII",
    fvc::domainIntegrate
    (
        beta*D
        * (tensor::I + fvc::grad(X))
        &
        (tensor::I + fvc::grad(X))().T()
    )  / vol
);
