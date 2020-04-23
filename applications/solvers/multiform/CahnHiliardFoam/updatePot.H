{

    //const volScalarField alpha1(alpha*scalar(2)-scalar(1));

    // Double-well potential
    pot = W0 * alpha * ( sqr(alpha)-scalar(1) ) ;
    // potImp = pos(W0 * ( scalar(3)*sqr(alpha)-scalar(1) ));
    potImp = W0;

    // solvability condition
    // pot -= fvc::domainIntegrate(pot-mu)/vol;

    meanpot = (fvc::domainIntegrate(pot)/vol);

}
