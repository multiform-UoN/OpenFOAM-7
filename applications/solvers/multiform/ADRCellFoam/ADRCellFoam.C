/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019
     \\/     M anipulation  | Matteo Icardi, Federico Municchi
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


Application
    ADRCellFoam

Description
    Compute apparent diffusivity, principal eigenvalue and first order correction
    for an advection diffusion reaction cell problem.
    Output creates the effectiveTransportProperties dictionary.

Developers
Matteo Icardi, Nottingham (2019)
Federico Municchi, Nottingham (2019)


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const surfaceScalarField phi("phi",fvc::flux(U));
    #include "CourantNo.H"

    dimensionedScalar lambda1("eigenvalue", dimTime, 1.0);
    dimensionedScalar lambda2 = lambda1;

    volScalarField T0adj(T);
    volScalarField T0(T);

    label powerIter(0);

    dimensionedScalar residualT("residualT",dimless,scalar(0.));
    dimensionedScalar residualTadj("residualTadj",dimless,scalar(0.));

    dimensionedScalar normT0("normT",T.dimensions(),scalar(0.));
    dimensionedScalar normT0adj("normTadj",Tadj.dimensions(),scalar(0.));

    dimensionedScalar lambdaRes("lambdaRes",dimless,scalar(0.));

    runTime++;
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Create adjoint flux
    surfaceScalarField phiAdj(-phi);
    

    do
    {
        Info<<endl<<"Spectral cell loop, iteration: " << powerIter <<endl;
        // Spectral cell problem

        for(label corr(0);corr<nNonOrthCorrCell;corr++)
        {
            fvScalarMatrix TEqn
            (
                fvm::div(phi,T)
              - fvm::laplacian(DT, T)
              ==
                fvm::Sp(R,T)
              + T0/lambda1
            );

            TEqn.relax();

            //Spectral cell problem Adjoint
            fvScalarMatrix TEqn2
            (
                  fvm::div(phiAdj,Tadj)
                - fvm::laplacian(DT, Tadj)
               ==
                  fvm::Sp(R,Tadj)
                + T0adj/lambda2
            );

            TEqn2.relax();

            if(corr==0)
            {
                residualT = TEqn.solve().max().initialResidual();
                residualTadj = TEqn2.solve().max().initialResidual();
            }
            else
            {
                TEqn.solve();
                TEqn2.solve();
            }

        }

        const dimensionedScalar normT0adj
        (
            sqrt(fvc::domainIntegrate(magSqr(T0adj))/totalVol)
        );
        const dimensionedScalar normT0
        (
            sqrt(fvc::domainIntegrate(magSqr(T0))/totalVol)
        );

        const dimensionedScalar normT
        (
            sqrt(fvc::domainIntegrate(magSqr(T))/totalVol)
        );

        const dimensionedScalar normTadj
        (
            sqrt(fvc::domainIntegrate(magSqr(Tadj))/totalVol)
        );

        dimensionedScalar lambda10 = lambda1;

        lambda1 = lambda1*
        (
            fvc::domainIntegrate(T*T0)
            /
            fvc::domainIntegrate(T0*T0)
        );

        Info << "Lambda " << 1.0/lambda1.value() << endl;

        lambda2 = lambda2*
        (
            fvc::domainIntegrate(Tadj*T0adj)
            /
            fvc::domainIntegrate(T0adj*T0adj)
        );

        Info << "Lambda_adj " << 1.0/lambda2.value() << endl;

        T *=  normT0/normT;
        Tadj *=  normT0adj/normTadj;

        T.correctBoundaryConditions();
        Tadj.correctBoundaryConditions();

        T0adj = Tadj;
        T0 = T;

        lambdaRes = mag(lambda10-lambda1)/lambda10;

        Info << "Residual on T "    <<  residualT.value() << endl;
        Info << "Residual on Tadj "    <<  residualTadj.value() << endl;
        Info << "Residual on Lambda "    <<  lambdaRes.value() << endl;


    } while (
        (
            lambdaRes.value() > powerIterationTolerance
            // ||
            // residualTadj.value() > powerIterationTolerance
        )
        &&
        (
            ++powerIter < maxPowerIter
        )
    );

    const dimensionedScalar normTTadj
    (
        (fvc::domainIntegrate(T*Tadj)/totalVol)
    );

    const volScalarField beta(T*Tadj);

    Uplus =
        (
            U + (DT*fvc::grad(Tadj)/Tadj)
            - (DT*fvc::grad(T)/T)
        );

    Uplus.correctBoundaryConditions();
 

    Info<< "\nCalculating closure variable distribution\n" << endl;


    const volScalarField DD
    (
        "DD",
        beta*DT 
    );
    const surfaceScalarField phiplus("phi",fvc::flux(beta*Uplus));

    const dimensionedVector Ueff ("Ueff", fvc::domainIntegrate(Uplus) / totalVol);

    label dispersionIter(0);

    dimensionedScalar initialResidual(0);
    //- Solve for dispersion
    do
    {

        Info<<endl<<"Dispersion loop, iteration: " << dispersionIter <<endl;

        for(label corr(0);corr<nNonOrthCorr;corr++)
        {
            fvVectorMatrix XiEqn
            (
                 fvm::div(phiplus, Xi)
               - fvm::laplacian(DD, Xi)
                 ==
                 beta*(Ueff-Uplus) + fvc::grad(DD)
            );

            XiEqn.relax();

            if(corr==0)
            {
                initialResidual = XiEqn.solve().max().initialResidual();
            }
            else
            {
                XiEqn.solve();
            }
        }

        Info << "Residual dispersion " << initialResidual.value() << endl;

    } while (
        (initialResidual.value() > dispersionTolerance)
        &&
        (++dispersionIter < maxDispersionIter)
    );

    Xi-=fvc::domainIntegrate(Xi)/totalVol;
    Xi.correctBoundaryConditions();
    runTime.writeAndEnd();


    volTensorField gradXi(fvc::grad(Xi));
    const dimensionedTensor DstdiffSymm
    (
        "DeffSymm" ,
        fvc::domainIntegrate(DD * (I+gradXi)&((I + gradXi)().T()))
         / totalVol
    );

    const dimensionedTensor Dstdiff
    (
        "Deff",
        fvc::domainIntegrate(DD * (I + gradXi) - beta*(U*Xi))
         / totalVol
    );

    const dimensionedVector firstCorrectCoeff
    (
        "CI",
        fvc::domainIntegrate(T*Xi)/totalVol
    );

    Info<< endl << "End\n" << "ExecutionTime = "
        << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    
   // Info<< "is div null? "<< fvc::div(phiplus);

    Info<<"**********************************************************"<<endl
        <<"                    EFFECTIVE PARAMETERS                "<<endl
        <<"\n   Dispersion tensor Symm : "
        << DstdiffSymm.value()<<endl
        <<"\n   Dispersion tensor Symm (trace) : "
        << tr(DstdiffSymm.value())<<endl
        <<"\n   Dispersion tensor : "
        << Dstdiff.value()<<endl
        <<"\n   Dispersion tensor  (trace) : "
        << tr(Dstdiff.value())<<endl
        <<"\n   Lambda " << scalar(1.0)/lambda1.value() << endl
        <<"\n   Lambda_adj " << scalar(1.0)/lambda2.value() << endl
        <<"\n   First order corrector coefficient : "
        << firstCorrectCoeff.value()
        <<"\n   Effective velocity : "
        << Ueff.value();
    
    //- Write dictionary for dispersion tensor
    Info<<"\nCreating effectiveTransportProperties dictionary..."<<endl;
    IOdictionary effectiveTransportProperties
    (
        IOobject
        (
            "effectiveTransportProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    effectiveTransportProperties.add
    (
        "Deff", 
        Dstdiff
    );
    
    effectiveTransportProperties.add
    (
        "DeffSymm", 
        DstdiffSymm
    );

    effectiveTransportProperties.add
    (
        "Eigenvalue",
        scalar(1.0)/lambda1.value() 
    );

    effectiveTransportProperties.add
    (
        "CI",
        firstCorrectCoeff
    );

    effectiveTransportProperties.add
    (
        "Ueff",
        Ueff
    );

    effectiveTransportProperties.regIOobject::write();

    return 0;
}


// ************************************************************************* //
