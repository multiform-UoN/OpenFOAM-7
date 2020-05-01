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
    EqRetCellFoam

Description
    Compute upscaled dynamicRetention model given an isotherm

Developers
Matteo Icardi, Nottingham (2019)
Federico Municchi, Nottingham (2019)


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"
#include "dynamicRetentionFvPatchScalarField.H"
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


    while (runTime.loop())
    {
        Info << "\n\nTime = " << runTime.timeName() << nl << endl;

        const dimensionedScalar normT0
        (
            (fvc::domainIntegrate(T0)/totalVol)
        );

        fluxPatch.SSp( - dSeqdC / lambda1 );
        fluxPatch.SSu( (normS) / lambda1 );

        for(label corr(0);corr<nNonOrthCorrCell;corr++)
        {
            fvScalarMatrix TEqn
            (
                fvm::ddt(T)
              + fvm::div(phi,T)
              - fvm::laplacian(DT, T)
              ==
              (
              -  normS
              +  dSeqdC * T0
              )
              /lambda1
              +
              A
            );

            TEqn.relax();


            if(corr==0)
            {
                residualT = TEqn.solve().max().initialResidual();
            }
            else
            {
                TEqn.solve();
            }

        }

        T += mag(T);
        T /= scalar(2);

        const dimensionedScalar normT
        (
            (fvc::domainIntegrate(T)/totalVol)
        );


        dimensionedScalar lambda10 = lambda1;

        dimensionedScalar rescaling( normT/normT0 );

        if (  (-  normS +  dSeqdC * normT0).value() <0 )
        {
          rescaling = normT0/normT;
        }
        lambda1 = lambda1 * rescaling;

        Info << "Lambda " << scalar(1)/lambda1.value() << endl;

        // Rescaling to preserve the initial average C and S
        T *=  normT0/normT;
        T.correctBoundaryConditions();
        fluxPatch.rescaleS(normS);
        T0 = T;

        lambdaRes = mag((lambda10-lambda1)/lambda10);

        Info << "Residual on T "    <<  residualT.value() << endl;
        Info << "Residual on Lambda "    <<  lambdaRes.value() << endl;

        if (lambdaRes.value() < powerIterationTolerance)
        {
          Info << "Converged!\n";
          break;
        }

    }

    Info << "End\n";

    return 0;
}


// ************************************************************************* //
