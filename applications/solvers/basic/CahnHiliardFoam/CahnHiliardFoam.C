/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\    /   O peration     | Website:  https://openfoam.org
\\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
\\/     M anipulation  |
-------------------------------------------------------------------------------
License
This file is part of OpenFOAM.

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
CahnHiliardFoam

Description
Solves Cahn Hilliard Model

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  #include "setRootCaseLists.H"

  #include "createTime.H"
  #include "createMesh.H"

  simpleControl simple(mesh);

  #include "createFields.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nCalculating Cahn-Hiliard solution\n" << endl;

  label muRefCell = 0;
  scalar muRefValue = 0.0;

  #include "updatePot.H"
  #include "updateMu.H"
  #include "computeEnergy.H"

  while (simple.loop(runTime))
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    const volScalarField mu0(mu);

    while (simple.correctNonOrthogonal())
    {

      for(int i=1; i<10; i++)
      {
      fvScalarMatrix alphaEqn
      (
        - fvm::laplacian(epsilon*epsilon, alpha)
        ==
        - fvm::Sp(pot_imp,alpha) - (pot-pot_imp*alpha) + mu
      );

      alphaEqn.solve();

      #include "updatePot.H"

      #include "computeEnergy.H"

      }


      for(int i=1; i<10; i++)
      {
      fvScalarMatrix muEqn
      (
        fvc::ddt(alpha) - fvm::laplacian(M/scalar(2), mu)
        ==
        //fvOptions(mu)
        M/scalar(2)*fvc::laplacian(mu0)
      );

      //fvOptions.constrain(muEqn);
      //muEqn.setReference(muRefCell, muRefValue);
      muEqn.solve();
      //fvOptions.correct(mu);
      }
      #include "computeEnergy.H"

    }

    runTime.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
    << nl << endl;
  }

  Info<< "End\n" << endl;

  return 0;
}


// ************************************************************************* //
