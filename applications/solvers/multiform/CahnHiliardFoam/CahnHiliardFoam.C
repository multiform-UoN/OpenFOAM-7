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
Solves Cahn-Hilliard equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
  #include "setRootCaseLists.H"

  #include "createTime.H"
  #include "createMesh.H"

  pimpleControl pimple(mesh);

  #include "createFields.H"


  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nCalculating Cahn-Hiliard solution\n" << endl;

  Info << "Initial Energy calculation" << endl;

  #include "computeEnergy.H"
  #include "updateMu.H"

  Info << "Start time loop" << endl;


  while (runTime.loop())
  {
    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Shall we update the non-linearity here?
    #include "updatePot.H"
    // Or just the solvability condition
    pot -= fvc::domainIntegrate(pot-mu)/vol;


    // -- DEBUG
    if (debugCH)
    {
      mu.write();
      alpha.write();
    }
    // -- DEBUG

    mu.storePrevIter();
    alpha.storePrevIter();

    while (pimple.loop())
    {

      // Shall we update the non-linearity here?
      // #include "updatePot.H"
      // Or just the solvability condition
      // pot -= fvc::domainIntegrate(pot-mu)/vol;

      // -- DEBUG
      if (debugCH)
      {
        runTime++;
        volScalarField laplalpha
        (
          "laplalpha",
          epsSq*(thetaAlpha-scalar(1))*fvc::laplacian(alpha)
        );
        laplalpha.write();

        Info << fvc::domainIntegrate(mu-pot) << endl;
        pot.write();
        volScalarField mupot("mupot",mu-pot);
        mupot.write();
      }
      // -- DEBUG

      while (pimple.correctNonOrthogonal())
      {

        fvScalarMatrix alphaEqn
        (
          - fvm::laplacian(epsSq, alpha)
          // + fvm::Sp(mag(potImp)+scalar(1),alpha) // semi-implicit potential
          // + fvm::Sp(scalar(1)/thetaAlpha,alpha) // add to diagonal
          ==
          (thetaAlpha-scalar(1))*epsSq*fvc::laplacian(alpha.prevIter()) // use previous time?
          // + (mag(potImp)+scalar(1))*alpha // semi-implicit potential
          // + (scalar(1)/thetaAlpha)*alpha // remove from diagonal
          - thetaAlpha*pot
          + thetaAlpha*mu.prevIter() // use previous time?
          // This term is to try to conserve mass
          - fvc::domainIntegrate(alpha-alpha.prevIter())/vol
        );

        alphaEqn.relax();
        alphaEqn.solve();

        // Shall we update the non-linearity here?
        // #include "updatePot.H"
        // Or just the solvability condition
        // pot -= fvc::domainIntegrate(pot-mu)/vol;


        // -- DEBUG
        if (debugCH)
        {
          Info << "Integral alpha " << fvc::domainIntegrate(alpha) << endl;
        }
        // -- DEBUG

      }

      // -- DEBUG
      if (debugCH)
      {
        volScalarField laplmu
        (
          "laplmu",
          (thetaMu-scalar(1))*M*fvc::laplacian(mu)
        );
        laplmu.write();
        volScalarField ddtAlpha
        (
          "ddtAlpha",
          thetaMu*(alpha-alpha.prevIter())/runTime.deltaT()
        );
        ddtAlpha.write();
        mu.write();
        alpha.write();
        runTime++;
      }
      // -- DEBUG

      while (pimple.correctNonOrthogonal())
      {
        fvScalarMatrix muEqn
        (
          thetaMu*fvc::ddt(alpha)
          -
          fvm::laplacian(M, mu)
          ==
          (thetaMu-scalar(1))*M*fvc::laplacian(mu.prevIter()) // use previous time?
        );

        // muEqn.mySetReference(0,scalar(0),true);
        muEqn.relax();
        muEqn.solve();

        // mu defined up to a constant
        // current implementation of setReference is not appropriate
        // for the moment just rescaling to have zero mean
        mu -= fvc::domainIntegrate(mu)/vol;
        mu.correctBoundaryConditions();

        // -- DEBUG
        if (debugCH)
        {
          Info << "Integral mu " << fvc::domainIntegrate(mu) << endl;
        }
        // -- DEBUG
      }

      // -- DEBUG
      if (debugCH)
      {
        mu.write();
        alpha.write();
      }
      // -- DEBUG


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
