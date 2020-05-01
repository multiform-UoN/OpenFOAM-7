{
  label pnpiter=1;
  // --- Coupling loop
      while (pimple.correct())
      {
          Info << "Fluid PNP iteration " << pnpiter << endl;
          pnpiter++;

            //-Update charge density
        
            rho *= scalar(0); //- Simply reset to zero

            if (!electroNeutrality){

                forAll(species,sp)
                {
                    dimensionedScalar Z(species[sp].lookup("Z"));
                    rho += Fy_*cPtrL[sp]*Z;
                }

            }
            

            //- Solve Poisson
            while (pimple.correctNonOrthogonal())
            {
                V.storePrevIter();  // WHY WE NEED THIS?
                fvScalarMatrix VEqn
                (
                  - fvm::laplacian(epsilon_,V)
                  ==
                    rho
                );

                // Set the reference value for V to RefCell & RefValue, do nothing if not
                VEqn.setReference(VRefCell, VRefValue);

                VEqn.relax();
                VEqn.solve();
                V.relax();
            }

            //- Update Nerst-Planck flux
            phiNP = -fvc::flux(e_*fvc::grad(V)/(k_*T))*(!diffusionDominant);

            //- Solve species
            forAll(species,sp)
            {

                //- Collect pointers and coefficients
                volScalarField&   C(cPtrL[sp]);
                dimensionedScalar D(species[sp].lookup("D"));
                dimensionedScalar Z(species[sp].lookup("Z"));

                //- Evaluate specific Nerst-Plank flux
                phiNP *= (Z*D);

                scalar CoNum = 0.0;
                scalar meanCoNum = 0.0;

                {
                    scalarField sumPhi
                    (
                        fvc::surfaceSum(mag(phi+phiNP))().primitiveField()
                    );

                    CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

                    meanCoNum =
                        0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
                }

                Info<< "Courant Number mean: " << meanCoNum
                    << " max: " << CoNum << endl;

                //Check if electro-neutrality is wanted
                if (electroNeutrality && sp == species.size()-1){


                    //Look up the first species
                    const dimensionedScalar& Z_ini(species[0].lookup("Z"));
                    const volScalarField& C_ini(cPtrL[0]);

                    //Set concentration C to first species
                    C  = -Z_ini*C_ini/Z;

                    //Loop through all species except final species
                    for (int i=1;i <= species.size()-2;i++)
                    {
                        //Look up the next species
                        const dimensionedScalar& Z_iter(species[i].lookup("Z"));
                        const volScalarField& C_iter(cPtrL[i]);

                        //Add next species to C
                        C -=  Z_iter*C_iter/Z;

                    }

                    C.correctBoundaryConditions();

                }else{

                    //- Non-orthogonal correction loop
                    while (pimple.correctNonOrthogonal())
                    {

                        C.storePrevIter();
                        fvScalarMatrix CEqn
                        (
                            fvm::ddt(C)
                          + fvm::div(phi,C,"div(phi,C)")
                          + fvm::div(phiNP,C,"div(phiNP,C)")
                          - fvm::laplacian(D,C,"laplacian(D,C)")
                        );

                        //constrainFluxes(CEqn);

                        CEqn.relax();
                        CEqn.solve();

                        //- Regularise solution
                        {
                            //- Force solution lower bound
                            C = (mag(C) + C)/scalar(2.0); // ?? Unclear on what this is -R
                            C.correctBoundaryConditions();

                            //- Relax
                            C.relax();
                        }
                    }
                    
                }


                //- Rescale Nerst-Planck flux
                phiNP /= (Z*D);
            }
        }
        Info<< "For fluid region " << fluidRegions[i].name()<<":\n";
        forAll(species,sp)
        {
      		  Info << "Average concentration specie " << sp << " = "
               	 << Foam::gSum(cPtrL[sp]().field()*mesh.V())/Foam::gSum(mesh.V())
			   	       << endl;
//            cPtrL[sp]().write();
        }

}
