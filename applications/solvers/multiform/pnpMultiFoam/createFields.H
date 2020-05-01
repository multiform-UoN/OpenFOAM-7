
Info<< "Reading electrokineticProperty files\n" << endl;

//Pointer lists for electroKineticProperty dictionaries
PtrList<IOdictionary> SolidElectroKineticProperties(solidRegions.size());
PtrList<IOdictionary> FluidElectroKineticProperties(fluidRegions.size());

//Pointer lists for relative and effective dielectrics per region
PtrList<dimensionedScalar> epsilonR_Solid(solidRegions.size());
PtrList<dimensionedScalar> epsilon_Solid(solidRegions.size());

PtrList<dimensionedScalar> epsilonR_Fluid(fluidRegions.size());
PtrList<dimensionedScalar> epsilon_Fluid(fluidRegions.size());

// Run through all solid regions
forAll(solidRegions, i){

	// Construct word for filename of region i
	word ElectroKinPropSolFile("electrokineticProperties" + solidsNames[i]); 

	// Set ith electrokineticProperties to dictionary of file
	SolidElectroKineticProperties.set(

		i, 
		new IOdictionary(

			IOobject(

				ElectroKinPropSolFile,
				runTime.constant(),
				solidRegions[i],
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		)
	);

// Construct dielectric values from dictionaries 
//(CHANGED TO NON-CONSTANT TO ALLOW VARYING VALUES BASED ON CONCENTRATION,ETC!!!!)
epsilonR_Solid.set(i,new dimensionedScalar(SolidElectroKineticProperties[i].lookup("epsilonR")));
epsilon_Solid.set(i,new dimensionedScalar(epsilonR_Solid[i] * epsilonZero_));

}

forAll(fluidRegions, i){

	word ElectroKinPropFluFile("electrokineticProperties" + fluidNames[i]); 

	FluidElectroKineticProperties.set(

		i,
		new IOdictionary(

			IOobject(

				ElectroKinPropFluFile,
				runTime.constant(),
				fluidRegions[i],
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		)
	);

// Construct dielectric values from dictionaries 
//(CHANGED TO NON-CONSTANT ALLOW VARYING VALUES BASED ON CONCENTRATION,ETC!!!!)
epsilonR_Fluid.set(i, new dimensionedScalar(FluidElectroKineticProperties[i].lookup("epsilonR"))); //relative dielectric
epsilon_Fluid.set(i,new dimensionedScalar(epsilonR_Fluid[i] * epsilonZero_));

}

Info<< "Creating fluid/solid fields\n" << endl;

#include "createFluidFields.H"
#include "createSolidFields.H"
