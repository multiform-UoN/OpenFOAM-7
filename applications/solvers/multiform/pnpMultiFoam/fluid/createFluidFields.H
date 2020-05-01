//PtrList<volScalarField> fieldsfluid(fluidRegions.size());

PtrList<volScalarField> VFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiNPFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());

//List of electro-neutrality and diffusion dominant flags
List<bool> electroNeutralityFluid(fluidRegions.size());
List<bool> diffusionDominantFluid(fluidRegions.size());

//Lists of reference cells and values
List<label> VRefCellFluid(fluidRegions.size());
List<scalar> VRefValueFluid(fluidRegions.size());

//Pointer list towards each fluid regions concentration pointer list
PtrList<PtrList<volScalarField>> cPtrLFluid(fluidRegions.size());

//Pointer list of species dictionary list for fluid regions
PtrList<PtrList<dictionary>> speciesFluid(fluidRegions.size());


// Populate fluid field pointer lists
forAll(fluidRegions, i){

    //Look up species from dictionary of fluid region i
    PtrList<entry> specEntries(

            FluidElectroKineticProperties[i].lookup("species")
        );

    //Set List entry to diffusion dominant flag in dictionary
    diffusionDominantFluid[i] = readBool(FluidElectroKineticProperties[i].lookup("diffusionDominant"));

    //Set List entry to electro-neutrality flag in dictionary
    electroNeutralityFluid[i] = readBool(FluidElectroKineticProperties[i].lookup("electroNeutrality"));

    //Pre-allocate PtrList of PtrLists for species dictionaries and concentrations
    cPtrLFluid.set(
        i,
        new PtrList<volScalarField>(specEntries.size())
    );

    speciesFluid.set(
        i,
        new PtrList<dictionary>(specEntries.size())
    );

    Info<< "*** Reading fluid mesh properties for region "
        << fluidRegions[i].name() << nl << endl;


    //VFluid := fluid electric potential field
    Info<< "    Adding to VFluid\n" << endl;
    VFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "V",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), fluidRegions[i]
        )
    );

    VRefCellFluid[i] = 0;
    VRefValueFluid[i] = 0.0;

    //Update RefCellFluid and RefValueFluid if need to
    setRefCell(
        VFluid[i],
        fluidRegions[i].solutionDict().subDict("couplingControls"),
        VRefCellFluid[i],
        VRefValueFluid[i]
    );

    //TFluid := fluid temperature field
    Info<<"     Adding to TFluid\n"<< endl;
    TFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), fluidRegions[i]
        )
    );

    //phiNPFluid := Nernst-Planck fluid face flux field
    Info<<"     Adding to phiNPFluid\n"<< endl;
/*
    //Set phiNPFluid (i.e the electrostatic flux term) to zero if diffusion dominant
    if (diffusionDominant)
    {
        phiNPFluid.set(
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiNP",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i],
                dimensionedScalar("Temp",dimensionSet(0,1,0,0,0,0,0),0)
            )
        );   
    }

    else*/
    phiNPFluid.set(
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phiNP",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),  -fvc::flux(e_*fvc::grad(VFluid[i])/(k_*TFluid[i]))
        )
    );


    //UFluid := fluid velocity field
    Info<<"     Adding to UFluid\n"<< endl;
    UFluid.set(
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),  fluidRegions[i]
        )
    );

    //rhoFluid := fluid charge density field
    Info<<"     Adding to rhoFluid\n"<< endl;
    rhoFluid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(
                "rho",
                dimensionSet(0,-3,1,0,0,1,0),
                scalar(0)
            )
        )
    );


    //- Set all species objects in list
    for(int sp=0; sp<specEntries.size(); sp++)
    {

        cPtrLFluid[i].set(
            sp,
            new volScalarField
            (
                IOobject
                (
                    specEntries[sp].keyword(),
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )

        );

        speciesFluid[i].set(sp,new dictionary(specEntries[sp].dict()));

    }
}
