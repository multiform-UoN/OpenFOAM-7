//PtrList<volScalarField> fieldssolid(solidRegions.size());

PtrList<volScalarField> VSolid(solidRegions.size());
PtrList<volScalarField> TSolid(solidRegions.size());
PtrList<surfaceScalarField> phiNPSolid(solidRegions.size());
PtrList<volVectorField> USolid(solidRegions.size());
PtrList<volScalarField> rhoSolid(solidRegions.size());

//List of electro-neutrality and diffusion dominant flags
List<bool> electroNeutralitySolid(solidRegions.size());
List<bool> diffusionDominantSolid(solidRegions.size());

//Lists of reference cells and values
List<label> VRefCellSolid(solidRegions.size());
List<scalar> VRefValueSolid(solidRegions.size());

//Pointer list towards each solid regions concentration pointer list
PtrList<PtrList<volScalarField>> cPtrLSolid(solidRegions.size());

//Pointer list of species dictionary list for solid regions
PtrList<PtrList<dictionary>> speciesSolid(solidRegions.size());


// Populate solid field pointer lists
forAll(solidRegions, i){

    //Look up species from dictionary of solid region i
    PtrList<entry> specEntries(

            SolidElectroKineticProperties[i].lookup("species")
    );

    //Set List entry to diffusion dominant flag in dictionary
    diffusionDominantSolid[i] = readBool(SolidElectroKineticProperties[i].lookup("diffusionDominant"));

    //Set List entry to electro-neutrality flag in dictionary
    electroNeutralitySolid[i] = readBool(SolidElectroKineticProperties[i].lookup("electroNeutrality"));

    //Pre-allocate PtrList of PtrLists for species dictionaries and concentrations
    cPtrLSolid.set(
        i,
        new PtrList<volScalarField>(specEntries.size())
    );

    speciesSolid.set(
        i,
        new PtrList<dictionary>(specEntries.size())
    );

    Info<< "*** Reading solid mesh properties for region "
        << solidRegions[i].name() << nl << endl;


    //VSolid := solid electric potential field
    Info<< "    Adding to VSolid\n" << endl;
    VSolid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "V",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), solidRegions[i]
        )
    );

    VRefCellSolid[i] = 0;
    VRefValueSolid[i] = 0.0;

    //Update RefCellSolid and RefValueSolid if need to
    setRefCell(
        VSolid[i],
        solidRegions[i].solutionDict().subDict("couplingControls"),
        VRefCellSolid[i],
        VRefValueSolid[i]
    );

    //TSolid := solid temperature field
    Info<<"     Adding to TSolid\n"<< endl;
    TSolid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ), solidRegions[i]
        )
    );

    //phiNPSolid := Nernst-Planck solid face flux field
    Info<<"     Adding to phiNPSolid\n"<< endl;
/*
    if (diffusionDominant)
    {
        phiNPSolid.set(
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiNP",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar("Temp",dimensionSet(0,1,0,0,0,0,0),0)
            )
        );        
    }

    else

    {*/
    phiNPSolid.set(
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phiNP",
                runTime.timeName(),
                solidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),  -fvc::flux(e_*fvc::grad(VSolid[i])/(k_*TSolid[i]))
        )
    );

    //USolid := solid velocity field
    Info<<"     Adding to USolid\n"<< endl;
    USolid.set(
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),  solidRegions[i]
        )
    );


    //rhoSolid := solid charge density field
    Info<<"     Adding to rhoSolid\n"<< endl;
    rhoSolid.set(
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                solidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i],
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

        cPtrLSolid[i].set(
            sp,
            new volScalarField
            (
                IOobject
                (
                    specEntries[sp].keyword(),
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )

        );

        speciesSolid[i].set(sp,new dictionary(specEntries[sp].dict()));

    }

}
