PtrList<volScalarField> fieldssolid(solidRegions.size());
PtrList<fv::options> solidFvOptions(solidRegions.size());

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh properties for region "
        << solidRegions[i].name() << nl << endl;


    Info<< "    Adding to fields solid\n" << endl;
    fieldssolid.set
    (
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

    Info<< "    Adding fvOptions\n" << endl;
    solidFvOptions.set
    (
        i,
        new fv::options(solidRegions[i])
    );

}
