/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  7
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      C.Cl;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [ 0 -3 0 0 1 0 0 ];

c_bulk 25.6e-6; //I think this is 25.6e-3 (from Wang Yixian email chain)

internalField   uniform $c_bulk;

boundaryField
{
    Inner_Pip_Wall
    {
        type            noFlux;//zeroGradient;
        value           uniform 0;
    }

    Outer_Pip_Wall
    {
        type            noFlux;//zeroGradient;
        value           uniform 0;
    }

    Inner_Electrode
    {
        type            noFlux;//fixedValue;
        value           uniform $c_bulk;
    }

    Planar_Electrode
    {
        type            noFlux;//fixedValue;
        value           $internalField;
    }

    Outer_Walls
    {
        type            fixedValue;
        value           uniform $c_bulk;//internalField;
    }

#includeEtc "caseDicts/setConstraintTypes"

}

// ************************************************************************* //
