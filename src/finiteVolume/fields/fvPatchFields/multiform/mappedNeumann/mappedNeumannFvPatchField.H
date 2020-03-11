/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

Class
    Foam::mappedNeumannFvPatchField

Description
    This boundary condition sets the gradient according to the the gradient in the mapped patch

    The sample mode is set by the underlying mapping engine, provided by the
    mappedPatchBase class.

Usage
    \table
        Property     | Description             | Required    | Default value
        fieldName    | name of field to be mapped | no       | this field name
        setAverage   | flag to activate setting of average value | yes |
        average      | average value to apply if \c setAverage = yes | yes |
        interpolationScheme | type of interpolation scheme | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            mapped;
        fieldName       T;
        setAverage      no;
        average         0;
        interpolationScheme cell;
        value           uniform 0;
    }
    \endverbatim

    When employing the \c nearestCell sample mode, the user must also specify
    the interpolation scheme using the \c interpolationScheme entry.

    In case of interpolation (where scheme != cell) the limitation is that
    there is only one value per cell.  For example, if you have a cell with two
    boundary faces and both faces sample into the cell, both faces will get the
    same value.

Note
    It is not possible to sample internal faces since volume fields are not
    defined on faces.

See also
    Foam::mappedPatchBase
    Foam::interpolation
    Foam::fixedValueFvPatchField

SourceFiles
    mappedNeumannFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef mappedNeumannFvPatchField_H
#define mappedNeumannFvPatchField_H

#include "fixedGradientFvPatchFields.H"
#include "newMappedPatchFieldBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class mappedNeumannFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class mappedNeumannFvPatchField
:
    public fixedGradientFvPatchField<Type>,
    public newMappedPatchFieldBase<Type>
{

protected:

    // Protected Member Functions

        const mappedPatchBase& mapper
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>& iF
        );


public:

    //- Runtime type information
    TypeName("mappedNeumann");


    // Constructors

        //- Construct from patch and internal field
        mappedNeumannFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        mappedNeumannFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mappedNeumannFvPatchField
        //  onto a new patch
        mappedNeumannFvPatchField
        (
            const mappedNeumannFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        mappedNeumannFvPatchField
        (
            const mappedNeumannFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedNeumannFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        mappedNeumannFvPatchField
        (
            const mappedNeumannFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new mappedNeumannFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "mappedNeumannFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
