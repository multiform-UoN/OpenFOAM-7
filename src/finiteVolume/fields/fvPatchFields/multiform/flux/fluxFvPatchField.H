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



Class
    Foam::fluxFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a base class for 'flux' type boundary
    conditions. This simply allows to specify the boundary flux for conservation
    laws:
        \f[
            \frac{d u}{ d t} = \nabla \cdot  J
        \f]
    Allowing to specify the patch normal flux \f[$ J_p \f]$ as a linear function
    of the internal cell values:

        \f[
            J_p = a u_c + b
        \f]

    Notice that suitable values of \[$ u_p \]$ surface normal gradients have to
    be specified and are not calculated within this class.

    \vartable
        u_p   | patch values
        u_c   | patch internal cell values
        J_p   | patch flux
        a     | flux internal coefficient
        b     | flux boundary coefficient
    \endvartable


Usage
    Basically, this BC requires the user to specify everything. Use the coded
    version to implement specific expressions for each coefficient.
    \table
        Property     | Description               | Required    | Default value
        fluxICoeffs  | flux internal coefficient | yes         |
        fluxBCoeffs  | flux boundary coefficient | yes         |
        value        | patch value               | yes         |
        snGrad       | surface normal gradient   | yes         |
    \endtable

SourceFiles
    fluxFvPatchField.C

Author
    Fedrico Municchi, Nottingham (2019)

\*---------------------------------------------------------------------------*/

#ifndef fluxFvPatchField_H
#define fluxFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class fluxFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fluxFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        //- Internal coefficient for implicit flux
        scalarField fluxICoeffs_;

        //- Boundary coefficient for explicit flux
        Field<Type> fluxBCoeffs_;

        //- Surface normal gradient (optional)
        Field<Type> snGrad_;

        //- Switch to ensure correct usage of this boundary condition
        mutable Switch fluxCorrected_;

public:

    //- Runtime type information
    TypeName("flux");


    // Constructors

        //- Construct from patch and internal field
        fluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given fluxFvPatchField onto a new patch
        fluxFvPatchField
        (
            const fluxFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy
        fluxFvPatchField
        (
            const fluxFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fluxFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fluxFvPatchField
        (
            const fluxFvPatchField<Type>&,
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
                new fluxFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return false: this patch field is not altered by assignment
            virtual bool assignable() const
            {
                return true;
            }


        // Return defining fields

            virtual scalarField& fluxInternalCoeffs()
            {
                return fluxICoeffs_;
            }

            virtual const scalarField& fluxInternalCoeffs() const
            {
                return fluxICoeffs_;
            }

            virtual Field<Type>& fluxBoundaryCoeffs()
            {
                return fluxBCoeffs_;
            }

            virtual const Field<Type>& fluxBoundaryCoeffs() const
            {
                return fluxBCoeffs_;
            }

            //- Returns surface normal gradient as reference.
            //  This is necessary for the coded version to allow modification
            //  of snGrad_.
            virtual Field<Type>& snGradRef()
            {
                return snGrad_;
            }

            virtual const Field<Type>& snGradRef() const
            {
                return snGrad_;
            }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<Type>> snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            virtual void updateCoeffs();

            virtual void fluxCorrected() const
            {
                fluxCorrected_ = true;
            }

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

        virtual void operator=(const UList<Type>&) {}

        virtual void operator=(const fvPatchField<Type>&) {}
        virtual void operator+=(const fvPatchField<Type>&) {}
        virtual void operator-=(const fvPatchField<Type>&) {}
        virtual void operator*=(const fvPatchField<scalar>&) {}
        virtual void operator/=(const fvPatchField<scalar>&) {}

        virtual void operator+=(const Field<Type>&) {}
        virtual void operator-=(const Field<Type>&) {}

        virtual void operator*=(const Field<scalar>&) {}
        virtual void operator/=(const Field<scalar>&) {}

        virtual void operator=(const Type&) {}
        virtual void operator+=(const Type&) {}
        virtual void operator-=(const Type&) {}
        virtual void operator*=(const scalar) {}
        virtual void operator/=(const scalar) {}

        //- All these functions return zero because standard boundary
        //  discretisation of differential operators is disabled.
        //  This is somehow equivalent to skipping this patch when
        //  calculating the boundary terms for differential operators.
        virtual tmp<Field<Type>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const
        {
            return tmp<Field<Type>>
            (
                new Field<Type>(this->size(), Zero)
            );
        }

        virtual tmp<Field<Type>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const
        {
            return tmp<Field<Type>>
            (
                new Field<Type>(this->size(), Zero)
            );
        }

        virtual tmp<Field<Type>> gradientInternalCoeffs() const
        {
            return tmp<Field<Type>>
            (
                new Field<Type>(this->size(), Zero)
            );
        }

        virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
        {
            return tmp<Field<Type>>
            (
                new Field<Type>(this->size(), Zero)
            );
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fluxFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
