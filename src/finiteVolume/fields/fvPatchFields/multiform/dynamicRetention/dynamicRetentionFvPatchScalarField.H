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
    Foam::dynamicRetentionFvPatchScalarField

Description
    This boundary condition provide a reversible reactive flux for a scalar
    concentration field.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    \endverbatim


SourceFiles
    dynamicRetentionFvPatchScalarField.C

Authors:
    Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/

#ifndef dynamicRetentionFvPatchScalarField_H
#define dynamicRetentionFvPatchScalarField_H

#include "RobinFvPatchFields.H"
#include "retentionRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class dynamicRetentionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class dynamicRetentionFvPatchScalarField
:
    public RobinFvPatchScalarField
{
private:

    //- Boundary concentration
    scalarField S_;

    //- Boundary concentration
    //  at previous time step
    scalarField S0_;

    //- Effective RobinK_
    scalarField RobinKeff_;

    //- Reversible reaction rate Kd
    scalarField Kd_;

    //- Switch to activate shear induced detachment
    Switch  shearDetachment_;

    //- Effective RobinF_
    scalarField RobinFeff_;

    //- Current time index.
    //  Used to trigger the storing of the old-time value
    mutable label timeIndex_;

    //- Pointer to retention rate
    autoPtr<retentionRates::retentionRate> rR_;

    //- Pointer to detachment rate
    autoPtr<retentionRates::retentionRate> dR_;

	//- Pointer to shear induced detachment rate
    autoPtr<retentionRates::retentionRate> sdR_;

    scalar SSp_;
    scalar SSu_;

public:

    //- Runtime type information
    TypeName("dynamicRetention");


    // Constructors

        //- Construct from patch and internal field
        dynamicRetentionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dynamicRetentionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given dynamicRetentionFvPatchScalarField
        //  onto a new patch
        dynamicRetentionFvPatchScalarField
        (
            const dynamicRetentionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        dynamicRetentionFvPatchScalarField
        (
            const dynamicRetentionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicRetentionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dynamicRetentionFvPatchScalarField
        (
            const dynamicRetentionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new dynamicRetentionFvPatchScalarField(*this, iF)
            );
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
                 const fvPatchField<scalar>&,
                 const labelList&
             );


        // Member functions

            void SSp(dimensionedScalar source)
            {
                SSp_ = source.value();
            }

            void SSu(dimensionedScalar source)
            {
                SSu_ = source.value();
            }

            scalar Savg() const
            {
                return gSum(this->patch().magSf()*S_)/gSum(this->patch().magSf());
            }

            void rescaleS(scalar avg)
            {
                S_ = S_*avg/Savg();
            }

            virtual scalarField& RobinF()
            {
                return RobinFeff_;
            }

            virtual const scalarField& RobinF() const
            {
                return RobinFeff_;
            }

            virtual scalarField& RobinK()
            {
                return RobinKeff_;
            }

            virtual const scalarField& RobinK() const
            {
                return RobinKeff_;
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            // //- Evaluate the patch field
            // virtual void updateCoeffs();

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
