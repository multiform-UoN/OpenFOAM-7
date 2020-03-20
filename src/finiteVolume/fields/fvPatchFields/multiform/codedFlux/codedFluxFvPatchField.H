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
    Foam::codedFluxFvPatchField

Description
    Constructs on-the-fly a new boundary condition (derived from
    fluxFvPatchField) which is then used to evaluate.

Usage
    Example:
    \verbatim
    <patchName>
    {
        type            codedFlux;

        value           uniform (0 0 0);
        fluxICoeffs     uniform (0 0 0);
        fluxBCoeffs     uniform 1;

        name    rampedFlux;   // name of generated BC

        code
        #{
            this->fluxBoundaryCoeffs() =
                vector(1, 0, 0)
               *min(10, 0.1*this->db().time().value());
            this->fluxInternalCoeffs() = Zero;
            this->snGradRef() = Zero;
            operator==(1.0);
        #};

        // codeInclude
        //#{
        //    #include "fvCFD.H"
        //#};

        // codeOptions
        //#{
        //    -I$(LIB_SRC)/finiteVolume/lnInclude
        //#};
    }
    \endverbatim

    A special form is if the 'code' section is not supplied. In this case
    the code gets read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry

    \verbatim
    <patchName>
    {
        code
        #{
            this->fluxBoundaryCoeffs() = min(10, 0.1*this->db().time().value());
            this->fluxInternalCoeffs() = Zero;
            this->snGradRef() = Zero;
            operator==(1.0);
        #};
    }
    \endverbatim

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    codedFluxFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef codedFluxFvPatchField_H
#define codedFluxFvPatchField_H

#include "fluxFvPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dynamicCode;
class dynamicCodeContext;
class IOdictionary;

/*---------------------------------------------------------------------------*\
                   Class codedFluxFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedFluxFvPatchField
:
    public fluxFvPatchField<Type>,
    public codedBase
{
    // Private static data

        //- The keywords associated with source code
        static const wordList codeKeys_;


    // Private Data

        //- Dictionary contents for the boundary condition
        mutable dictionary dict_;

        const word name_;

        mutable autoPtr<fluxFvPatchField<Type>> redirectPatchFieldPtr_;


    // Private Member Functions

        const IOdictionary& dict() const;

        //- Set the rewrite vars controlling the Type
        static void setFieldTemplates(dynamicCode& dynCode);

        //- Get the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;

        //- Return a description (type + name) for the output
        virtual string description() const;

        //- Clear the ptr to the redirected object
        virtual void clearRedirect() const;

        //- Get the dictionary to initialize the codeContext
        virtual const dictionary& codeDict() const;

        //- Get the keywords associated with source code
        virtual const wordList& codeKeys() const;


public:

    // Static Data Members

        //- Name of the C code template to be used
        static const word codeTemplateC;

        //- Name of the H code template to be used
        static const word codeTemplateH;


    //- Runtime type information
    TypeName("codedFlux");


    // Constructors

        //- Construct from patch and internal field
        codedFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedFluxFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedFluxFvPatchField
        //  onto a new patch
        codedFluxFvPatchField
        (
            const codedFluxFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        codedFluxFvPatchField
        (
            const codedFluxFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new codedFluxFvPatchField<Type>(*this)
            );
        }

        //- Copy constructor setting internal field reference
        codedFluxFvPatchField
        (
            const codedFluxFvPatchField<Type>&,
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
                new codedFluxFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Get reference to the underlying patchField
        const fluxFvPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        //  This is only needed to set the updated() flag of the name
        //  to false.
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "codedFluxFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
