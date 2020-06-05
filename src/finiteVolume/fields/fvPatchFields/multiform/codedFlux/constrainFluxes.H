#include "fluxFvPatchFields.H"
namespace Foam{

template<class Type>
void constrainFluxes(fvMatrix<Type>& eqn)
{
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        psiBf(eqn.psi().boundaryField());

    forAll(psiBf, patchi)
    {
        if (isA<fluxFvPatchField<Type>>(psiBf[patchi]))
        {
            const fluxFvPatchField<Type>& fluxPatch
            (
                refCast<const fluxFvPatchField<Type>>
                (
                    psiBf[patchi]
                )
            );

            eqn.internalCoeffs()[patchi] =
                fluxPatch.patch().magSf() * fluxPatch.fluxInternalCoeffs();

            eqn.boundaryCoeffs()[patchi] =
                fluxPatch.patch().magSf() * fluxPatch.fluxBoundaryCoeffs();

            fluxPatch.fluxCorrected();
        }
    }
}
}
