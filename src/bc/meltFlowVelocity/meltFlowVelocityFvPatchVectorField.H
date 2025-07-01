/*---------------------------------------------------------------------------*\
Description
    This boundary condition provides an inflow condition for boundaries
    between two phases of a phase-change material. The melt velocity is
    assumed to be given at the boundary of the field pointMotionU
\*---------------------------------------------------------------------------*/

#ifndef meltFlowVelocityFvPatchVectorField_H
#define meltFlowVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

namespace Foam
{

class meltFlowVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
private:
    const scalar rhoRatio_;
    const vector URef_;

public:

    //- Runtime type information
    TypeName("meltFlowVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        meltFlowVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        meltFlowVelocityFvPatchVectorField
        (
            const meltFlowVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fieldMapper&
        );

        //- Disallow copy without setting internal field reference
        meltFlowVelocityFvPatchVectorField
        (
            const meltFlowVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        meltFlowVelocityFvPatchVectorField
        (
            const meltFlowVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new meltFlowVelocityFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given fvPatchField onto this fvPatchField
            virtual void map(const fvPatchVectorField&, const fieldMapper&);

            //- Reset the fvPatchField to the given fvPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const fvPatchVectorField&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif
