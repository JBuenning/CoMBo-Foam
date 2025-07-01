/*---------------------------------------------------------------------------*\
Description
    Foam::twoPhaseStefanMeltVelocityPointPatchVectorField

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseStefanMeltVelocityPointPatchVectorField_H
#define twoPhaseStefanMeltVelocityPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "surfaceVelocityTools.H"


namespace Foam
{

class twoPhaseStefanMeltVelocityPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
private:
    surfaceVelocityTools surfaceVelocityTools_;
    const word TName_;
    const scalar kappa_;
    const scalar rhoTimesH_;
    const vector URef_;
    bool externally_updated_;

public:

    //- Runtime type information
    TypeName("twoPhaseStefanMeltVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        twoPhaseStefanMeltVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  twoPhaseStefanMeltVelocityPointPatchVectorField
        //  onto a new patch
        twoPhaseStefanMeltVelocityPointPatchVectorField
        (
            const twoPhaseStefanMeltVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        twoPhaseStefanMeltVelocityPointPatchVectorField
        (
            const twoPhaseStefanMeltVelocityPointPatchVectorField&,
            const DimensionedField<vector, pointMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<vector>> clone
        (
            const DimensionedField<vector, pointMesh>& iF
        ) const
        {
            return autoPtr<pointPatchField<vector>>
            (
                new twoPhaseStefanMeltVelocityPointPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map the given pointPatchField onto this pointPatchField
            virtual void map
            (
                const pointPatchField<vector>&,
                const fieldMapper&
            );

            //- Reset the pointPatchField to the given pointPatchField
            //  Used for mesh to mesh mapping
            virtual void reset(const pointPatchField<vector>&);


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif
