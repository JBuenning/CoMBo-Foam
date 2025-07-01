/*---------------------------------------------------------------------------*\
Description
    Foam::onePhaseStefanMeltVelocityPointPatchVectorField

\*---------------------------------------------------------------------------*/

#ifndef onePhaseStefanMeltVelocityPointPatchVectorField_H
#define onePhaseStefanMeltVelocityPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "surfaceVelocityTools.H"


namespace Foam
{

class onePhaseStefanMeltVelocityPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
private:
    surfaceVelocityTools surfaceVelocityTools_;
    const word TName_;
    const scalar kappaOverRhoH_;
    const vector URef_;

public:

    //- Runtime type information
    TypeName("onePhaseStefanMeltVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        onePhaseStefanMeltVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  stefanMeltVelocityPointPatchVectorField
        //  onto a new patch
        onePhaseStefanMeltVelocityPointPatchVectorField
        (
            const onePhaseStefanMeltVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        onePhaseStefanMeltVelocityPointPatchVectorField
        (
            const onePhaseStefanMeltVelocityPointPatchVectorField&,
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
                new onePhaseStefanMeltVelocityPointPatchVectorField
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
                const pointPatchFieldMapper&
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
