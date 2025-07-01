/*---------------------------------------------------------------------------*\
Description
    Boundary condition for the mesh velocity of a free surface. Forces the mesh
    at the boundary to have an equal velocity to the fluid

\*---------------------------------------------------------------------------*/

#ifndef kinematicSurfaceVelocityPointPatchVectorField_H
#define kinematicSurfaceVelocityPointPatchVectorField_H

#include "fixedValuePointPatchField.H"
#include "surfaceVelocityTools.H"


namespace Foam
{

class kinematicSurfaceVelocityPointPatchVectorField
:
    public fixedValuePointPatchField<vector>
{
private:
    surfaceVelocityTools surfaceVelocityTools_;
    const word UName_;

public:

    //- Runtime type information
    TypeName("kinematicSurfaceVelocity");


    // Constructors

        //- Construct from patch, internal field and dictionary
        kinematicSurfaceVelocityPointPatchVectorField
        (
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  stefanMeltVelocityPointPatchVectorField
        //  onto a new patch
        kinematicSurfaceVelocityPointPatchVectorField
        (
            const kinematicSurfaceVelocityPointPatchVectorField&,
            const pointPatch&,
            const DimensionedField<vector, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        kinematicSurfaceVelocityPointPatchVectorField
        (
            const kinematicSurfaceVelocityPointPatchVectorField&,
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
                new kinematicSurfaceVelocityPointPatchVectorField
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
