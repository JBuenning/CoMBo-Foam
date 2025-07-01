/*---------------------------------------------------------------------------*\
Description
    Tools for the programming of mesh velocity boundary conditions that
    calculate mesh point velocities from face velocities.
    Handles face-to-point interpolation and laplace smoothing
Usage
    Adds the following properties to the boudnary condition
    \table
        Property                    | Description           | Required | Default value
        URef                        | velocity the mesh velocity should be \\
          relative to (useful for moving objects)           | no       | (0 0 0)
        linearUpwindBlendingFactor  | blending between linear and upwind for \\
          the face-to-point interpolation (1 corresponds to full linear and \\
          0 corresponds to full upwind interpolation)       | yes      |
        restrictOnPatches           | don't allow point movement in the normal \\
          direction of the specified patches                | no       | ()
        laplaceSmoothing            | apply lapalace smoothing? | no   | no
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef surfaceVelocityTools_H
#define surfaceVelocityTools_H

#include "autoPtr.H"
#include "fieldMapper.H"
#include "pointFields.H"
#include "pointPatch.H"
#include "polyPatch.H"
#include "scalarField.H"
#include "vectorField.H"
#include "surfaceFields.H"
#include "surfaceVelocityInterpolation.H"


namespace Foam
{

class surfaceVelocityTools
{
private:
    const scalar linearUpwindBlendingFactor_;
    const bool laplaceSmoothing_;
    label timeIndex_;
    vectorField laplacePointCorrection_;
    autoPtr<UpwindSurfaceVelocityInterpolation> upwindInterp_;

    void updateLaplaceCorrection(const polyPatch& pp);
    vectorField boundaryConformingPointNormals(const polyPatch& pp) const;

public:

    // Constructors

        surfaceVelocityTools
        (
            const pointPatch&,
            const dictionary&
        );

        surfaceVelocityTools
        (
            const surfaceVelocityTools&,
            const pointPatch&
        );

        surfaceVelocityTools
        (
            const surfaceVelocityTools&
        );


    // Member Functions

        void map
        (
            const surfaceVelocityTools&,
            const fieldMapper&
        );

        void reset(const surfaceVelocityTools&);


        vectorField pointVelocityFromFaceVelocity
        (
            const scalarField& faceNormalVelocity,
            const polyPatch& pp
        );

        void write(Ostream&) const;
};


} // End namespace Foam

#endif
