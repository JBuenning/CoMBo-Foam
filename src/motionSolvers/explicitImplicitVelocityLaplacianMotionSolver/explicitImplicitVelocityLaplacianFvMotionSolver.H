/*---------------------------------------------------------------------------*\
Description
    Mesh motion solver for an fvMesh.  Based on solving the cell-centre
    Laplacian for the motion velocity.

\*---------------------------------------------------------------------------*/

#ifndef explicitImplicitVelocityLaplacianFvMotionSolver_H
#define explicitImplicitVelocityLaplacianFvMotionSolver_H

#include "velocity_pointMeshMover.H"
#include "fvMotionSolver.H"


namespace Foam
{

// Forward class declarations
class motionDiffusivity;

namespace fvMotionSolvers
{

class explicitImplicitVelocityLaplacianFvMotionSolver
:
    public fvMotionSolver,
    public pointMeshMovers::velocity
{
    // Private Data

        //- Cell-centre motion field
        mutable volVectorField cellMotionU_;

        //- Diffusivity type and data stream
        ITstream diffusivityType_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        mutable label timeIndex_;
        scalar theta_;


public:

    //- Runtime type information
    TypeName("explicitImplicitVelocityLaplacian");


    // Constructors


        //- Construct from fvMesh and dictionary
        explicitImplicitVelocityLaplacianFvMotionSolver
        (
            fvMesh&,
            const dictionary&
        );

        //- Construct from polyMesh and dictionary
        explicitImplicitVelocityLaplacianFvMotionSolver
        (
            const polyMesh&,
            const dictionary&
        );

        //- Disallow default bitwise copy construction
        explicitImplicitVelocityLaplacianFvMotionSolver
        (
            const explicitImplicitVelocityLaplacianFvMotionSolver&
        ) = delete;


    //- Destructor
    ~explicitImplicitVelocityLaplacianFvMotionSolver();


    // Member Functions

        //- Return reference to the cell motion velocity field
        volVectorField& cellMotionU()
        {
            return cellMotionU_;
        }

        //- Return const reference to the cell motion velocity field
        const volVectorField& cellMotionU() const
        {
            return cellMotionU_;
        }

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> newPoints();

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const explicitImplicitVelocityLaplacianFvMotionSolver&) = delete;
};

} // End namespace fvMotionSolvers
} // End namespace Foam

#endif
