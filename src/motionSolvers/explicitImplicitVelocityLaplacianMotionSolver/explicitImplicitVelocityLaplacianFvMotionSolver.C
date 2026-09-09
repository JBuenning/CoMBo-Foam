#include "explicitImplicitVelocityLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMotionSolvers
{
    defineTypeNameAndDebug(explicitImplicitVelocityLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        fvMeshMover,
        explicitImplicitVelocityLaplacianFvMotionSolver,
        fvMesh
    );

    addToRunTimeSelectionTable
    (
        pointMeshMover,
        explicitImplicitVelocityLaplacianFvMotionSolver,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::explicitImplicitVelocityLaplacianFvMotionSolver::explicitImplicitVelocityLaplacianFvMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    fvMotionSolver(mesh),
    pointMeshMovers::velocity(mesh, dict, typeName),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMotionSolver::mesh(),
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    diffusivityType_(dict.lookup("diffusivity")),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMotionSolver::mesh(), diffusivityType_)
    ),
    timeIndex_(fvMotionSolver::mesh().time().timeIndex()),
    theta_(dict.lookup<scalar>("theta"))
{
    // make sure old time is saved and available in first iteration
    pointMotionU_.oldTime();
}


Foam::fvMotionSolvers::explicitImplicitVelocityLaplacianFvMotionSolver::explicitImplicitVelocityLaplacianFvMotionSolver
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    explicitImplicitVelocityLaplacianFvMotionSolver(mesh.poly(), dict)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMotionSolvers::explicitImplicitVelocityLaplacianFvMotionSolver::~explicitImplicitVelocityLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::fvMotionSolvers::explicitImplicitVelocityLaplacianFvMotionSolver::newPoints()
{

    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(mesh().points());

    diffusivityPtr_->correct();
    pointMotionU_.boundaryFieldRef().updateCoeffs();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellMotionU_,
            "laplacian(diffusivity,cellMotionU)"
        )
    );


    volPointInterpolation::New(mesh()).interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    tmp<pointField> tcurPoints;

    if (timeIndex_ != mesh().time().timeIndex())
    {
        timeIndex_ = mesh().time().timeIndex();
        tcurPoints = tmp<pointField>
        (
            mesh().points()
            + mesh().time().deltaTValue() *
            (theta_ * pointMotionU_.primitiveField() + (1. - theta_) * pointMotionU_.oldTime().primitiveField())
        );
    }
    else
    {
        tcurPoints = tmp<pointField>
        (
            mesh().oldPoints()
            + mesh().time().deltaTValue() *
            (theta_ * pointMotionU_.primitiveField() + (1. - theta_) * pointMotionU_.oldTime().primitiveField())
        );
    }

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


//void Foam::explicitImplicitVelocityLaplacianFvMotionSolver::movePoints(const pointField& p)
//{
//    // Movement of pointMesh and volPointInterpolation already
//    // done by polyMesh,fvMesh
//}


void Foam::fvMotionSolvers::explicitImplicitVelocityLaplacianFvMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    pointMeshMovers::velocity::topoChange(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityType_.rewind();
    diffusivityPtr_ = motionDiffusivity::New
    (
        mesh(),
        diffusivityType_
    );
}


void Foam::fvMotionSolvers::explicitImplicitVelocityLaplacianFvMotionSolver::mapMesh
(
    const polyMeshMap& map
)
{
    pointMeshMovers::velocity::mapMesh(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityType_.rewind();
    diffusivityPtr_ = motionDiffusivity::New
    (
        mesh(),
        diffusivityType_
    );
}
