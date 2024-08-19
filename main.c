/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   PROGRAM DESCRIPTION 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  Solving the linear advection equation in 1D:

  \(\frac{\delta u}{\delta x} = \frac{\delta u}{\delta t}, ~ 0 < x < 10\)

  With constant (Dirichlet) boundary conditions:

  \(u = 0 ~ \text{for} ~ x = 10, ~ x = 0\)

  Example of Usage:
   ./main -infile line.gmsh -outfile -output.vtk -log_view
*/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   HELP MESSAGE 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static char help[] =
    "Solving the linear advection equation in 1D\n";

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   INCLUDE STATEMENTS 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#include <petscdmplex.h>

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   USER CONTEXT 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
typedef struct {
  // Example user structure
  char     infile[PETSC_MAX_PATH_LEN];  // Input mesh filename 
  char     outfile[PETSC_MAX_PATH_LEN]; // Dump/reload mesh filename 
} AppCtx;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   FUNCTION DECLARATIONS 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options);
static PetscErrorCode CreateFiniteElements(DM dm, PetscFE fe, PetscInt dim);
static PetscErrorCode SaveDMtoVTK(AppCtx *ctx, DM dm); 
static PetscErrorCode PLACEHOLDER();

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   MAIN PROGRAM
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
int main(int argc, char **argv) {

  // Declare variables
  AppCtx      ctx;    // user program context
  DM          dm;     // mesh data managment
  PetscFE fe;         // finite element
  PetscInt dim = 1;   // dimension of the problem
  Vec u;

  // Initialize PETSc code
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  
  // Set default user-defined parameters
  PetscCall(PetscStrcpy(ctx.infile, "input.gmsh"));
  PetscCall(PetscStrcpy(ctx.outfile, "output.vtk"));

  // Process user options
  PetscCall(ProcessOptions(PETSC_COMM_WORLD, &ctx));
  
  // Read in 3D mesh from file
  PetscCall(DMPlexCreateFromFile(PETSC_COMM_WORLD, ctx.infile, NULL, PETSC_TRUE,
                                 &dm));

  // Create finite element
  CreateFiniteElements(dm, fe, dim);
  
  // Create timestepping solver context

  // Set initial conditions

  // Set boundary conditions

  // Run solver
  DMCreateGlobalVector(dm, &u);
  
  // attach DM to solver object
  //KSPSetDM(ksp, dm);

  // Save results
  SaveDMtoVTK(&ctx, dm); 

  // Free objects from memory
  
  // End main program
  PetscCall(PetscFinalize());
  return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   FUNCTIONS 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Process Options 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscBool flg;

  PetscFunctionBeginUser;
  
  PetscOptionsBegin(comm, "", "My Project's Options", "");

  // get the input file name
  PetscCall(PetscOptionsString("-infile", "The input mesh file", "",
                               options->infile, options->infile,
                               sizeof(options->infile), &flg));
  // check if it was successfully passed as an argument
  PetscCheck(flg, comm, PETSC_ERR_USER_INPUT, "-infile needs to be specified");
  // get the output file name
  PetscCall(PetscOptionsString("-outfile", "The output mesh file", "",
                               options->outfile, options->outfile,
                               sizeof(options->outfile), &flg));
  // check if it was sucessfully passed as an argument
  PetscCheck(flg, comm, PETSC_ERR_USER_INPUT, "-outfile needs to be specified");

  // End function
  PetscOptionsEnd();
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Create Finite Element 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode CreateFiniteElements(DM dm, PetscFE fe, PetscInt dim) {
  PetscBool simplex;
  
  PetscFunctionBeginUser;

  // check if DMPlex contains simplexes
  DMPlexIsSimplex(dm, &simplex);
  // Create finite element
  PetscFECreateDefault(PETSC_COMM_SELF, dim, dim, simplex, "finite_el_", -1 , &fe);
  // set one field (labled "0") with dof declared above in the DMPLEX
  DMSetField(dm, 0, NULL, (PetscObject)fe);

  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Save DM to VTK 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode SaveDMtoVTK(AppCtx *ctx, DM dm) {
  PetscViewer viewer; // object to view the mesh and solution

  PetscFunctionBeginUser;

  // give the mesh a name
  PetscCall(PetscObjectSetName((PetscObject)dm, "plexA"));
  // open the VTK file 
  PetscCall(PetscViewerVTKOpen(PETSC_COMM_WORLD, ctx->outfile, FILE_MODE_WRITE, &viewer));
  // write DM to VTK file
  PetscCall(DMView(dm, viewer));
  // destroy viewer
  PetscCall(PetscViewerDestroy(&viewer));

  PetscFunctionReturn(PETSC_SUCCESS);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   Example Function 
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static PetscErrorCode PLACEHOLDER() {
  PetscFunctionBeginUser;
  // ...
  PetscFunctionReturn(PETSC_SUCCESS);
}

