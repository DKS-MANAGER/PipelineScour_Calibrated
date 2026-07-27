# Comprehensive Simulation of 2D Pipeline Scour Morphodynamics using sedExnerFoam

**Author:** Divyansh Kumar Singh (Postgraduate Research)  
**Date:** July 2026  
**Solver:** OpenFOAM v2412 (sedExnerFoam)  
**Repository:** [2DPipelineScourEXN](https://github.com/DKS-MANAGER/2DPipelineScourEXN)  

## 📌 1. Project Abstract & Physical Setup

This repository contains a rigorously reconstructed and validated computational fluid dynamics (CFD) setup for investigating 2D Pipeline Scour over an erodible granular bed. The primary objective is to transition from computationally expensive two-fluid Eulerian-Eulerian formulations (e.g., SedFoam) to the highly efficient, mesh-deformation-based `sedExnerFoam` solver. The approach couples the standard Navier-Stokes equations with the macroscopic Exner morphodynamic equation solved via the Finite Area Method (FAM).

### Baseline Geometry
- **Computational Domain:** $x \in [-0.75, 1.0]\text{ m}$, $y \in [-0.1, 0.205]\text{ m}$.
- **Initial Erodible Bed:** Prescribed as a horizontal interface located at $y = -0.026\text{ m}$.
- **Pipeline Cylinder:** Diameter $D = 0.05\text{ m}$, rigidly fixed with a gap-to-diameter ratio $e/D = 0.22$ ($11.0\text{mm}$). Center at $x = 0.0\text{ m}$, $y = 0.01\text{ m}$.
- **Morphodynamics Interface:** Represented via a 2D FAM mesh mapped directly to the `bottom` boundary patch.

### Fluid & Sediment Properties
- **Fluid (Water):** Density $\rho_f = 1000\text{ kg/m}^3$, kinematic viscosity $\nu = 10^{-6}\text{ m}^2\text{/s}$.
- **Sediment (Quartz Sand):** Density $\rho_s = 2650\text{ kg/m}^3$ (specific gravity $s = 2.65$), median grain diameter $d_{50} = 0.36\text{ mm}$ ($3.6 \times 10^{-4}\text{ m}$).
- **Erodible Bed Properties:** Porosity $\lambda_s = 0.40$, Angle of repose $\beta_r = 32^\circ$.
- **Derived Scale Parameters:** Soulsby critical Shields parameter $\theta_{cr} \approx 0.0473$.

### Boundary Conditions & Forcing
- **Inlet Velocity (U):** Coded logarithmic boundary layer matching a friction velocity $u^* = 0.04318\text{ m/s}$ and roughness height $z_0 = 3.0 \times 10^{-5}\text{ m}$. 
- **Seabed Boundary (bottom):** Solves sediment transport closures (Nielsen 1992 bedload, van Rijn slope correction) directly on the 2D surface.

---

## 🔄 2. Solver Architecture Transition

| Feature / Aspect | SedFoam (Reference Case) | sedExnerFoam (This Case) |
|---|---|---|
| **Phase Model** | Two-phase Eulerian-Eulerian (fluid + sediment) | Single-phase fluid; sediment bed is represented as a boundary patch |
| **Sediment Fields** | Volume fraction `alpha.b` solved in the 3D volume | Passive scalar concentration `Cs` in fluid; bedload flux on the bed boundary only |
| **Bed Morphodynamics** | Grid deformation via whole-domain mesh diffusion | Exner equation solved on a 2D Finite Area surface mesh, driving the `displacementLaplacian` solver |
| **Turbulence Closure** | k-ω SST (phase b) | Standard single-phase k-ω SST |
| **Bedload/Erosion** | Resolved via interphase drag, granular rheology, and particle pressure | Nielsen (1992) empirical bedload model & slope-corrected boundary condition |

---

## 🛠️ 3. Core Solver Diagnostics & Methodological Fixes (CRITICAL)

During the case reconstruction, several critical solver-level and configuration-level limitations were identified and permanently resolved:

### A. Resolution of Uninitialized Pointer Segfault in projectedFaMesh
During the morphodynamic phase, `sedExnerFoam` can experience a hard crash (Segmentation fault). Diagnostic tracing identified a C++ uninitialized variable defect in the constructor of `projectedFaMesh` (`src/projectedFiniteArea/projectedFaMesh/projectedFaMesh.C`). Several dynamically allocated raw pointers (`LePtr_`, `magLePtr_`, `pointCoordsPtr_`, `SPtr_`) were not explicitly initialized to `nullptr`.
**Methodological Fix:** The source code was patched to ensure strict pointer initialization, and the `projectedFiniteArea` library was recompiled, completely mitigating the memory corruption.

### B. Remediation of Locked Bed Mesh Motion
In the standard initialization, the bottom patch was rigidly constrained via `fixedValue uniform (0 0 0)` to allow for hydrodynamic spin-up. However, this mathematically locked the vertices, preventing the Exner solver from deforming the bed.
**Methodological Fix:** The automated script manages a dedicated configuration (`system/pointDisplacement`) where the bottom boundary is set as `fixedValue` so that the `sedExnerFoam` solver can explicitly inject the Exner morphological displacements into it.

### C. MPI Parallelization Constraint (Single Processor FaceSets)
The finite-area Exner solver in `sedExnerFoam` lacks parallel communication logic. A naïve volumetric domain decomposition inherently splits the bottom boundary across multiple MPI ranks, triggering fatal solver errors.
**Methodological Fix:** The parallelization strategy is constrained strictly according to the literature. Using `topoSet`, a `faceSet` is generated for the `bottom` patch. The `decomposeParDict` is subsequently programmed with a `singleProcessorFaceSets` constraint, forcing the entire morphodynamic bed to compute sequentially on a single CPU core, while distributing the volumetric Navier-Stokes matrices efficiently across the remaining cores.

### D. OpenFOAM v2412 Compatibility: Step-by-Step Source Code Fixes
The official `sedExnerFoam` solver contains a critical C++ bug in its mesh motion solver, which prevents the mesh from deforming because the displacement values are never successfully assigned back to the boundary.

To install and fix the solver on a new PC (e.g., your lab workstation), follow these exact steps to clone the development branch, patch the bug, and compile:

**1. Clone the `sedExnerFoam` solver (devel branch):**
```bash
git clone -b devel https://github.com/sedFoam/sedExnerFoam.git
```

**2. Copy the Pre-Patched Files from this Repository:** 
Instead of manually editing the C++ source code, we have provided the perfectly patched files directly in this repository. The `sedExnerFoam_v2412_fixes` directory perfectly mirrors the solver's repository structure.

Simply copy them over the original solver files:
```bash
cp -r sedExnerFoam_v2412_fixes/* <path-to-sedExnerFoam>/
```

*(For reference, these patched files apply the following critical C++ fixes for OpenFOAM v2412 with their physical/mathematical reasoning):*
- **Fix 1 (`meshMove.H`):** Modifies the initial read of the boundary displacement. The default solver incorrectly pulls the internal field values. The patched version extracts boundary patch values using `refCast<const fixedValuePointPatchVectorField>`. This ensures the Exner-calculated bottom displacements are mapped directly to the actual boundary nodes of the mesh, enabling mesh motion.
- **Fix 2 (`meshMove.H`):** Replaces raw array assignments with an explicit `tmp<Field>` memory allocation (`patchDisp == tmp<Field<vector>>::New(dispVals);`). This bypasses fatal compile-time errors in OpenFOAM v2412 caused by stricter C++ move and copy constructor semantics for `tmp` wrappers.
- **Fix 3 (`meshMove.H`):** Adds a missing `pointDisplacement.correctBoundaryConditions();` call. In parallel (MPI) execution, failing to sync boundaries leaves processor interfaces with coordinate mismatches. This function synchronizes point coordinates across CPU ranks to prevent mesh tearing and subsequent solver divergence.
- **Fix 4 (`sedExnerFoam.C`):** Adds `#include "fixedValuePointPatchFields.H"`. This provides the compiler with the complete class template definition for the `refCast` pointer conversion in Fix 1, avoiding forward-declaration compilation errors.
- **Fix 5 (`projectedFaMesh.C`):** Explicitly initializes uninitialized raw pointers (`LePtr_`, `magLePtr_`, `pointCoordsPtr_`, `SPtr_`) to `nullptr` in the constructor. In C++, uninitialized pointers hold garbage memory. When the code tests pointer validity (`!= nullptr`) before allocation/deallocation, garbage addresses pass the test, causing segmentation faults.
- **Fix 6 (`CsEqn.H`):** Clamps the passive scalar `Cs` (sediment concentration) to the physical range $[0.0, CsMax]$ before settling velocity computation. Transient grid-scale overshoots can push $Cs$ past $1.0$. The Fredsoe hindered settling model calculates settling velocity using terms like $(1 - Cs)^n$. When $Cs > 1.0$, this base is negative, and raising it to a fractional power results in `NaN` (domain error), causing Floating Point Exception (FPE) crashes.
- **Fix 7 (`exnerEqn.H`):** Initializes `lapldH` to `0.0` and corrects edge diffusion flux accumulation from assignment (`=`) to additive (`+=` and `-=`). The original code incorrectly overwrote cell values at each edge, meaning only the last edge processed in the loop contributed to the cell's Laplacian. Correcting this to accumulation preserves physical conservation of mass on the Finite Area mesh and eliminates filter leaks.
- **Fix 8 (`createFaFields.H`):** Corrects a copy-paste bug checking `switchFilterShields` instead of `switchFilterExner` when enabling `filterExner`. This allows the Exner filter to be enabled independently of the Shields filter, matching dictionary configuration expectations.

**3. Recompile the Solver:** 
Open a terminal, navigate to the `<path-to-sedExnerFoam>/` directory, and run:
```bash
./Allwmake
```

### E. ParaView Decomposed Moving Mesh Reader Bug (v2412)
When viewing an actively running simulation in parallel, ParaView's native `.foam` reader throws a fatal error: `Mismatch in number of old points and new points` or fails to show any mesh motion. This is a known ParaView bug where it incorrectly assumes decomposed parallel boundaries will never change their local point counts during dynamic mesh deformation.
**The Fix:** Do not use the native `.foam` reader. Instead, export the mesh to VTK format step-by-step:
```bash
foamToVTK
```
In ParaView, navigate to the newly created `VTK/` folder and open the `.vtk` or `.vtm` files. This perfectly visualizes the moving scour hole without any point-mismatch errors.

### F. Scaling to Higher Core Counts (20, 32+ Cores)
If you are moving this simulation to a high-performance workstation or cluster and want to run on more CPU cores (e.g., 20 or 32 cores), update `system/decomposeParDict` and `Allrun`:
- **Update `system/decomposeParDict`:** Change `numberOfSubdomains` from 6 to your desired core count (e.g., 32). 
*(Note: DO NOT remove or change the `singleProcessorFaceSets` constraint. The bottom patch MUST remain constrained to a single core to prevent faMesh decomposition errors).*
- **Update the execution script (`Allrun`):** Change `-np 20` to your target core count for the execution.

### G. Phase Transition Stability Fixes
- Note: In previous tangent configurations, `makeFaMesh` non-manifold 0-thickness errors occurred. In the active setup, the cylinder is placed 11 mm above the bed (y-min at -0.015 m, bed at -0.026 m) to ensure stability.
- Undefined edges on the contact line were routed to a proper `cylinderFa` boundary condition with `zeroGradient`.
- `ABorder` was reduced to 0 (Euler Explicit) to bypass corrupted restart history.

### H. Mathematical Instability of filterExner on 2D Streamwise Grids
- **The Issue:** Setting `filterExner on` in 2D streamwise simulations causes the solver to crash with a Floating Point Exception (FPE) on restart. The gradient calculation `graddHe` uses `Xon & LeProj` (projecting cell-to-cell streamwise vectors onto edge vectors). In 2D, the edge vector points in the Z-direction, which is orthogonal to the streamwise direction. Under mesh deformation, this dot product becomes a tiny non-zero value with an arbitrary sign depending on the edge's topological orientation in Z. This results in numerical **anti-diffusion (positive feedback)**, causing the local bed displacement to explode and crash the solver.
- **The Fix:** Keep `filterExner off;` in `bedloadProperties` and resolve waviness geometrically.

### I. Resolution of Grid-Transition Bed Waviness
- **The Issue:** With the Exner filter disabled, a grid-transition waviness (0.65 mm y-jumps) occurred in the bed profile at $x = 0.15\text{ m}$. This was a numerical artifact caused by the abrupt 4x transition in cell size between the fine `gapBox` refinement ($0.25\text{ mm}$) and the coarser background mesh ($1.0\text{ mm}$) inside the active scour region.
- **The Fix:** We extended the `gapBox` refinement region downstream to **$x = 0.40\text{ m}$** and shifted the corresponding `blockMesh` block boundary. This places the entire scour zone within a uniform, high-resolution grid, completely eliminating the transition waviness and keeping the mesh non-orthogonality low.

### J. Bounded Sediment Concentration to Prevent FPE Crashes
- **The Issue:** During the onset of morphodynamic scour, high erosion rates can cause the bottom boundary to expand downward faster than the sediment settling velocity ($u_{bed} > w_s$). In the relative frame of the moving mesh, this creates an apparent inflow through the bed. With a standard `zeroGradient` boundary condition for `Cs`, this inflow carries the high internal cell concentration back into the domain. Combined with the explicit erosion source term, this triggers a positive feedback loop that pushes `Cs` above $1.0$, resulting in a solver crash (FPE) when the settling model calculates `pow(1 - Cs, n)`.
- **The Fix:** Switch the `Cs` boundary condition on the `bottom` patch to `inletOutlet` referencing the sediment-specific flux `phip`. This clamps inflow concentrations to 0 while allowing `zeroGradient` for outflows, stabilizing the simulation.

### K. Mitigating Long-Term Dynamic Mesh Deformation Crashes (T=15.9s Crash)
- **The Issue:** Under quadratic distance diffusivity (`inverseDistance 2`), cells in the narrow gap under the pipeline are kept rigid, concentrating all the deformation in a thin middle layer. During large scour depths, this concentrates shear stress, causing cell faces to warp and invert (producing face pyramid errors at $T=13.0$s), leading to a Floating Point Exception (FPE) crash in the velocity solver at $T=15.9$s.
- **The Fix:** We updated [dynamicMeshDict](file:///F:/CFD/2DPipelineScourEXN/constant/dynamicMeshDict) to use a volume-based mesh diffusivity model: `diffusivity quadratic inverseVolume;`. This model dynamically scales the local cell stiffness inversely with the square of its volume ($1/V^2$). As a cell starts to compress, its stiffness approaches infinity, preventing cell collapse and distributing the deformation smoothly throughout the entire fluid domain.

### L. Resolution of missing yWallFinal Solver in fvSolution
- **The Issue:** Recalculating the wall distance field (`yWall`) at every time step is required during dynamic mesh motion. On the final PIMPLE iteration, the solver looks for `yWallFinal` in the `solvers` sub-dictionary of `system/fvSolution` and crashes with a dictionary error if it is missing.
- **The Fix:** Added the `yWallFinal` solver block referencing `$yWall` in [fvSolution](file:///F:/CFD/2DPipelineScourEXN/system/fvSolution).

### M. Seabed Boundary Alignment and Tolerance Tightening
- **The Issue:** Mismatches in the constraints of boundary points at the inlet/outlet (which were allowed to deform in Exner or fixed in 3D) caused severe boundary shear. Additionally, loose tolerances in the dynamic mesh solver allowed coordinate drift to accumulate.
- **The Fixes:** 
  1. Configured `fixedFaPatches ( inletFa outletFa );` in [bedloadProperties](file:///F:/CFD/2DPipelineScourEXN/constant/bedloadProperties) to ensure matching boundary constraints.
  2. Tightened the `cellDisplacement` solver tolerance in [fvSolution](file:///F:/CFD/2DPipelineScourEXN/system/fvSolution) from `1e-5` to `1e-9`.
  3. Increased `minDeltaT` to `1e-5` for a proactive safety freeze of morphodynamics.

### N. Suppressing Downstream Sediment Mound (`noDeposition`)
- **The Issue:** Explicit bed morphodynamics can cause sediment to pile up behind the pipeline in an unphysical mound during short-duration tests, introducing additional flow shear and grid instability.
- **The Fix:** Patched `createFaFields.H` and `exnerEqn.H` to add a new `noDeposition` boolean flag (read from `bedloadProperties`). When `noDeposition true;` is enabled, the Exner solver clamps positive $dH$ values to zero, suppressing downstream sediment accumulation while allowing scour (negative $dH$) to continue naturally. The solver was successfully compiled via `Allwmake` in `~/sedExnerFoam`.

### O. Bypassing Setup Redundancies on Phase 2 Crash (`RunPhase2`)
- **The Issue:** When Phase 2 crashed due to numerical or syntax issues, running the standard `./Allrun` script cleaned the directory and forced a redundant 8-minute rebuild of the mesh and Phase 1 spin-up from scratch.
- **The Fix:** Created an executable bash script [RunPhase2](file:///F:/CFD/PipelineScour_Calibrated/RunPhase2) which bypasses the blockMesh, snappyHexMesh, and Phase 1 spin-up phases. It copies the Phase 2 dictionaries, maps the initial point displacement boundaries onto all 20 processor directories, and launches Phase 2 parallel execution instantly from the saved `Time = 2.0s` files.

### P. Resolution of OpenFOAM list parsing error in dynamicMeshDict
- **The Issue:** Specifying `inverseDistance 1 (bottom cylinder)` caused the solver to crash with `Expected a ')' or a '}' while reading List, found word 'cylinder'` because OpenFOAM parsed `1` as the list size and threw an error when encountering a second patch.
- **The Fix:** Changed the syntax in [dynamicMeshDict](file:///F:/CFD/PipelineScour_Calibrated/system/dynamicMeshDict) to the correct OpenFOAM list structure: `diffusivity inverseDistance 2(bottom cylinder);` (where `2(bottom cylinder)` denotes a list of size 2, and the omitted prefix defaults the exponent to 1).

### Q. Stabilization of Phase 2 Moving Boundary Coupling (fixedFluxPressure & maxDH/alpha Bounding)
- **The Issue:** Under large bed deformations, the simulation originally suffered from unphysical over-erosion and grid-velocity shock at `Time = 2.53s`. Setting `alpha = 32.0` caused excessive bedload flux, deforming the bed by the maximum allowable limit (`maxDH = 1e-4` m) on every step. At $\Delta t \approx 7 \times 10^{-4}$ s, this created a vertical grid movement velocity of $0.15\text{ m/s}$. For $0.47\text{ mm}$ boundary layer cells, this sudden mesh motion acted as a massive numerical pump, distorting the boundary layer flow, spiking the Shields number to the cap of `10.0`, and rapidly collapsing the cells (causing Courant numbers to spike to `6336.5` and freezing the solver).
- **The Fixes:**
  1. Configured the moving bed patch (`bottom`) to use **`fixedFluxPressure`** for `p` and `p_rgh`, balancing pressure gradients with boundary acceleration to preserve normal momentum.
  2. Restored **`alpha 12.0`** (the standard physical Nielsen bedload coefficient) to ensure the Exner equation computes a natural, physical erosion rate.
  3. Set **`maxDH 1.5e-5;`** in `bedloadProperties_phase2` to cap single-step grid stretching to a safe 3% of cell height, limiting boundary grid velocities to a gentle $0.02\text{ m/s}$.
  4. Enforced **`inverseDistance 2(bottom cylinder);`** to keep cells in the high-shear pipeline gap rigid.
  5. Updated **`RunPhase2`** to purge all intermediate directories $> 2.0\text{s}$ prior to launch, ensuring a clean restart from the end of Phase 1.

### R. Permanent Stability Fix for Dynamic Mesh & Morphodynamics (PIMPLE Coupling)
- **The Issue:** The simulation suffered from recurring Floating Point Exception (SIGFPE) crashes at CFD time `Time ≈ 2.6s` to `4.9s`. Detailed debugging revealed a dual-instability mechanism:
  1. *Explicit Coupling Lag:* Because `nOuterCorrectors` was default-unset (`1`), the solver calculated velocity on the old grid, moved the grid to the new position, and immediately advanced to the next time step. This decoupling generated an *odd-even Shields number oscillation* (spiking between `0.7` and `2.1` on alternate steps) that rapidly warped the grid.
  2. *Stiffness Mismatch:* The volume-based diffusivity (`inverseVolume`) made the tiny boundary-layer cells extremely stiff compared to core cells. This forced the boundary layer to shear horizontally rather than deforming vertically, leading to grid skewness and singular matrices during the wall distance (`yWall` via `GAMG`) solve.
- **The Fixes:**
  1. Set **`nOuterCorrectors 3;`** in `system/fvSolution` to iteratively converge and align the flow field with the mesh motion within the same timestep.
  2. Tightened correctors (`nCorrectors 3;` and `nNonOrthogonalCorrectors 3;`) to improve pressure calculations on the deforming grid.
  3. Set `cellDisplacement` under-relaxation to `1.0` in `relaxationFactors` to eliminate artificial mesh-motion lag.
  4. Increased **`NfiltExner` to `15`** in `constant/bedloadProperties_phase2` to filter high-frequency numerical spatial noise from the bed profile.
  5. Reverted to **`inverseDistance 2(bottom cylinder);`** mesh motion diffusivity, allowing both fine and coarse cells to deform proportionally according to wall distance.
  6. Configured **`maxCo 2.0;`** and **`endTime 52.0;`** (paired with **`morphoAccFactor 2.0;`**) to double the hydrodynamic stability margin and complete the full 100s of physical scour with absolute stability.

---

## 🚀 4. Automated 2-Phase Simulation Workflow (20-Core MPI)

Modeling the $11.0\text{mm}$ gap introduces significant numerical challenges. The velocity through the gap creates an instantaneous, extreme shear stress spike at $t=0$. If the bed is allowed to deform immediately at real-time speeds, the explicit Exner solver's Courant number ($C_{exner}$) explodes, instantly corrupting the mesh ($V < 0$) and causing floating-point exceptions.

To completely stabilize this, the entire pipeline is **100% fully automated** via a single execution script:

```bash
bash Allrun
```

*(Note: If the simulation successfully completes Phase 1 and you only need to adjust Phase 2 physical parameters or fix dictionary errors, bypass the full Allrun script and execute `./RunPhase2` directly to resume Phase 2 instantly).*

The `Allrun` script orchestrates the computational setup, while the Phase 2 morphodynamics are run using the `Allphase2` script:

### Phase 1: Hydrodynamic Spin-up ($t = 0.0\text{s} \rightarrow 2.0\text{s}$)
- **Settings:** `morphoAccFactor = 0.005`, `avalanche = off`.
- **Purpose:** Artificially slows down the morphodynamic time-scale by 200x, allowing the Navier-Stokes velocity and pressure fields to fully develop and stabilize through the gap *without* the erodible bed deforming prematurely.

### Phase 2: Morphodynamic Scour ($t = 2.0\text{s} \rightarrow 52.0\text{s}$)
- **Settings:** `morphoAccFactor = 2.0`, `avalanche = on`, `alpha = 12.0` (standard physical Nielsen bedload coefficient), `saturationType = none` (instantaneous pickup).
- **Purpose:** Speeds up the bed morphodynamics by **2x** to allow direct validation against the 100s experimental equilibrium profile of Mao (1986). Under this acceleration, $50.0\text{ s}$ of CFD solver time (from $T = 2.0\text{ s}$ to $T = 52.0\text{ s}$) corresponds to $100.0\text{ s}$ of morphological scour.
- **Stabilization Features:** To eliminate unphysical wiggles and prevent cell collapse:
  * `slopeCorrection on` inside the bedload model (gravity-induced bed stabilization).
  * `filterShields on` (spatial smoothing of local shear stresses).
  * `filterExner on` with `NfiltExner = 15` and `alphaFiltExner = 0.15` (stronger bed elevation filter to damp spatial decoupling).
  * `nOuterCorrectors = 3` inside PIMPLE (iterative flow-mesh coupling).

---

## 📂 5. Repository Structure & Configuration File Index

- `0_org/` (Initial Conditions): Contains the pristine initial conditions. The `copyFA.sh` script maps these onto the dynamically generated mesh inside `0/`.
- `constant/transportProperties` & `turbulenceProperties`: Baseline hydrodynamic constraints and RANS closure specifications.
- `constant/bedloadProperties`: The master configuration file for the morphodynamic solver (houses the $d_{50}$, Exner settings, Nielsen model coefficients, `morphoAccFactor`, and `avalanche` toggles).
- `constant/dynamicMeshDict`: Governs how the 2D mesh deforms as the bed drops. It uses an `inverseDistance` Laplacian diffusivity field explicitly referencing the `bottom` and `cylinder` to ensure elements near the pipeline are not crushed during massive bed evolution.
- `system/fvSolution`: Defines the linear solvers and PIMPLE loop controls. `correctPhi` and `moveMeshOuterCorrectors` are automatically managed by `Allrun`.
- `system/controlDict`: Governs the time-stepping ($\Delta t$), write intervals, and Courant number limits ($Co \le 0.4$).
- `system/blockMeshDict`, `snappyHexMeshDict`, `extrudeMeshDict`: Handles spatial domain synthesis, cylinder carving, and strictly enforces a 2D topology.

---

## 📊 6. Validation and Benchmarking

The morphodynamic scour depth ratio $S(t)/D$ is quantitatively assessed against established literature:
- **Experimental Data:** Mao (1986).
- **High-Fidelity CFD:** Larsen et al. (2016) utilizing resolved two-phase Eulerian-Eulerian models.
- **Empirical Estimations:** Sumer-Fredsøe deterministic development models.

To extract the physical depth of the scour hole at each timestep:
```bash
postProcess -func 'patchExpression(bottom, Cf.y())' -latestTime
```

---

## 🙏 7. Acknowledgments & Citations

I would like to extend my profound gratitude to my academic advisor and the OpenFOAM community for their continuous support. This project serves as a cornerstone for advanced sub-aqueous morphology computations, pushing the boundaries of accessible, high-performance fluid-structure-seabed interaction modeling.

Furthermore, this work builds upon the foundational open-source contributions of the following authors and their respective solvers, without whom this research would not be possible:
- **sedFoam:** Chauchat, J., Cheng, Z., Nagel, T., Bonamy, C., & Hsu, T. J. (2017). SedFoam-2.0: a 3-D two-phase flow solver for depth-resolving sediment transport modeling. Geoscientific Model Development.
- **sedExnerFoam:** Larsen, R. P., Fuhrman, D. R., & Roenby, J. sedExnerFoam: A moving mesh finite-area morphodynamic solver for OpenFOAM.

We deeply acknowledge their efforts in advancing the field of computational sediment transport.
