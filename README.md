## Purpose
 
Turbulence plays a key role in most engineering and science applications, hence the need for reliable mathematical tools that can capture the main aspects of turbulent flows without an excessive demand on computational power so a broad community can assess their suitability in a number of cases. The Spalart-Allmaras model has proved to be reliable in various situations, especially in aerodynamic flows and, due
to its single transport equation, the hardware requirements are not prohibitively, nonetheless, some authors have pointed out some problems at the edge of boundary layers and wakes where the transported variable may become negative on under-resolved grids. The negative Spalart-Allmaras model is proposed to treat adequately the aforementioned situation by introducing a new equation that is solved instead of the standard one in those regions where the turbulent variable $`\tilde{\nu} < 0`$.

Implementation of the Negative Spalart-Allmaras model is carried out and successfully tested on five cases on OpenFOAM v2012 release (it is expected to work without problem in v2206):

- **2D Incompressible boundary layer:** Reynolds per unit length $`Re = 5 \times 10^{6}`$.

- **2D Incompressible NACA0012 airfoil:** Reynolds based on chord length $`Re_c = 3.95 \times 10^{6}`$ and various angles of attack $`\alpha=0^\circ,4.09^\circ,10.08^\circ,15.14^\circ`$.

- **2D Compressible NACA0012 airfoil:** Reynolds based on chord length $`Re_c = 9 \times 10^{6}`$ and angle of attack $`\alpha=1.86^\circ`$.

- **2D Backward-facing step:** Reynolds based on step height $`Re_H=36000`$.

- **3D Bump-in-channel:** Reynolds per unit length $`Re = 3 \times 10^{6}`$.

## Compilation

The following process assumes the model will be saved in the user directory ```$WM_PROJECT_USER_DIR/```.

- Move the ```negSpalartAllmaras.C``` and ```negSpalartAllmaras.H``` files (location: ```Source_Code/turbulenceModels/RAS/negSpalartAllmaras/```) to ```$WM_PROJECT_USER_DIR/src/TurbulenceModels/turbulenceModels/RAS/negSpalartAllmaras/``` on your local system.

- Move the ```compressible``` and ```incompressible``` folders (location: ```Source_Code```) to ```$WM_PROJECT_USER_DIR/src/TurbulenceModels/``` on your local system.

- Go to the ```incompressible``` folder and create the ```lnInclude``` folder, then compile:

```bash
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/incompressible
wmakeLnInclude -u ../turbulenceModels
wmake 
```
- Go to the ```compressible``` folder and create the ```lnInclude``` folder, then compile:

```bash
cd $WM_PROJECT_USER_DIR/src/TurbulenceModels/compressible
wmakeLnInclude -u ../turbulenceModels
wmake 
```

Compilation is complete for both incompressible and compressible solvers. Remember to add the following lines at the end of ```controlDict``` file according to the case:

```cpp
// Incompressible
libs ( " libmyIncompressibleTurbulenceModels.so");
// Compressible
libs ( " libmyCompressibleTurbulenceModels.so");
```

In the ```turbulenceProperties``` file set the new model as follows (for compressible cases it may be important to set Ct3 = 0),

```cpp
RAS
{
    RASModel        negSpalartAllmaras;

    turbulence      on;
 // Ct3             0;      (Compressible cases) 
    printCoeffs     on;
}
```