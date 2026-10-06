# pimpleCPMFoam

`pimpleCPMFoam` is a custom OpenFOAM solver for transient incompressible PIMPLE flow with bipolar ion transport and electrostatic coupling. It solves the gas flow, electric potential, positive/negative ion number densities, and the resulting electric body force in the momentum equation. It is especially useful for CPM/electrode cases where ion deposition changes the charged plate potential and where the ion Courant number should control the time step.

This solver has been checked only with OpenFOAM-v2406 through OpenFOAM-v2506.

## Features

- Based on `pimpleFoam` for transient incompressible Newtonian flow with optional dynamic mesh support.
- Adds Poisson's equation for `phiE`, electric field reconstruction `E = -grad(phiE)`, and drift-diffusion-recombination transport for `nP` and `nN`.
- Supports ion deposition accounting on patches named `CPM` and `walls`; the `CPM` patch can update a fixed-value `phiE` boundary according to deposited current and `Ccap`.

## Compilation

Load an OpenFOAM-v2406, v2412, or v2506 environment first:

```bash
source /path/to/OpenFOAM-v2506/etc/bashrc
```

Then build from this directory:

```bash
wmake
```

The executable is written to:

```text
$FOAM_USER_APPBIN/pimpleCPMFoam
```

To clean generated build files:

```bash
wclean
```

## Required Case Setup

Prepare a normal incompressible PIMPLE case and add the electrostatic/ion fields.

Required fields in the initial time directory:

- `U`
- `p`
- `phiE`
- `nP`
- `nN`

Optional fields that the solver can read if already present, otherwise it calculates and writes them:

- `E`
- `rho`
- `phiEF`

Add these entries to `constant/physicalProperties`:

```text
T       T       [0 0 0 1 0 0 0] 300;
rho0    rho0    [1 -3 0 0 0 0 0] 1.2;
nu      nu      [0 2 -1 0 0 0 0] 1.5e-05;
Ccap    Ccap    [-1 -2 4 0 0 2 0] 1e-12;
muP     muP     [-1 0 2 0 0 1 0] 1.5e-04;
muN     muN     [-1 0 2 0 0 1 0] 1.5e-04;
DP      DP      [0 2 -1 0 0 0 0] 1e-05;
DN      DN      [0 2 -1 0 0 0 0] 1e-05;
beta    beta    [0 3 -1 0 0 0 0] 1e-12;
```

Use dimensions and values appropriate for your case.

## Usage

Run it like a standard OpenFOAM solver:

```bash
pimpleCPMFoam
```

For parallel cases:

```bash
decomposePar
mpirun -np <N> pimpleCPMFoam -parallel
reconstructPar
```

The solver looks for patches named `CPM` and `walls`. If either patch is missing, it prints a warning and skips the corresponding deposition accounting. For `CPM` potential updates, set `phiE` on that patch as a `fixedValue` boundary condition.

The time-step selection can be controlled in `system/controlDict`:

```text
adjustTimeStep  yes;
adjustTimeType  ionBased;  // or phi
maxCo           0.5;
```

Use `ionBased` when the drift velocity from the electric field is the restrictive transport speed. Use `phi` when the ordinary flow Courant number should control the time step.

## Notes

- The electric body force enters the momentum equation as `-e*(nP-nN)*E/rho0`.
- Ion densities are clipped to non-negative values after transport solves.
- The deposition model subtracts only outward positive deposited amounts from adjacent cells, which helps avoid negative ion number densities near absorbing patches.
