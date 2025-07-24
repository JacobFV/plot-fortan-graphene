# Twisted Bilayer Graphene Band Structure Calculator

This is a Python conversion of the original Fortran 90 code for calculating electronic band structures of twisted bilayer graphene (TBG). The code implements a tight-binding model with lattice relaxation to compute the band structure along high-symmetry paths in the Brillouin zone.

## Original Fortran Code

The original code consisted of three main files:
- `LapackRoutines.f90` - LAPACK wrapper routines for matrix diagonalization
- `Hamiltonian.f90` - Main physics module with Hamiltonian construction and utilities
- `plotBands.f90` - Main program for band structure calculation

## Python Conversion

The Python version maintains the same algorithmic structure with the following files:
- `lapack_routines.py` - Python equivalents using scipy.linalg
- `hamiltonian.py` - Main physics module with all functions and constants
- `plot_bands.py` - Main script for running calculations

## Installation

0. (Optional) Create a Python virtual env: 
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

1. Install Python dependencies:
```bash
pip install -r requirements.txt
```

2. Run the calculation:
```bash
python plot_bands.py
```

## Dependencies

- **numpy**: For numerical arrays and mathematical operations
- **scipy**: For LAPACK routines (matrix diagonalization)
- **joblib**: For parallel processing (replaces OpenMP from Fortran)
- **matplotlib**: For optional plotting utilities

## Key Features

### Physics Parameters
- **NTHETA**: Commensurate twist angle parameter (default: 5)
- **NUMK**: Number of k-points along each direction (default: 30)
- **NRELAX**: Enable lattice relaxation (1) or not (0)
- **PRESSURE**: Hydrostatic pressure parameter for tuning

### Computational Features
- Automatic parallelization using joblib (replaces OpenMP)
- Partial matrix diagonalization for efficiency
- Support for different twist angles with pre-computed relaxation parameters
- Output compatible with original Fortran format

## Usage

The main script `plot_bands.py` will:
1. Set up the moiré superlattice geometry
2. Generate atomic coordinates using Wigner-Seitz construction
3. Apply lattice relaxation (if enabled)
4. Sample k-points along Γ-K-M-Γ path
5. Calculate band structure using parallel processing
6. Output results to a data file

### Customization

To modify parameters, edit the constants in `hamiltonian.py`:

```python
NTHETA = 5      # Twist angle parameter
NUMK = 30       # K-point sampling
NRELAX = 1      # Enable relaxation
PRESSURE = 1.0  # Hydrostatic pressure
```

### Output

The program generates a file named `Bands-I{NTHETA}-numk{NUMK}-relax{NRELAX}.dat` containing:
- Column 1: k-point coordinate along path
- Columns 2+: Band energies relative to neutrality point

## Performance Notes

- The Python version uses joblib for parallelization, which may have different performance characteristics than OpenMP
- Matrix operations use optimized BLAS/LAPACK through scipy
- Memory usage is similar to the original Fortran code
- For large systems, consider adjusting the number of parallel jobs in `plot_bands.py`

## Differences from Original

1. **Parallelization**: Uses joblib instead of OpenMP
2. **Matrix operations**: scipy.linalg instead of direct LAPACK calls
3. **I/O**: Python file handling instead of Fortran
4. **Memory management**: Automatic Python garbage collection
5. **Indexing**: 0-based indexing (Python) vs 1-based (Fortran)

## Validation

The Python version should produce numerically identical results to the original Fortran code when using the same parameters. Small differences may arise due to:
- Different LAPACK implementations
- Compiler optimizations
- Floating-point precision handling

## Extending the Code

To add new twist angles or modify the physics:

1. **New twist angles**: Add relaxation coefficients to the `relaxation()` function in `hamiltonian.py`
2. **Modified hopping**: Edit the `ft()` and `ftperp()` functions
3. **Different k-paths**: Modify `sample_points_bands()` function
4. **Additional observables**: Extend the diagonalization to compute eigenvectors

## Troubleshooting

- **Memory errors**: Reduce NUMK or NTHETA for smaller systems
- **Slow performance**: Check joblib installation and CPU core count
- **Import errors**: Ensure all dependencies are installed with correct versions
- **Numerical issues**: Verify LAPACK installation through scipy

For questions or issues with the conversion, please check that the original Fortran parameters are correctly mapped to the Python constants. 