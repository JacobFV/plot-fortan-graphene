#!/usr/bin/env python3
"""
Python version of plotBands.f90 - Electronic band structure calculation
for twisted bilayer graphene.

Compile equivalent:
python plot_bands.py

Original Fortran compile commands:
gfortran -lm bilayer LapackRoutines.f90 Hamiltonian.f90 plotBands.f90 -llapack -lblas
ifort -qopenmp -o test LapackRoutines.f90 Hamiltonian.f90 plotBands.f90 -qmkl -lpthread -lm
"""

import numpy as np
import os
import sys
from joblib import Parallel, delayed
from hamiltonian import (
    NUMK, NTHETA, NDIM, NUMB, NCB, NRELAX, PI, VKD, A1, A2,
    wigner_seitz, complex_h0, relaxation, sample_points_bands, plot_bands_gamma
)
from lapack_routines import diagonalize


def calculate_bands_for_k(vk, coord, tn, il, iu):
    """
    Calculate bands for a single k-point.
    
    Parameters:
    -----------
    vk : np.ndarray
        k-vector
    coord : np.ndarray
        Atomic coordinates
    tn : np.ndarray
        Lattice vectors
    il, iu : int
        Band index range
        
    Returns:
    --------
    evals : np.ndarray
        Eigenvalues for this k-point
    """
    # Build Hamiltonian matrix
    zh = complex_h0(coord, vk, tn)
    
    # Allocate eigenvalues array
    evals = np.zeros(NDIM)
    
    # Diagonalize
    eigenvals = diagonalize(zh, evals, 'N', il, iu)[0]
    
    return eigenvals


def main():
    """Main calculation routine."""
    
    print("Starting band structure calculation...")
    print(f"NTHETA = {NTHETA}")
    print(f"NUMK = {NUMK}")
    print(f"NDIM = {NDIM}")
    
    # Setup geometry - Brillouin zone
    a_moire = 3.0 * NTHETA**2 + 3.0 * NTHETA + 1.0
    
    # Reciprocal lattice vectors of the superlattice
    vq1 = VKD / a_moire * ((3*NTHETA + 1) * A1 + A2)
    vq12 = VKD / a_moire * ((3*NTHETA + 2) * A2 - A1)
    
    # Angle of rotation
    cs = 1.0 - 1.0 / (2.0 * a_moire)
    ang = 0.5 * np.arccos(cs)
    ang_mat = np.array([[np.cos(ang), -np.sin(ang)], 
                        [np.sin(ang), np.cos(ang)]])
    sn = np.sqrt(1.0 - cs**2)
    
    # Rotate reciprocal lattice vectors
    vq1 = ang_mat @ vq1
    vq12 = ang_mat @ vq12
    
    # Setup G vectors
    gn = np.zeros((8, 2))
    gn[0, :] = vq1
    gn[1, :] = vq12
    gn[2, :] = vq12 - vq1
    gn[3, :] = -gn[0, :]
    gn[4, :] = -gn[1, :]
    gn[5, :] = -gn[2, :]
    gn[6, :] = gn[0, :]
    gn[7, :] = gn[1, :]
    
    print(f'Twist Angle: {2*ang*180/PI:.6f} degrees')
    
    # Setup geometry - real space
    # Moire lattice parameters
    t1 = NTHETA * A1 + (NTHETA + 1) * A2
    t2 = -(NTHETA + 1) * A1 + (2*NTHETA + 1) * A2
    t3 = t2 - t1
    
    tn = np.array([t1, t2, t3, -t1, -t2, -t3]).T  # Shape (2, 6)
    
    # Generate Wigner-Seitz coordinates
    coord = wigner_seitz(t1, t2, t3, cs, sn)
    
    # Rotate cell to symmetrize around x-axis
    coord_rotated = np.zeros_like(coord)
    for n in range(NDIM):
        coord_rotated[n, :2] = ang_mat @ coord[n, :2]
        coord_rotated[n, 2] = coord[n, 2]
    
    coord = coord_rotated
    
    # Rotate lattice vectors
    t1 = ang_mat @ t1
    t2 = ang_mat @ t2
    t3 = t2 - t1
    
    tn = np.array([t1, t2, t3, -t1, -t2, -t3]).T
    
    # Apply relaxation if requested
    if NRELAX == 1:
        print("Applying lattice relaxation...")
        relaxation(coord, gn)
    
    # Setup band calculation parameters
    ndim_ev = NUMB
    il = (NDIM - 2*ndim_ev + 2*NCB) // 2 + 1
    iu = (NDIM + 2*NCB) // 2
    n_np = (NDIM // 2 - il) + 1
    
    ncount = NUMK // 3 + NUMK // 2 + NUMK // 6 + 1
    
    print(f"Band range: {il} to {iu}")
    print(f"Number of k-points: {ncount}")
    
    # Sample k-points for band structure
    points_bz, points_back = sample_points_bands(NUMK, ncount, vq1, vq12)
    
    # Initialize bands array
    bands_tb = np.zeros((ndim_ev, ncount))
    
    # Prepare k-points for parallel calculation
    k_points = []
    for icount in range(ncount):
        vk = (points_bz[icount, 0] * vq1 + points_bz[icount, 1] * vq12) / NUMK
        k_points.append((vk, icount))
    
    print("Starting parallel band calculation...")
    
    # Parallel calculation of bands
    def calc_single_k(k_data):
        vk, icount = k_data
        eigenvals = calculate_bands_for_k(vk, coord, tn[:, :3], il, iu)
        return icount, eigenvals
    
    # Use joblib for parallel processing (replacement for OpenMP)
    results = Parallel(n_jobs=-1, verbose=1)(
        delayed(calc_single_k)(k_data) for k_data in k_points
    )
    
    # Collect results
    for icount, eigenvals in results:
        bands_tb[:, icount] = eigenvals[:ndim_ev]
    
    # Calculate Fermi level (neutrality point)
    fmu = (bands_tb[n_np - 1, points_back[NUMK//3, NUMK//3]] + 
           bands_tb[n_np, points_back[NUMK//3, NUMK//3]]) / 2.0
    
    print(f'Neutrality point (Fermi level): {fmu:.6f}')
    
    # Generate output filename
    parameters = f'-I{NTHETA}-numk{NUMK}-relax{NRELAX}.dat'
    filename = 'Bands' + parameters
    
    print(f"Writing results to {filename}")
    
    # Plot bands
    plot_bands_gamma(bands_tb, fmu, points_back, ndim_ev, NUMK, ncount, vq1, vq12, filename)
    
    print("Calculation completed successfully!")


if __name__ == "__main__":
    main() 