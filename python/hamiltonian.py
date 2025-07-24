import numpy as np
from typing import Tuple, Optional
from lapack_routines import diagonalize
import multiprocessing as mp
from joblib import Parallel, delayed

# Constants and parameters (equivalent to Fortran parameters)
# Precision
DP = np.float64

# Simulation parameters
NUMK = 30  # number of k-points, use multiples of 6
NTHETA = 5  # commensurate twist angle 
NDIM = 4 * (NTHETA**2 + (NTHETA+1)**2 + (NTHETA+1)*NTHETA)  # dimension of Ham
NUMB = 8  # total bands to be plot (maximal value is ndim)
NCB = 4  # conduction bands to be plot (maximal value is ndim/2)
NRELAX = 1  # include lattice relaxation by nrelax=1 (0 otherwise)

# Hamiltonian parameters
FPHASE = 0.0
PRESSURE = 1.0  # tune "twist angle" by hydrostatic pressure, e.g., 1.096
DELTA_R = 0.1
A0 = 1.0 / np.sqrt(3.0)
R0 = 0.184
D0 = 1.35772 / PRESSURE
TPERP = 1.0
MASS1 = 0.0  # in eV
TZ = 1.35772 / PRESSURE
A1 = np.array([0.5, np.sqrt(3.0) * 0.5])
A2 = np.array([-0.5, np.sqrt(3.0) * 0.5])

# General parameters
PI = np.pi
SQ3 = np.sqrt(3.0)
VKD = 4.0 * PI / 3.0


def wigner_seitz(t1: np.ndarray, t2: np.ndarray, t3: np.ndarray, 
                cs: float, sn: float) -> np.ndarray:
    """
    Generate Wigner-Seitz coordinates for the supercell.
    
    Parameters:
    -----------
    t1, t2, t3 : np.ndarray
        Lattice vectors
    cs, sn : float
        Cosine and sine of rotation angle
        
    Returns:
    --------
    coord : np.ndarray
        Coordinates array of shape (ndim, 3)
    """
    nrad = 3 * NTHETA
    coord = np.zeros((NDIM, 3))
    
    # Initialize counters
    nind = 0
    na1 = na2 = nb1 = nb2 = 0
    
    r_max = 3.0 * NTHETA**2 + 3.0 * NTHETA + 1.0
    
    # Layer 1, site A
    for n1 in range(-nrad, nrad + 1):
        for n2 in range(-nrad, nrad + 1):
            r_temp1 = abs(n1 * (3.0 * NTHETA + 1) + n2 * (3.0 * NTHETA + 2)) + 1e-6
            r_temp2 = abs(n1 - n2 * (3.0 * NTHETA + 1)) + 1e-6
            r_temp3 = abs(n1 * (3.0 * NTHETA + 2) + n2) + 1e-6
            
            if r_temp1 < r_max and r_temp2 < r_max and r_temp3 < r_max:
                na1 += 1
                coord[nind, 0] = n1 * A1[0] + n2 * A2[0]
                coord[nind, 1] = n1 * A1[1] + n2 * A2[1]
                coord[nind, 2] = 0.0
                nind += 1
    
    # Layer 1, site B
    for n1 in range(-nrad, nrad + 1):
        for n2 in range(-nrad, nrad + 1):
            bn1 = n1 + 1.0/3.0
            bn2 = n2 + 1.0/3.0
            
            r_temp1 = abs(bn1 * (3.0 * NTHETA + 1) + bn2 * (3.0 * NTHETA + 2)) + 1e-6
            r_temp2 = abs(bn1 - bn2 * (3.0 * NTHETA + 1)) + 1e-6
            r_temp3 = abs(bn1 * (3.0 * NTHETA + 2) + bn2) + 1e-6
            
            if r_temp1 < r_max and r_temp2 < r_max and r_temp3 < r_max:
                nb1 += 1
                coord[nind, 0] = bn1 * A1[0] + bn2 * A2[0]
                coord[nind, 1] = bn1 * A1[1] + bn2 * A2[1]
                coord[nind, 2] = 0.0
                nind += 1
    
    # Include B atom at zone boundary (Layer 1)
    coord[nind, 0] = (t2[0] + t3[0]) / 3.0
    coord[nind, 1] = (t2[1] + t3[1]) / 3.0
    coord[nind, 2] = 0.0
    nind += 1
    
    # Layer 2, site A
    for n1 in range(-nrad, nrad + 1):
        for n2 in range(-nrad, nrad + 1):
            an1 = n1 * (cs - sn/SQ3) - 2*n2 * sn/SQ3
            an2 = n2 * (cs + sn/SQ3) + 2*n1 * sn/SQ3
            
            r_temp1 = abs(an1 * (3.0 * NTHETA + 1) + an2 * (3.0 * NTHETA + 2)) + 1e-6
            r_temp2 = abs(an1 - an2 * (3.0 * NTHETA + 1)) + 1e-6
            r_temp3 = abs(an1 * (3.0 * NTHETA + 2) + an2) + 1e-6
            
            if r_temp1 < r_max and r_temp2 < r_max and r_temp3 < r_max:
                na2 += 1
                coord[nind, 0] = an1 * A1[0] + an2 * A2[0]
                coord[nind, 1] = an1 * A1[1] + an2 * A2[1]
                coord[nind, 2] = TZ
                nind += 1
    
    # Layer 2, site B
    for n1 in range(-nrad, nrad + 1):
        for n2 in range(-nrad, nrad + 1):
            bn1 = (n1 + 1.0/3.0) * (cs - sn/SQ3) - 2*(n2 + 1.0/3.0) * sn/SQ3
            bn2 = (n2 + 1.0/3.0) * (cs + sn/SQ3) + 2*(n1 + 1.0/3.0) * sn/SQ3
            
            r_temp1 = abs(bn1 * (3.0 * NTHETA + 1) + bn2 * (3.0 * NTHETA + 2)) + 1e-6
            r_temp2 = abs(bn1 - bn2 * (3.0 * NTHETA + 1)) + 1e-6
            r_temp3 = abs(bn1 * (3.0 * NTHETA + 2) + bn2) + 1e-6
            
            if r_temp1 < r_max and r_temp2 < r_max and r_temp3 < r_max:
                nb2 += 1
                coord[nind, 0] = bn1 * A1[0] + bn2 * A2[0]
                coord[nind, 1] = bn1 * A1[1] + bn2 * A2[1]
                coord[nind, 2] = TZ
                nind += 1
    
    # Include B atom at zone boundary (Layer 2)
    coord[nind, 0] = (t1[0] + t2[0]) / 3.0
    coord[nind, 1] = (t1[1] + t2[1]) / 3.0
    coord[nind, 2] = TZ
    nind += 1
    
    return coord


def ft(r: float) -> float:
    """In-plane hopping function."""
    if r > DELTA_R:
        return -2.7 * np.exp((A0 - r) / R0)
    else:
        return 0.0


def ftperp(r: float) -> float:
    """Inter-layer hopping function."""
    pd = (D0 / r)**2
    return (-2.7 * np.exp((A0 - r) / R0) * (1.0 - pd) + 
            0.48 * np.exp((D0 * PRESSURE - r) / R0) * pd) * TPERP


def complex_h0(coord: np.ndarray, vk: np.ndarray, tn: np.ndarray) -> np.ndarray:
    """
    Build the complex Hamiltonian matrix.
    
    Parameters:
    -----------
    coord : np.ndarray
        Atomic coordinates (ndim, 3)
    vk : np.ndarray
        k-vector (2,)
    tn : np.ndarray
        Lattice vectors (2, 3)
        
    Returns:
    --------
    zh : np.ndarray
        Complex Hamiltonian matrix (ndim, ndim)
    """
    zh = np.zeros((NDIM, NDIM), dtype=complex)
    
    for i in range(NDIM):
        for j in range(i + 1):  # Only lower triangle
            zphase = np.exp(-1j * np.dot(vk, coord[i, :2] - coord[j, :2]) * FPHASE)
            
            rij = np.linalg.norm(coord[i, :] - coord[j, :])
            coord_diff = (coord[i, 2] - coord[j, 2])**2
            
            # Intra-layer coupling
            if coord_diff < DELTA_R and rij > DELTA_R:
                zh[i, j] = ft(rij) * zphase
            
            # Inter-layer coupling
            if coord_diff > DELTA_R:
                zh[i, j] = ftperp(rij) * zphase
            
            # Coupling to adjacent Wigner-Seitz cells
            for nt in range(3):
                # Positive direction
                coord_adj = coord[j, :2] - tn[:, nt]
                rij_adj = np.linalg.norm(coord[i, :] - np.array([*coord_adj, coord[j, 2]]))
                zi_vk_tn = 1j * np.dot(vk, tn[:, nt])
                
                if coord_diff < DELTA_R:
                    zh[i, j] += ft(rij_adj) * np.exp(-zi_vk_tn) * zphase
                else:
                    zh[i, j] += ftperp(rij_adj) * np.exp(-zi_vk_tn) * zphase
                
                # Negative direction
                coord_adj = coord[j, :2] + tn[:, nt]
                rij_adj = np.linalg.norm(coord[i, :] - np.array([*coord_adj, coord[j, 2]]))
                zi_vk_tn = -1j * np.dot(vk, tn[:, nt])
                
                if coord_diff < DELTA_R:
                    zh[i, j] += ft(rij_adj) * np.exp(-zi_vk_tn) * zphase
                else:
                    zh[i, j] += ftperp(rij_adj) * np.exp(-zi_vk_tn) * zphase
    
    # Add mass terms
    for i in range(NDIM // 4):
        zh[i, i] += MASS1
        j = i + NDIM // 4
        zh[j, j] -= MASS1
    
    # Make Hermitian (fill upper triangle)
    for i in range(NDIM):
        for j in range(i + 1, NDIM):
            zh[i, j] = np.conj(zh[j, i])
    
    return zh


def relaxation(coord: np.ndarray, gn: np.ndarray):
    """
    Apply lattice relaxation to coordinates.
    
    Parameters:
    -----------
    coord : np.ndarray
        Atomic coordinates (ndim, 3) - modified in place
    gn : np.ndarray
        Reciprocal lattice vectors (8, 2)
    """
    # Relaxation coefficients (hardcoded for different ntheta values)
    uqx = np.zeros((5, 5))
    uqy = np.zeros((5, 5))
    
    # Load relaxation parameters based on NTHETA
    if NTHETA == 2:
        uqx[1, 0] = -0.0007312154977677281
        uqy[1, 0] = 0.0004714203951774816
        uqx[2, 0] = -1.154897989107238e-06
        uqy[2, 0] = 7.447553144970968e-07
        uqx[2, 1] = -1.586938089071483e-06
        uqy[2, 1] = 0.0
        # ... (more coefficients would be added here)
    elif NTHETA == 3:
        uqx[1, 0] = -0.001435896542049262
        uqy[1, 0] = 0.0008971722346957707
        uqx[2, 0] = -4.431771204029964e-06
        uqy[2, 0] = 2.770018809111695e-06
        uqx[2, 1] = -5.980337522539451e-06
        uqy[2, 1] = 0.0
        # ... (more coefficients)
    elif NTHETA == 5:
        uqx[1, 0] = -0.00352744756912305
        uqy[1, 0] = 0.002140730764691146
        uqx[2, 0] = -2.68010445377339e-05
        uqy[2, 0] = 1.627423402644026e-05
        uqx[2, 1] = -3.516363993145865e-05
        uqy[2, 1] = 0.0
        uqx[3, 0] = -2.729706294022432e-07
        uqy[3, 0] = 1.657616285518053e-07
        uqx[3, 1] = -4.604302998995516e-07
        uqy[3, 1] = 1.078115792758269e-07
        uqx[3, 2] = -4.604302998995516e-07
        uqy[3, 2] = -1.078115792758269e-07
        uqx[4, 0] = -3.147029295689003e-09
        uqy[4, 0] = 1.911153328301163e-09
        uqx[4, 1] = -6.629538953061466e-09
        uqy[4, 1] = 2.259053162518054e-09
        uqx[4, 2] = -7.46551582954267e-09
        uqy[4, 2] = 0.0
        uqx[4, 3] = -6.629538953061466e-09
        uqy[4, 3] = -2.259053162518054e-09
    # Add more NTHETA cases as needed...
    
    # Apply relaxation to layer 1
    for i in range(NDIM // 2):
        rx = coord[i, 0]
        ry = coord[i, 1]
        
        delta = np.zeros(2)
        
        # Calculate displacement
        for nn in range(0, 5, 2):  # 0, 2, 4 (equivalent to 1, 3, 5 in Fortran)
            for n in range(2, 6):  # 2, 3, 4, 5 (equivalent to Fortran 2 to 5)
                for m in range(1, n):
                    vk = np.array([
                        (n-1) * gn[nn, 0] + (m-1) * gn[nn+2, 0],
                        (n-1) * gn[nn, 1] + (m-1) * gn[nn+2, 1]
                    ])
                    
                    ux = (np.cos(nn * PI / 3.0) * uqx[n-1, m-1] - 
                          np.sin(nn * PI / 3.0) * uqy[n-1, m-1])
                    uy = (np.sin(nn * PI / 3.0) * uqx[n-1, m-1] + 
                          np.cos(nn * PI / 3.0) * uqy[n-1, m-1])
                    
                    delta[0] -= 2 * ux * np.sin(vk[0] * rx + vk[1] * ry)
                    delta[1] -= 2 * uy * np.sin(vk[0] * rx + vk[1] * ry)
        
        # Apply displacement
        coord[i, 0] = rx + delta[0] / 2
        coord[i, 1] = ry + delta[1] / 2
        
        # Mirror for layer 2
        coord[i + NDIM // 2, 0] = -coord[i, 0]
        coord[i + NDIM // 2, 1] = coord[i, 1]


def sample_points_bands(numk: int, ncount: int, g1: np.ndarray, g2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Sample k-points for band structure calculation.
    
    Parameters:
    -----------
    numk : int
        Number of k-points
    ncount : int
        Total count of points
    g1, g2 : np.ndarray
        Reciprocal lattice vectors
        
    Returns:
    --------
    points_bz : np.ndarray
        BZ points array
    points_back : np.ndarray
        Back-mapping array
    """
    points_bz = np.zeros((ncount, 8), dtype=int)
    points_back = np.zeros((numk, numk), dtype=int)
    
    # K-points
    vk1 = (g1 + g2) / 3.0
    vm1 = g1 / 2.0
    
    # Transformation matrices
    a = np.array([[g1[0]/numk, g2[0]/numk], [g1[1]/numk, g2[1]/numk]])
    b = np.linalg.inv(a)
    
    icount = 0
    
    # Gamma to K path
    for ivk2 in range(numk//3 + 1):
        ivk1 = ivk2
        points_bz[icount, 0] = ivk1
        points_bz[icount, 1] = ivk2
        points_back[ivk1, ivk2] = icount
        icount += 1
    
    # K to M path
    for ivk in range(numk//6 - 1, 0, -1):
        vk = vm1 + (vk1 - vm1) * ivk / (numk//6)
        ivk_indices = np.round(b @ vk).astype(int)
        points_bz[icount, 0] = ivk_indices[0]
        points_bz[icount, 1] = ivk_indices[1]
        points_back[ivk_indices[0], ivk_indices[1]] = icount
        icount += 1
    
    # M to Gamma path
    for ivk1 in range(numk//2, -1, -1):
        ivk2 = 0
        points_bz[icount, 0] = ivk1
        points_bz[icount, 1] = ivk2
        points_back[ivk1, ivk2] = icount
        icount += 1
    
    return points_bz, points_back


def plot_bands_gamma(bands_tb: np.ndarray, fmu: float, points_back: np.ndarray, 
                    ndim_ev: int, numk: int, ncount: int, vq1: np.ndarray, vq12: np.ndarray,
                    filename: str):
    """
    Plot band structure data.
    
    Parameters:
    -----------
    bands_tb : np.ndarray
        Band data
    fmu : float
        Fermi level
    points_back : np.ndarray
        Back-mapping array
    ndim_ev : int
        Number of eigenvalues
    numk : int
        Number of k-points
    ncount : int
        Total count
    vq1, vq12 : np.ndarray
        Reciprocal lattice vectors
    filename : str
        Output filename
    """
    # K-points
    vk1 = (vq1 + vq12) / 3.0
    vm1 = vq1 / 2.0
    
    # Path lengths
    akg = np.linalg.norm(vk1)
    agm = np.linalg.norm(vm1)
    amk = np.linalg.norm(vk1 - vm1)
    at = akg + agm + amk
    
    # Transformation matrices
    a = np.array([[vq1[0]/numk, vq12[0]/numk], [vq1[1]/numk, vq12[1]/numk]])
    b = np.linalg.inv(a)
    
    with open(filename, 'w') as f:
        # Gamma to K
        for ivk2 in range(numk//3 + 1):
            ivk1 = ivk2
            k_coord = ivk2 * akg / (numk//3) / at
            f.write(f"{k_coord:16.8f}   ")
            for i in range(ndim_ev):
                band_val = bands_tb[i, points_back[ivk1, ivk2]] - fmu
                if i < ndim_ev - 1:
                    f.write(f"{band_val:16.8f}   ")
                else:
                    f.write(f"{band_val:16.8f}\n")
        
        # K to M
        for ivk in range(numk//6, 0, -1):
            vk = vm1 + (vk1 - vm1) * ivk / (numk//6)
            ivk_indices = np.round(b @ vk).astype(int)
            k_coord = (akg + (numk//6 - ivk) * amk / (numk//6)) / at
            f.write(f"{k_coord:16.8f}   ")
            for i in range(ndim_ev):
                band_val = bands_tb[i, points_back[ivk_indices[0], ivk_indices[1]]] - fmu
                if i < ndim_ev - 1:
                    f.write(f"{band_val:16.8f}   ")
                else:
                    f.write(f"{band_val:16.8f}\n")
        
        # M to Gamma
        for ivk1 in range(numk//2, -1, -1):
            ivk2 = 0
            k_coord = 1.0 - ivk1 * agm / (numk//2) / at
            f.write(f"{k_coord:16.8f}   ")
            for i in range(ndim_ev):
                band_val = bands_tb[i, points_back[ivk1, ivk2]] - fmu
                if i < ndim_ev - 1:
                    f.write(f"{band_val:16.8f}   ")
                else:
                    f.write(f"{band_val:16.8f}\n") 