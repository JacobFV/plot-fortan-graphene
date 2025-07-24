import numpy as np
from scipy.linalg import eigh, eig
from typing import Tuple, Optional, Union

def diagonalize(A: np.ndarray, 
                eigenvalues: Optional[np.ndarray] = None,
                ev: str = 'V', 
                il: Optional[int] = None, 
                iu: Optional[int] = None,
                d: Optional[str] = None) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """
    Diagonalize a Hermitian matrix using LAPACK routines.
    
    Parameters:
    -----------
    A : np.ndarray
        Input Hermitian matrix (complex or real)
    eigenvalues : np.ndarray, optional
        Output array for eigenvalues
    ev : str
        'V' to compute eigenvectors, 'N' for eigenvalues only
    il, iu : int, optional
        Index range for partial diagonalization (1-based, inclusive)
    d : str, optional
        Use divide-and-conquer algorithm if specified
        
    Returns:
    --------
    eigenvals : np.ndarray
        Eigenvalues
    eigenvecs : np.ndarray or None
        Eigenvectors (if ev='V') or None
    """
    
    if ev not in ['V', 'N']:
        raise ValueError("ev must be either 'V' or 'N'")
    
    # Convert 1-based indices to 0-based
    if il is not None and iu is not None:
        il_0 = il - 1
        iu_0 = iu - 1
        subset_by_index = [il_0, iu_0]
    else:
        subset_by_index = None
    
    # Determine if we need eigenvectors
    eigvals_only = (ev == 'N')
    
    # Handle different data types
    if np.iscomplexobj(A):
        # Complex matrix
        if subset_by_index is not None:
            # Partial diagonalization
            eigenvals, eigenvecs = eigh(A, eigvals_only=eigvals_only, 
                                      subset_by_index=subset_by_index)
        else:
            # Full diagonalization
            eigenvals, eigenvecs = eigh(A, eigvals_only=eigvals_only)
    else:
        # Real matrix
        if subset_by_index is not None:
            eigenvals, eigenvecs = eigh(A, eigvals_only=eigvals_only,
                                      subset_by_index=subset_by_index)
        else:
            eigenvals, eigenvecs = eigh(A, eigvals_only=eigvals_only)
    
    # Handle eigenvalues-only case
    if eigvals_only:
        eigenvecs = None
    
    # Update the input matrix with eigenvectors if requested
    if not eigvals_only and subset_by_index is not None:
        # For partial diagonalization, update only the relevant columns
        n_eigs = iu - il + 1
        A[:, :n_eigs] = eigenvecs
    elif not eigvals_only:
        # For full diagonalization, replace entire matrix
        A[:] = eigenvecs
    
    return eigenvals, eigenvecs


def zdiag_red(A: np.ndarray, eigenvalues: np.ndarray, ev: str, il: int, iu: int):
    """
    Diagonalize complex double precision matrix (partial).
    
    Parameters:
    -----------
    A : np.ndarray (complex128)
        Input/output matrix
    eigenvalues : np.ndarray (float64)
        Output eigenvalues
    ev : str
        'V' or 'N'
    il, iu : int
        Index range (1-based)
    """
    eigenvals, eigenvecs = diagonalize(A, eigenvalues, ev, il, iu)
    eigenvalues[:len(eigenvals)] = eigenvals
    return eigenvals


def cdiag_red(A: np.ndarray, eigenvalues: np.ndarray, ev: str, il: int, iu: int):
    """
    Diagonalize complex single precision matrix (partial).
    
    Parameters:
    -----------
    A : np.ndarray (complex64)
        Input/output matrix
    eigenvalues : np.ndarray (float32)
        Output eigenvalues
    ev : str
        'V' or 'N'
    il, iu : int
        Index range (1-based)
    """
    eigenvals, eigenvecs = diagonalize(A, eigenvalues, ev, il, iu)
    eigenvalues[:len(eigenvals)] = eigenvals
    return eigenvals


def zdiag_full(A: np.ndarray, eigenvalues: np.ndarray, ev: str):
    """
    Full diagonalization of complex double precision matrix.
    
    Parameters:
    -----------
    A : np.ndarray (complex128)
        Input/output matrix
    eigenvalues : np.ndarray (float64)
        Output eigenvalues
    ev : str
        'V' or 'N'
    """
    eigenvals, eigenvecs = diagonalize(A, eigenvalues, ev)
    eigenvalues[:] = eigenvals
    return eigenvals


def cdiag_full(A: np.ndarray, eigenvalues: np.ndarray, ev: str):
    """
    Full diagonalization of complex single precision matrix.
    
    Parameters:
    -----------
    A : np.ndarray (complex64)
        Input/output matrix
    eigenvalues : np.ndarray (float32)
        Output eigenvalues
    ev : str
        'V' or 'N'
    """
    eigenvals, eigenvecs = diagonalize(A, eigenvalues, ev)
    eigenvalues[:] = eigenvals
    return eigenvals


def zdiag(A: np.ndarray, eigenvalues: np.ndarray, ev: str, d: Optional[str] = None):
    """
    Diagonalize complex double precision matrix.
    
    Parameters:
    -----------
    A : np.ndarray (complex128)
        Input/output matrix
    eigenvalues : np.ndarray (float64)
        Output eigenvalues
    ev : str
        'V' or 'N'
    d : str, optional
        Use divide-and-conquer if specified
    """
    eigenvals, eigenvecs = diagonalize(A, eigenvalues, ev, d=d)
    eigenvalues[:] = eigenvals
    return eigenvals


def cdiag(A: np.ndarray, eigenvalues: np.ndarray, ev: str, d: Optional[str] = None):
    """
    Diagonalize complex single precision matrix.
    
    Parameters:
    -----------
    A : np.ndarray (complex64)
        Input/output matrix
    eigenvalues : np.ndarray (float32)
        Output eigenvalues
    ev : str
        'V' or 'N'
    d : str, optional
        Use divide-and-conquer if specified
    """
    eigenvals, eigenvecs = diagonalize(A, eigenvalues, ev, d=d)
    eigenvalues[:] = eigenvals
    return eigenvals 