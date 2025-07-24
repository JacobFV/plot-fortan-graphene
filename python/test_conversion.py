#!/usr/bin/env python3
"""
Test script to verify the Python conversion of the Fortran band structure code.
This runs a simplified version to check that all modules import and basic functions work.
"""

import numpy as np
import sys

def test_imports():
    """Test that all modules can be imported."""
    print("Testing imports...")
    try:
        import lapack_routines
        print("✓ lapack_routines imported successfully")
    except ImportError as e:
        print(f"✗ Error importing lapack_routines: {e}")
        return False
    
    try:
        import hamiltonian
        print("✓ hamiltonian imported successfully")
    except ImportError as e:
        print(f"✗ Error importing hamiltonian: {e}")
        return False
    
    try:
        from hamiltonian import (
            NUMK, NTHETA, NDIM, NUMB, NCB, NRELAX,
            wigner_seitz, complex_h0, sample_points_bands
        )
        print("✓ hamiltonian functions imported successfully")
    except ImportError as e:
        print(f"✗ Error importing hamiltonian functions: {e}")
        return False
    
    return True


def test_basic_functions():
    """Test basic function calls."""
    print("\nTesting basic functions...")
    
    from hamiltonian import (
        NTHETA, NDIM, A1, A2, VKD, PI, 
        wigner_seitz, ft, ftperp, sample_points_bands
    )
    
    # Test hopping functions
    try:
        t_val = ft(0.5)
        t_perp_val = ftperp(1.0)
        print(f"✓ Hopping functions work: ft(0.5) = {t_val:.6f}, ftperp(1.0) = {t_perp_val:.6f}")
    except Exception as e:
        print(f"✗ Error in hopping functions: {e}")
        return False
    
    # Test Wigner-Seitz generation (small system)
    try:
        # Simple lattice vectors
        t1 = np.array([1.0, 0.0])
        t2 = np.array([0.5, np.sqrt(3)/2])
        t3 = t2 - t1
        cs = 0.9
        sn = np.sqrt(1 - cs**2)
        
        coord = wigner_seitz(t1, t2, t3, cs, sn)
        print(f"✓ Wigner-Seitz generation works: created {coord.shape[0]} atoms")
    except Exception as e:
        print(f"✗ Error in Wigner-Seitz generation: {e}")
        return False
    
    return True


def test_diagonalization():
    """Test matrix diagonalization."""
    print("\nTesting matrix diagonalization...")
    
    from lapack_routines import diagonalize
    
    try:
        # Create a small test Hermitian matrix
        n = 4
        H = np.random.random((n, n)) + 1j * np.random.random((n, n))
        H = H + H.conj().T  # Make Hermitian
        
        # Test full diagonalization
        eigenvals, eigenvecs = diagonalize(H.copy(), ev='V')
        print(f"✓ Full diagonalization works: {len(eigenvals)} eigenvalues")
        
        # Test partial diagonalization
        H_copy = H.copy()
        eigenvals_partial, _ = diagonalize(H_copy, ev='N', il=1, iu=2)
        eigenvals_partial = np.asarray(eigenvals_partial)  # Ensure it's an array
        if eigenvals_partial.ndim == 0:  # If it's a scalar, make it 1D
            eigenvals_partial = np.array([eigenvals_partial])
        print(f"✓ Partial diagonalization works: {len(eigenvals_partial)} eigenvalues")
        
    except Exception as e:
        print(f"✗ Error in diagonalization: {e}")
        return False
    
    return True


def test_k_point_sampling():
    """Test k-point sampling."""
    print("\nTesting k-point sampling...")
    
    from hamiltonian import sample_points_bands, NUMK
    
    try:
        # Simple reciprocal lattice vectors
        g1 = np.array([1.0, 0.0])
        g2 = np.array([0.5, np.sqrt(3)/2])
        
        ncount = NUMK // 3 + NUMK // 2 + NUMK // 6 + 1
        points_bz, points_back = sample_points_bands(NUMK, ncount, g1, g2)
        
        print(f"✓ K-point sampling works: {points_bz.shape[0]} k-points generated")
        
    except Exception as e:
        print(f"✗ Error in k-point sampling: {e}")
        return False
    
    return True


def main():
    """Run all tests."""
    print("=" * 50)
    print("Testing Python conversion of Fortran band structure code")
    print("=" * 50)
    
    tests = [
        test_imports,
        test_basic_functions,
        test_diagonalization,
        test_k_point_sampling
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                print("Test failed!")
        except Exception as e:
            print(f"Test crashed with error: {e}")
    
    print("\n" + "=" * 50)
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("✓ All tests passed! The conversion appears to be working correctly.")
        print("\nYou can now run the main calculation with:")
        print("python plot_bands.py")
    else:
        print("✗ Some tests failed. Please check the error messages above.")
        print("Make sure all dependencies are installed:")
        print("pip install -r requirements.txt")
    
    print("=" * 50)
    
    return passed == total


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 