import logging
import numpy as np
from stips.utilities import overlapadd2
from scipy.signal import convolve2d


def test_convolution():
    """
    This function tests the overlapadd2
    convolve function.
    """
    # Default tests:
    # Construct inputs and correct output
    array = np.random.randn(33, 55)
    kernel = np.random.randn(4, 5)
    correct_output = convolve2d(array, kernel)

    # Test different L values
    assert(np.allclose(correct_output, overlapadd2(array, kernel, L=[12, 12])))
    assert(np.allclose(correct_output, overlapadd2(array, kernel, L=[12, 120])))
    assert(np.allclose(correct_output, overlapadd2(array, kernel, L=[90, 120])))

    # Swap kernel and array
    assert(np.allclose(correct_output, overlapadd2(kernel, array, L=[190, 220])))
    assert(np.allclose(correct_output, overlapadd2(kernel, array, L=[1, 1])))

    # Convolve without breaking down
    assert(np.allclose(correct_output, overlapadd2(kernel, array)))
    assert(np.allclose(correct_output, overlapadd2(array, kernel)))

    # Transpose and convolve
    assert(np.allclose(convolve2d(array.T, kernel.T),
                       overlapadd2(array.T, kernel.T, L=[190, 220])))

    # Second set of tests:
    # Construct inputs and correct output
    array = np.random.randn(33, 55) + 1j * np.random.randn(33, 55)
    kernel = np.random.randn(4, 5) + 1j * np.random.randn(4, 5)
    correct_output = convolve2d(array, kernel)
    assert(np.allclose(correct_output, overlapadd2(kernel, array)))
    assert(np.allclose(correct_output, overlapadd2(array, kernel)))

    # Test storage for the output (y):
    M = np.array(kernel.shape)
    Na = np.array(array.shape)
    y = np.zeros(M + Na - 1, dtype=array.dtype)
    overlapadd2(kernel, array, y=y)
    assert(np.allclose(correct_output, y))
    assert(np.allclose(overlapadd2(array, kernel), y))

    # Test Nfft:
    L = M * 100
    Nfft = 2 ** np.ceil(np.log2(L + M - 1)).astype(int)
    assert (np.allclose(correct_output, overlapadd2(array, kernel, Nfft=Nfft)))

    # Test logging:
    overlapadd2(array, kernel, verbose=True, logger=logging)
