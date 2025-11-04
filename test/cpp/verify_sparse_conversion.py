#!/usr/bin/env python3
"""
Verify expected CSR/CSC conversion results using scipy.
Used to generate test expectations for BIOMTable::ToCSR() and ToCSC().
"""

import numpy as np
from scipy.sparse import coo_matrix


def print_csr_csc(name, row, col, data, shape):
    """Print CSR and CSC representations of a COO matrix."""
    print(f"\n{'='*60}")
    print(f"Test Case: {name}")
    print(f"{'='*60}")
    print(f"Shape: {shape[0]} features (rows) × {shape[1]} samples (columns)")
    print(f"\nCOO Input (sorted by row, col):")
    print(f"  row indices:    {row}")
    print(f"  col indices:    {col}")
    print(f"  values:         {data}")

    # Create COO matrix
    coo = coo_matrix((data, (row, col)), shape=shape)

    # Print dense representation for visualization
    print(f"\nDense Matrix:")
    dense = coo.toarray()
    for i, row_data in enumerate(dense):
        print(f"  F{i}: {row_data}")

    # Convert to CSR (observation/matrix in BIOM)
    csr = coo.tocsr()
    print(f"\nCSR (observation/matrix):")
    print(f"  data:    {list(csr.data)}")
    print(f"  indices: {list(csr.indices)}")
    print(f"  indptr:  {list(csr.indptr)}")

    # Convert to CSC (sample/matrix in BIOM)
    csc = coo.tocsc()
    print(f"\nCSC (sample/matrix):")
    print(f"  data:    {list(csc.data)}")
    print(f"  indices: {list(csc.indices)}")
    print(f"  indptr:  {list(csc.indptr)}")


def main():
    print("Generating test expectations for BIOMTable sparse conversions")

    # Test 1: Simple 2x3 matrix
    # Matrix:
    #         S0   S1   S2
    #   F0    1.0  2.0  4.0
    #   F1    3.0  0.0  0.0
    row = np.array([0, 0, 0, 1], dtype=np.int32)
    col = np.array([0, 1, 2, 0], dtype=np.int32)
    data = np.array([1.0, 2.0, 4.0, 3.0])
    print_csr_csc("Simple 2x3 matrix", row, col, data, shape=(2, 3))

    # Test 2: Empty matrix
    row = np.array([], dtype=np.int32)
    col = np.array([], dtype=np.int32)
    data = np.array([])
    print_csr_csc("Empty matrix", row, col, data, shape=(0, 0))

    # Test 3: Single row (1 feature, 3 samples)
    row = np.array([0, 0, 0], dtype=np.int32)
    col = np.array([0, 1, 2], dtype=np.int32)
    data = np.array([1.0, 2.0, 3.0])
    print_csr_csc("Single row", row, col, data, shape=(1, 3))

    # Test 4: Single column (3 features, 1 sample)
    row = np.array([0, 1, 2], dtype=np.int32)
    col = np.array([0, 0, 0], dtype=np.int32)
    data = np.array([1.0, 2.0, 3.0])
    print_csr_csc("Single column", row, col, data, shape=(3, 1))

    # Test 5: 3x3 with sparse data
    # Matrix:
    #         S0   S1   S2
    #   F0    1.0  0.0  0.0
    #   F1    0.0  2.0  0.0
    #   F2    0.0  0.0  3.0
    row = np.array([0, 1, 2], dtype=np.int32)
    col = np.array([0, 1, 2], dtype=np.int32)
    data = np.array([1.0, 2.0, 3.0])
    print_csr_csc("3x3 diagonal", row, col, data, shape=(3, 3))

    # Test 6: Single value
    row = np.array([0], dtype=np.int32)
    col = np.array([0], dtype=np.int32)
    data = np.array([5.0])
    print_csr_csc("Single value", row, col, data, shape=(1, 1))


if __name__ == "__main__":
    main()
