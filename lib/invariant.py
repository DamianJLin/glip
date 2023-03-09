import sympy as sp


def dimension(matrix: sp.Matrix):
    dims = sp.shape(matrix)
    assert len(dims) == 2 and dims[0] == dims[1]
    return dims[0]


def determinant(matrix: sp.Matrix):
    return matrix.det()


def kobayashi(matrix: sp.Matrix):
    return (matrix.transpose() @ matrix.inverse_CH()).trace()


def mock_alexander(matrix: sp.Matrix, for_sort=False):
    t = sp.symbols('t')
    poly = sp.Poly(
        (t * matrix - matrix.transpose()).det()
    )

    if not for_sort:
        return poly

    # Return polynomial as tuple of (degree, coefficient) pairs.
    elif for_sort:
        return tuple(
            (poly.total_degree() - i, c) for i, c in enumerate(poly.all_coeffs())
        )
