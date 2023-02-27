import sympy as sp


def determinant(matrix: sp.Matrix):
    return matrix.det()


def kobayashi(matrix: sp.Matrix):
    return (matrix.transpose() @ matrix.inverse_CH()).trace()


def signature(matrix: sp.Matrix):
    raise NotImplementedError()
