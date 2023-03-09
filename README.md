# Mock Seifert Matrix

## Background

The Gordon-Litherland linking form [1] of a virtual knot (knot in thickened surface) is the map
```math
\mathscr{L}_F: H_1(F) \times H_1(F) \longrightarrow \mathbb{Z}
```
given by
```math
\mathscr{L}_F(\alpha, \beta) = \ell k(\tau \alpha, \beta).
```

The mock Seifert matrix $A$ of a virtual knot is the Gram matrix of the Gordon-Litherland linking form with repsect to some basis for $H_1(F)$. Hence, representing the same Gordon-Litherland linking form defines equivalence classes of mock Seifert matrices under unimodular congurence, i.e. the existence of unimodular $P$ such that
```math
A' = P^\top A P.
```

## Features

This package computes the mock Seifert matrix of an alternating virtual knot, andcan compute some invariants derived from it, namely:

- The determinant of the mock Seifert matrix, $\operatorname{det} S$.
- The dimension of the mock Seifert matrix, $\operatorname{dim} S$
- The Kobayashi invariant of the mock Seifert matrix, $\operatorname{tr}(S^\top S^{-1})$. See [2].
- The mock Alexander polynomial of the knot (yet to be implemented).

## Usage
Compute mock Seifert matrix:

`python gl.py <gauss code> <flags>`

The gauss code must be alternating.

Flags:
- `-i` prints invariants.
- `-v` verbose.
- `-vv` very verbose.
- `-s` compute instead the symmetrisation of the mock Seifert matrix, correposnding to the Gordon-Litherland pairing $\mathscr{G}_F$ (yet to be implemented).

Compute mock Seifert matrix and invariants for many knots:

`python catalogue.py <flags>`

Flags:
- `-s` compute instead the symmetrisation of the mock Seifert matrix, correposnding to the Gordon-Litherland pairing $\mathscr{G}_F$ (yet to be implemented).

## References


[1] https://arxiv.org/abs/2301.05946

[2] https://arxiv.org/abs/1904.04397
