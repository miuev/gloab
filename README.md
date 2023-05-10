## gloab - **g**raph **l**aplacians **o**f **a**toms and **b**onds

This repository contains code which can be used to compute global structural descriptors for arbitrary atomic systems.
Atomic position files are interpreted as matrices cataloging all atom-to-atom connections.
Using spectral graph theory, these matrices are then converted to scalars which quantify global connectivity.
The `atomic simulation environment` (`ASE`) is used to read position files, evaluate atomic distances, assign neighbors, and construct matrices. 
As a result, any file format that is understood by `ASE` can be used by `gloab`.
Algebraic manipulations are performed using `numpy`.

A minimal example of code use is as follows:

```python
from ase.io import read
from gloab.encode import Laplacian

atoms = read('/some/path/to/position/file')

system = Laplacian(atoms = atoms)

gc = system.global_connectivity
```

In this example, a global connectivity value that corresponds to `atoms` is assigned to the `gc` variable.
`atoms` is a prototypical `ASE` atoms object.
`system` is an instance that can be used to generate various connectivity properties which are representative of `atoms`.
As shown through the assignment of `gc`, connectivity properties are called as properties of the `system` instance.

`Laplacian` can be called with additional flags that change underlying behavior.
The first relevant flag is `basis`, which can be used to select either neighbor-based matrices or distance-based matrices as the basis for calculations.
These are respectively selected by using ```basis = adjacency``` or ```basis = distances```.

The second flag is `normalize`, which takes a Boolean argument.
If `True`, the generated Laplacian matrix is the symmetric normalized Laplacian.
If `False`, it is instead an unnormalized Laplacian.

The third flag is `radii`, which is used to control atomic radii assignments.
When using `adjacency` as a basis, `gloab` labels atom neighbors based on overlapping atomic radii.
These can be a default set of atomic radii (as shipped with `ASE`), or can be customized.
A slightly customized set of radii are used as the default in this repository, as described in a recent publication from our group **cite**.
Alternatively, a single universal radius can be used for all elements by supplying a single float to the `radii` flag when instantiating a new `Laplacian`.
If `radii` is `None`, then the default radii will be used.

Once instantiated, the Laplacian can be probed for the different properties summarized here:

| Property | Name |
| - | - |
| adjacency matrix | `adjacency` |
| degree matrix | `degree` |
| laplacian matrix | `laplacian` |
| distance matrix | `distances` |
| eigenvalue spectrum | `eig` | 
| eigenfunctions | `vec` |
| spectral gap | `gap` |
| algebraic connectivity | `algebraic_connectivity` |
| Laplacian energy like invariant | `lel` |
| global connectivity | `global_connectivity`|
| unscaled connectivity | `us_connectivity` |
