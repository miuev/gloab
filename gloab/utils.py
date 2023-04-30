import numpy as np
import ase
import ase.io
import ase.neighborlist
from ase.data import covalent_radii

# create a default radii array to use in get_bonds
default_radii = covalent_radii.copy()

# defining custom radii to produce expected metal oxide connectivity
# Oxygen
default_radii[8] = 0.89
# Titanium
default_radii[22] = 1.38
# Vanadium
default_radii[23] = 1.44
# Molybdenum
default_radii[42] = 1.7

def get_ase_atoms_obj(atoms):
    """Creates ase.Atoms object from geometry file
    - will return <atoms> if its already an ase.Atoms object
    Args:
    atoms (str): path to geometry file
    Returns:
    ase.Atoms: atoms object
    """
    # if atoms is a str, assume path to geometry file
    if isinstance(atoms, str):
        # read in geometry
        atoms = ase.io.read(atoms)

    # ensure atoms is an ase.Atoms object
    if not isinstance(atoms, ase.Atoms):
        raise TypeError('atoms must be an ase.Atoms object')

    return atoms


def get_bonds(atoms, radii = None, scale = 1):
    """
    Finds bonds between atoms based on bonding radii
    Args:
    atoms (ase.Atoms): atoms object
    KArgs:
    radii (np.ndarray): bonding radii of atoms
                        if bonding radii overlap, a bond is drawn
                        (Default: gloab.utils.default_radii)
    scale (float): scales bonding radii array
                   (Default: 1.0)
    """
    # remove periodic boundaries
    atoms = atoms.copy()

    # create atom-specific radii array
    if radii == None:
        spec_radii = default_radii[atoms.numbers]
    elif ((type(radii) == float) or (type(radii) == int)):
        spec_radii = np.ones(len(atoms.numbers)) * radii

    # create neighborlist object
    n = ase.neighborlist.NeighborList(spec_radii * scale, skin=0,
                                      self_interaction=False)
    n.update(atoms)
    if not n.nneighbors:
        return []

    bonds = np.zeros((n.nneighbors, 2), int)
    spot1 = 0
    for atomi in range(len(atoms)):
        # get neighbors of atomi
        neighs = n.get_neighbors(atomi)[0]

        # find second cutoff in matrix
        spot2 = spot1 + len(neighs)

        # add bonds to matrix
        bonds[spot1:spot2, 0] = atomi
        bonds[spot1:spot2, 1] = neighs

        # shift down matrix
        spot1 = spot2

        # once all bonds have been found break loop
        if spot1 == n.nneighbors:
            break

    return bonds
