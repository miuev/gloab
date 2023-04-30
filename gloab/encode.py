from ase.geometry import get_distances
from ase.data import atomic_masses
from ase.io.trajectory import TrajectoryWriter

import warnings

import numpy as np

import gloab.utils as utils

class Laplacian:
    """ reads an ASE atoms object and constructs its related network theoretical objects """
    
    def __init__(self, atoms, basis: str = 'adjacency', cutoff: float = None, normalize: bool = False, radii = None):
        """ build the network representation of atoms and bonds of interest

            ~args~
            atoms (ase.Atoms | str): ase.Atoms obj of molecule or path to geometry file
            basis (str): the type of laplacian that should be considered
                         'adjacency' - standard adjacency connectivity matrix
                         'distances' - inverse distance matrix w/ arb. cutoff
            
            cutoff (float): length in Angstroms below which pairwise distances are retained,
                            if 'None', then all pairwise distances are considered
            normalize (bool): flag to consider symmetric normalized laplacian
                              if False, no normalization
                              if True, normalization
            radii (list of floats): radii to consider for neighbor assignments (overlapping)
                                    if 'None', then default radii are used (utils.py)
                                    if float, then all radii are set to this value
        """        

        self.atoms = utils.get_ase_atoms_obj(atoms)
        self.bonds = utils.get_bonds(self.atoms, radii = radii)
        self.basis = basis
        self.cutoff = cutoff
        self.normalize = normalize
        self._adjacency = None
        self._degree = None
        self._laplacian = None
        self._distances = None
        self._eig = None
        self._vec = None

    @property
    def adjacency(self):
        
        self._adjacency = self.get_adjacency()

        return self._adjacency
    
    @property
    def degree(self):

        self._degree = self.get_degree()

        return self._degree

    @property
    def laplacian(self):
        
        self._laplacian = self.get_laplacian()
        
        return self._laplacian

    @property
    def distances(self):

        self._distances = self.get_distances()

        return self._distances

    @property
    def eig(self):
        
        self.get_spectrum()

        return self._eig

    @property
    def vec(self):
        
        self.get_spectrum()

        return self._vec

    @property
    def gap(self):
        spectral_gap = np.min(self.eig[self.eig != 0])
        return spectral_gap
        
    @property
    def algebraic_connectivity(self):
        return self.eig[1]

    @property
    def lel(self):
        energy = np.sqrt(self.eig).sum()
        return energy

    @property
    def global_connectivity(self):
        gc = (np.sqrt(self.algebraic_connectivity)+np.sqrt(self.gap))/self.lel
        return gc

    @property
    def us_connectivity(self):
        us = (np.sqrt(self.algebraic_connectivity)+np.sqrt(self.gap))
        return us

    def get_adjacency(self):
        """ creates adjacency matrix that represents the system given by self.atoms

            ~returns~
            np.ndarray: natom x natom matrix with 0's for no bonds, 1 for single bond,
                        and 2 for a periodic and a non-periodic bond, matrix is in
                        adjacency format
                        N.B. - this matrix is based off of the radii supplied in 
                        the initial object generation (e.g. default in gloab.utils or 
                        single float supplied when calling gloab.encode.Laplacian)

        """

        # create natom x natom matrix of 0's
        adj_mat = np.zeros((len(self.atoms), len(self.atoms))).astype(int)

        # check for multi-edge connections through periodic images
        # sort allows matching of mirrored pairs
        edge_count = np.unique(np.sort(self.bonds, axis=1), axis=0, return_counts=True)

        # set all bonded positions to edge connectivity
        adj_mat[edge_count[0][:, 0], edge_count[0][:, 1]] = edge_count[1]

        # must do both ways to make matrix symmetric
        adj_mat[edge_count[0][:, 1], edge_count[0][:, 0]] = edge_count[1]
        
        return adj_mat

    def get_degree(self):
        """ creates degree vector that represents the system given by self.atoms

        ~returns~
        np.ndarrray: natom array with entries corresponding to atoms numbers in self.atoms
                     each entry is the number of neighbors the specific atom has, according
                     to the overlapping radii (either gloab.utils.default_radii or supplied float)
        """
    
        adj_mat = self.adjacency

        # create degree vector for diagonal
        deg_vec = []
        for i in np.arange(len(self.atoms)):
            deg_vec.append(adj_mat[:,i].sum())
        deg_vec = np.array(deg_vec)

        return deg_vec

    def get_laplacian(self):
        """ creates laplacian matrix that represents the system given by self.atoms

            ~returns~
            np.ndarray: natom x natom matrix with off-diagonal elements being 0's for no bonds,
                        -1 for single bond, and -2 for a periodic and a non-periodic bond,
                        and on-diagonal elements being the positive sums of the corresponding
                        rows / columns, matrix is in laplacian format
                        N.B. - this matrix is based off of the radii supplied in 
                        the initial object generation (e.g. default in gloab.utils or 
                        single float supplied when calling gloab.encode.Laplacian)        
        """
        
        adj_mat = self.adjacency

        # create degree vector for diagonal
        deg_vec = []
        for i in np.arange(len(self.atoms)):
            deg_vec.append(adj_mat[:,i].sum())
        deg_vec = np.array(deg_vec)

        # create laplacian matrix
        lap_mat = np.diag(deg_vec) - adj_mat
        
        if self.normalize == True:
            deg_mat_root = np.nan_to_num(np.diag(np.power(deg_vec, -0.5)), posinf = 0)
            lap_mat = deg_mat_root @ lap_mat @ deg_mat_root

        return lap_mat

    def get_distances(self):
        """ creates matrix that contains distance between each atom in system

            ~returns~
            np.ndarray: natom x natom matrix with off-diagonal elements being distance
                        and on-diagonal elements being 0, matrix is in adjacency format
        """

        distances = get_distances(self.atoms.positions)[1]
        
        np.multiply

        return distances

    def get_spectrum(self):
        """ calculates the eigenvalues and eigenfunctions associates with a particular
            Laplacian that corresponds with the atoms and bonds of interest
            
            ~returns (directly from numpy)~
            w : (..., M) ndarray
                The eigenvalues in ascending order, each repeated according to
                its multiplicity.
            v : {(..., M, M) ndarray, (..., M, M) matrix}
                The column ``v[:, i]`` is the normalized eigenvector corresponding
                to the eigenvalue ``w[i]``.  Will return a matrix object if `a` is
                a matrix object.
        
        """

        if self.basis == 'adjacency':

            mat = self.laplacian

        elif self.basis == 'distances':

            dis_mat = self.distances
            
            # create degree vector for diagonal
            deg_vec = []
            for i in np.arange(len(self.atoms)):
                deg_vec.append(dis_mat[:,i].sum())
            deg_vec = np.array(deg_vec)
            
            # create laplacian matrix
            mat = np.diag(deg_vec) - dis_mat
            
            if self.normalize == True:
                deg_mat_root = np.diag(np.power(deg_vec, -0.5))
                mat = deg_mat_root @ mat @ deg_mat_root

        else:
            warnings.warn('Requested basis is not known. Please check you have entered an appropriate basis.')
            return

        eig, vec = np.linalg.eigh(mat)

        eig[eig < 1E-12] = 0
        self._eig = eig
        self._vec = vec

        return eig, vec