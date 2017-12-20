"""
@package biosimspace
@author  Lester Hedges
@brief   Common helper functions for the various "*_process.py" modules.
"""

from Sire.Base import *
from Sire.Mol import *

from operator import add, sub

def _compute_box_size(system, tol=0.3, buffer=0.1):
    """ Compute the box size and origin from the atomic coordinates.

        Keyword arguments:

        system -- A Sire molecular system.
        tol    -- The tolerance for determining whether the box is square
                  and whether the origin lies at (0, 0, 0).
        buffer -- The percentage by which to expand the box to account for
                  periodic wrapping.
    """

    # Store the list of molecule indices.
    mol_nums = system.molNums()

    # Initialise the min and max box size for each dimension.
    box_min = [1000000]  * 3
    box_max = [-1000000] * 3

    # Loop over all of the molecules.
    for num in mol_nums:

        # Loop over all atoms in the molecule.
        for atom in system[num].atoms():

            # Extract the atomic coordinates.
            try:
                coord = atom.property("coordinates")

            except UserWarning:
               raise

            # Check coordinates against the current min/max.
            for x in range(0, 3):

               if coord[x] < box_min[x]:
                   box_min[x] = coord[x]

               elif coord[x] > box_max[x]:
                   box_max[x] = coord[x]

    # Calculate the base length of the simulation box.
    box_size = list(map(sub, box_max, box_min))

    # Calculate the centre of the box.
    box_origin = [x * 0.5 for x in list(map(add, box_min, box_max))]

    # Store the base length with the maximum size.
    max_size = max(box_size)

    # Loop over all box dimensions.
    for x in range(0, 3):

        # Assume the box is square if the base lengths are similar.
        if box_size[x] > (1 - tol) * max_size:
            box_size[x] = max_size

        # Assume the origin is at zero if the centre of mass is
        # close  to (0, 0, 0)
        if box_origin[x] / box_size[x] < tol:
            box_origin[x] = 0

    # Add a buffer to the box size to account for atom wrapping.
    box_size = [x * (1 + buffer) for x in box_size]

    return (tuple(box_size), tuple(box_origin))
