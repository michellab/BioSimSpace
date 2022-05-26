######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2022
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# BioSimSpace is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with BioSimSpace. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

"""
A class for holding restraints.
"""

from .._SireWrappers import Atom
from ..Types import Length, Angle
from .._SireWrappers import System as _System
from ..Units.Length import nanometer
from ..Units.Angle import degree, radian
from ..Units.Energy import kj_per_mol

class Restraint():
    '''The Restraint class which holds the restraint information for the ABFE
    calculations. Currently only Boresch type restraint is supported.
    Boresch restraint is a set of harmonic restraints containing one bond, two
    angle and three dihedrals, which comes from three atoms in the ligand
    (l1, l2, l3) and three atoms in the protein (r1, r2, r3). The restraints
    are arranged in the format of atom1-atom2 (equilibrium value, force constant):
    Bonds: r1-l1 (r0, kr)
    Angles: r2-r1-l1 (thetaA0, kthetaA), r1-l1-l2 (thetaB0, kthetaB)
    Dihedrals: r3-r2-r1-l1 (phiA0, kphiA), r2-r1-l1-l2 (phiB0, kphiB), r1-l1-l2-l3 (phiC0, kphiC)
    '''
    def __init__(self, system, restraint_dict, type='Boresch'):
        if not isinstance(system, _System):
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            # Store a copy of solvated system.
            self._system = system.copy()

        if type.lower() == 'boresch':
            # Test if the atoms are of BioSimSpace._SireWrappers.Atom
            for key in ['r3', 'r2', 'r1', 'l1', 'l2', 'l3']:
                if not isinstance(restraint_dict['anchor_points'][key], Atom):
                    raise ValueError(f"restraint_dict['anchor_points']['{key}'] "
                                     f"must be of type "
                                     f"'BioSimSpace._SireWrappers.Atom'")

            # Test if the equilibrium length of the bond r1-l1 is a length unit
            # Such as angstrom or nanometer
            if not isinstance(restraint_dict['equilibrium_values']['r0'], Length):
                raise ValueError(
                    "restraint_dict['equilibrium_values']['r0'] must be of type 'BioSimSpace.Types.Length'")

            # Test if the equilibrium length of the angle and dihedral is a
            # angle unit such as radian or degree
            for key in ["thetaA0", "thetaB0", "phiA0", "phiB0", "phiC0"]:
                if not isinstance(restraint_dict['equilibrium_values'][key], Angle):
                    raise ValueError(
                        f"restraint_dict['equilibrium_values']['{key}'] must be "
                        f"of type 'BioSimSpace.Types.Angle'")

            # Test if the force constant of the bond r1-l1 is the correct unit
            # Such as kcal/mol/angstrom^2
            dim = restraint_dict['force_constants']['kr'].dimensions()
            if dim != (0, 0, 0, 1, -1, 0, -2):
                raise ValueError(
                    "restraint_dict['force_constants']['kr'] must be of type "
                    "'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Length^2'")

            # Test if the force constant of the angle and dihedral is the correct unit
            # Such as kcal/mol/rad^2
            for key in ["kthetaA", "kthetaB", "kphiA", "kphiB", "kphiC"]:
                dim = restraint_dict['force_constants'][key].dimensions()
                if dim != (-2, 0, 2, 1, -1, 0, -2):
                    raise ValueError(
                        f"restraint_dict['force_constants']['{key}'] must be of type "
                        f"'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Angle^2'")

            self._Boresch = restraint_dict

        else:
            raise NotImplementedError(f'Restraint type {type} not implemented '
                                      f'yet. Only boresch restraint is supported.')

    def toString(self, engine='Gromacs'):
        '''The method for convert the restraint to a format that could be used
        by MD Engines'''
        if engine.lower() == 'gromacs':
            atom_idx = '{:<10}'
            output = ['[ intermolecular_interactions ]',]
            output.append('[ bonds ]')
            output.append('#;ai         aj      type bA         kA         bB         kB')
            bond_string = '  {idx_0:<10} {idx_1:<10} 6 {eq0:<10} {fc0:<10} {eq1:<10} {fc1:<10}'
            # Bonds: r1-l1 (r0, kr)
            output.append(bond_string.format(
                idx_0=self._system.getIndex(
                    self._Boresch['anchor_points']['r1']) + 1,
                idx_1=self._system.getIndex(
                    self._Boresch['anchor_points']['l1']) + 1,
                eq0='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['r0'] / nanometer)),
                fc0='{:.2f}'.format(0),
                eq1='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['r0'] / nanometer)),
                fc1='{:.2f}'.format(
                    (self._Boresch['force_constants']['kr'] / (kj_per_mol / nanometer ** 2)).value()),
            ))
            output.append('[ angles ]')
            output.append('#;ai         aj         ak      type thA        fcA        thB        fcB')
            angle_string = '  {idx_0:<10} {idx_1:<10} {idx_2:<10} 1 {thA:<10} {fcA:<10} {thB:<10} {fcB:<10}'
            # Angles: r2-r1-l1 (thetaA0, kthetaA)
            output.append(angle_string.format(
                idx_0=self._system.getIndex(
                    self._Boresch['anchor_points']['r2']) + 1,
                idx_1=self._system.getIndex(
                    self._Boresch['anchor_points']['r1']) + 1,
                idx_2=self._system.getIndex(
                    self._Boresch['anchor_points']['l1']) + 1,
                thA='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['thetaA0'] / degree)),
                fcA='{:.2f}'.format(0),
                thB='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['thetaA0'] / degree)),
                fcB='{:.2f}'.format(
                    (self._Boresch['force_constants']['kthetaA'] / (kj_per_mol / radian * radian)).value()),
            ))
            # Angles: r1-l1-l2 (thetaB0, kthetaB)
            output.append(angle_string.format(
                idx_0=self._system.getIndex(
                    self._Boresch['anchor_points']['r1']) + 1,
                idx_1=self._system.getIndex(
                    self._Boresch['anchor_points']['l1']) + 1,
                idx_2=self._system.getIndex(
                    self._Boresch['anchor_points']['l2']) + 1,
                thA='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['thetaB0'] / degree)),
                fcA='{:.2f}'.format(0),
                thB='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['thetaB0'] / degree)),
                fcB='{:.2f}'.format(
                    (self._Boresch['force_constants']['kthetaB'] / (kj_per_mol / radian * radian)).value()),
            ))
            output.append('[ dihedrals ]')
            output.append('#;ai         aj         ak         al      type phiA       fcA        phiB       fcB')
            dihedral_string = '  {idx_0:<10} {idx_1:<10} {idx_2:<10} {idx_3:<10} 2 {phiA:<10} {fcA:<10} {phiB:<10} {fcB:<10}'
            # Dihedrals: r3-r2-r1-l1 (phiA0, kphiA)
            output.append(dihedral_string.format(
                idx_0=self._system.getIndex(
                    self._Boresch['anchor_points']['r3']) + 1,
                idx_1=self._system.getIndex(
                    self._Boresch['anchor_points']['r2']) + 1,
                idx_2=self._system.getIndex(
                    self._Boresch['anchor_points']['r1']) + 1,
                idx_3=self._system.getIndex(
                    self._Boresch['anchor_points']['l1']) + 1,
                phiA='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['phiA0'] / degree)),
                fcA='{:.2f}'.format(0),
                phiB='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['phiA0'] / degree)),
                fcB='{:.2f}'.format(
                    (self._Boresch['force_constants']['kphiA'] / (kj_per_mol / radian * radian)).value()),
            ))
            # Dihedrals: r2-r1-l1-l2 (phiB0, kphiB)
            output.append(dihedral_string.format(
                idx_0=self._system.getIndex(
                    self._Boresch['anchor_points']['r2']) + 1,
                idx_1=self._system.getIndex(
                    self._Boresch['anchor_points']['r1']) + 1,
                idx_2=self._system.getIndex(
                    self._Boresch['anchor_points']['l1']) + 1,
                idx_3=self._system.getIndex(
                    self._Boresch['anchor_points']['l2']) + 1,
                phiA='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['phiB0'] / degree)),
                fcA='{:.2f}'.format(0),
                phiB='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['phiB0'] / degree)),
                fcB='{:.2f}'.format(
                    (self._Boresch['force_constants']['kphiB'] / (kj_per_mol / radian * radian)).value()),
            ))
            # Dihedrals: r1-l1-l2-l3 (phiC0, kphiC)
            output.append(dihedral_string.format(
                idx_0=self._system.getIndex(
                    self._Boresch['anchor_points']['r1']) + 1,
                idx_1=self._system.getIndex(
                    self._Boresch['anchor_points']['l1']) + 1,
                idx_2=self._system.getIndex(
                    self._Boresch['anchor_points']['l2']) + 1,
                idx_3=self._system.getIndex(
                    self._Boresch['anchor_points']['l3']) + 1,
                phiA='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['phiC0'] / degree)),
                fcA='{:.2f}'.format(0),
                phiB='{:.3f}'.format(
                    (self._Boresch['equilibrium_values']['phiC0'] / degree)),
                fcB='{:.2f}'.format(
                    (self._Boresch['force_constants']['kphiC'] / (kj_per_mol / radian * radian)).value()),
            ))
            return '\n'.join(output)
        else:
            raise NotImplementedError(f'MD Engine {engine} not implemented '
                                      f'yet. Only Gromacs is supported.')