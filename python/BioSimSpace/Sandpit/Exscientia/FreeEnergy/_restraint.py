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

import numpy as np
from sire.legacy.Units import k_boltz as _k_boltz
from sire.legacy.Units import meter3 as _meter3
from sire.legacy.Units import nanometer3 as _nanometer3
from sire.legacy.Units import mole as _mole

from .._SireWrappers import Atom as _Atom
from ..Types import Length as _Length
from ..Types import Angle as _Angle
from ..Types import Temperature as _Temperature
from .._SireWrappers import System as _System
from ..Units.Length import nanometer as _nanometer
from ..Units.Area import nanometer2 as _nanometer2
from ..Units.Temperature import kelvin as _kelvin
from ..Units.Angle import degree as _degree
from ..Units.Angle import radian as _radian
from ..Units.Energy import kj_per_mol as _kj_per_mol
from ..Units.Energy import kcal_per_mol as _kcal_per_mol

class Restraint():
    '''The Restraint class which holds the restraint information for the ABFE
    calculations. Currently only Boresch type restraint is supported.

    Boresch restraint is a set of harmonic restraints containing one bond, two
    angle and three dihedrals, which comes from three atoms in the ligand
    (l1, l2, l3) and three atoms in the protein (r1, r2, r3). The restraints
    are represented in the format of:
    atom1-atom2-... (equilibrium value, force constant)

    The nomenclature:
    Bonds: r1-l1 (r0, kr)
    Angles: r2-r1-l1 (thetaA0, kthetaA), r1-l1-l2 (thetaB0, kthetaB)
    Dihedrals: r3-r2-r1-l1 (phiA0, kphiA), r2-r1-l1-l2 (phiB0, kphiB), r1-l1-l2-l3 (phiC0, kphiC)

    The restraint_dict has the following compact format.

    restraint_dict = {
        "anchor_points":{"r1": BioSimSpace._SireWrappers.Atom,
                         "r2": BioSimSpace._SireWrappers.Atom,
                         "r3": BioSimSpace._SireWrappers.Atom,
                         "l1": BioSimSpace._SireWrappers.Atom,
                         "l2": BioSimSpace._SireWrappers.Atom,
                         "l3": BioSimSpace._SireWrappers.Atom},
        "equilibrium_values":{"r0": BioSimSpace.Types.Length,
                              "thetaA0": BioSimSpace.Types.Angle,
                              "thetaB0": BioSimSpace.Types.Angle,
                              "phiA0": BioSimSpace.Types.Angle,
                              "phiB0": BioSimSpace.Types.Angle,
                              "phiC0": BioSimSpace.Types.Angle},
        "force_constants":{"kr": BioSimSpace.Types.Energy / BioSimSpace.Types.Area,
                           "kthetaA": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kthetaB": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kphiA": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kphiB": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area),
                           "kphiC": BioSimSpace.Types.Energy / (BioSimSpace.Types.Area * BioSimSpace.Types.Area)}}

    '''
    def __init__(self, system, restraint_dict, temperature, restraint_type='Boresch'):
        """Constructor.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.

           restraint_dict : dict
               The dict for holding the restraint.

           temperature : :class:`System <BioSimSpace.Types.Temperature>`
               The temperature of the system

           restraint_type : str
               The type of the restraint. (`Boresch`, )
        """
        if not isinstance(temperature, _Temperature):
            raise ValueError(
                "temperature must be of type 'BioSimSpace.Types.Temperature'")
        else:
            self.T = temperature

        if restraint_type.lower() == 'boresch':
            self._restraint_type = 'boresch'
            # Test if the atoms are of BioSimSpace._SireWrappers.Atom
            for key in ['r3', 'r2', 'r1', 'l1', 'l2', 'l3']:
                if not isinstance(restraint_dict['anchor_points'][key], _Atom):
                    raise ValueError(f"restraint_dict['anchor_points']['{key}'] "
                                     f"must be of type "
                                     f"'BioSimSpace._SireWrappers.Atom'")

            # Test if the equilibrium length of the bond r1-l1 is a length unit
            # Such as angstrom or nanometer
            if not isinstance(restraint_dict['equilibrium_values']['r0'], _Length):
                raise ValueError(
                    "restraint_dict['equilibrium_values']['r0'] must be of type 'BioSimSpace.Types.Length'")

            # Test if the equilibrium length of the angle and dihedral is a
            # angle unit such as radian or degree
            for key in ["thetaA0", "thetaB0", "phiA0", "phiB0", "phiC0"]:
                if not isinstance(restraint_dict['equilibrium_values'][key], _Angle):
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
        else:
            raise NotImplementedError(f'Restraint type {type} not implemented '
                                      f'yet. Only boresch restraint is supported.')

        self._restraint_dict = restraint_dict
        self.system = system

    @property
    def system(self):
        return self._system

    @system.setter
    def system(self, system):
        """Update the system object.

           Parameters
           ----------

           system : :class:`System <BioSimSpace._SireWrappers.System>`
               The molecular system.
        """
        if not isinstance(system, _System):
            raise TypeError("'system' must be of type 'BioSimSpace._SireWrappers.System'")
        else:
            if self._restraint_type == 'boresch':
                # Check if the ligand atoms are decoupled.
                # Find the decoupled molecule, assume that only one can be
                # decoupled.
                (decoupled_mol,) = system.getDecoupledMolecules()
                for key in ['l1', 'l2', 'l3']:
                    atom = self._restraint_dict['anchor_points'][key]
                    # Discussed in https://github.com/michellab/BioSimSpace/pull/337
                    if atom._sire_object.molecule().number() != decoupled_mol._sire_object.number():
                        raise ValueError(
                            f'The ligand atom {key} is not from decoupled moleucle.')
                for key in ['r1', 'r2', 'r3']:
                    atom = self._restraint_dict['anchor_points'][key]
                    if not atom in system:
                        raise ValueError(
                            f'The protein atom {key} is not in the system.')

            # Store a copy of solvated system.
            self._system = system.copy()

    def _gromacs_boresch(self):
        '''Format the Gromacs string for boresch restraint.'''
        # Format the atoms into index list
        def format_index(key_list):
            formated_index = []
            for key in key_list:
                formated_index.append('{:<10}'.format(
                    self._system.getIndex(
                        self._restraint_dict['anchor_points'][key]) + 1))
            return ' '.join(formated_index)

        parameters_string = '{eq0:<10} {fc0:<10} {eq1:<10} {fc1:<10}'

        # Format the parameters for the bonds
        def format_bond(equilibrium_values, force_constants):
            converted_equ_val = \
                self._restraint_dict['equilibrium_values'][
                    equilibrium_values] / _nanometer
            converted_fc = \
                self._restraint_dict['force_constants'][force_constants] / (
                            _kj_per_mol / _nanometer ** 2)
            return parameters_string.format(
                eq0='{:.3f}'.format(converted_equ_val),
                fc0='{:.2f}'.format(0),
                eq1='{:.3f}'.format(converted_equ_val),
                fc1='{:.2f}'.format(converted_fc),
            )

        # Format the parameters for the angles and dihedrals
        def format_angle(equilibrium_values, force_constants):
            converted_equ_val = \
                self._restraint_dict['equilibrium_values'][
                    equilibrium_values] / _degree
            converted_fc = \
                self._restraint_dict['force_constants'][force_constants] / (
                            _kj_per_mol / (_radian * _radian))
            return parameters_string.format(
                eq0='{:.3f}'.format(converted_equ_val),
                fc0='{:.2f}'.format(0),
                eq1='{:.3f}'.format(converted_equ_val),
                fc1='{:.2f}'.format(converted_fc),
            )

        # basic format of the Gromacs string
        master_string = '  {index} {func_type} {parameters}'

        def write_bond(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=6,
                parameters=format_bond(equilibrium_values,
                                       force_constants),
            )

        def write_angle(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=1,
                parameters=format_angle(equilibrium_values,
                                        force_constants),
            )

        def write_dihedral(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=2,
                parameters=format_angle(equilibrium_values,
                                        force_constants),
            )

        # Writing the string
        output = ['[ intermolecular_interactions ]', ]

        output.append('[ bonds ]')
        output.append(
            '; ai         aj      type bA         kA         bB         kB')
        # Bonds: r1-l1 (r0, kr)
        output.append(
            write_bond(('r1', 'l1'), 'r0', 'kr'))

        output.append('[ angles ]')
        output.append(
            '; ai         aj         ak      type thA        fcA        thB        fcB')
        # Angles: r2-r1-l1 (thetaA0, kthetaA)
        output.append(
            write_angle(('r2', 'r1', 'l1'), 'thetaA0', 'kthetaA'))
        # Angles: r1-l1-l2 (thetaB0, kthetaB)
        output.append(
            write_angle(('r1', 'l1', 'l2'), 'thetaB0', 'kthetaB'))

        output.append('[ dihedrals ]')
        output.append(
            '; ai         aj         ak         al      type phiA       fcA        phiB       fcB')
        # Dihedrals: r3-r2-r1-l1 (phiA0, kphiA)
        output.append(
            write_dihedral(('r3', 'r2', 'r1', 'l1'), 'phiA0', 'kphiA'))
        # Dihedrals: r2-r1-l1-l2 (phiB0, kphiB)
        output.append(
            write_dihedral(('r2', 'r1', 'l1', 'l2'), 'phiB0', 'kphiB'))
        # Dihedrals: r1-l1-l2-l3 (phiC0, kphiC)
        output.append(
            write_dihedral(('r1', 'l1', 'l2', 'l3'), 'phiC0', 'kphiC'))

        return '\n'.join(output)

    def toString(self, engine='Gromacs'):
        """The method for convert the restraint to a format that could be used
        by MD Engines.

           Parameters
           ----------

           engine : str
               The molecular dynamics engine used to generate the restraint.
               Available options currently is "GROMACS" only. If this argument
               is omitted then BioSimSpace will choose an appropriate engine
               for you.
        """
        if engine.strip().lower() == 'gromacs':
            if self._restraint_type == 'boresch':
                return self._gromacs_boresch()
            else:
                raise NotImplementedError(
                    f'Restraint type {self.restraint_type} not implemented '
                    f'yet. Only boresch restraint is supported.')
        else:
            raise NotImplementedError(f'MD Engine {engine} not implemented '
                                      f'yet. Only Gromacs is supported.')

    @property
    def correction(self):
        '''Give the free energy of removing the restraint.'''
        if self._restraint_type == 'boresch':
            K = _k_boltz * (_kcal_per_mol / _kelvin) / (_kj_per_mol / _kelvin) # Gas constant in kJ/mol/K
            V = ((_meter3 / 1000 / _mole) / _nanometer3).value() # standard volume in nm^3 (liter/N_A)

            T = self.T / _kelvin # Temperature in Kelvin
            r0 = self._restraint_dict['equilibrium_values']['r0'] / _nanometer # Distance in nm
            thA = self._restraint_dict['equilibrium_values']['thetaA0'] / _radian # Angle in radians
            thB = self._restraint_dict['equilibrium_values']['thetaB0'] / _radian  # Angle in radians

            K_r = self._restraint_dict['force_constants']['kr'] / (_kj_per_mol / _nanometer2) # force constant for distance (kJ/mol/nm^2)
            K_thA = self._restraint_dict['force_constants']['kthetaA'] / (_kj_per_mol / (_radian * _radian))  # force constant for angle (kJ/mol/rad^2)
            K_thB = self._restraint_dict['force_constants']['kthetaB'] / (_kj_per_mol / (_radian * _radian))  # force constant for angle (kJ/mol/rad^2)
            K_phiA = self._restraint_dict['force_constants']['kphiA'] / (_kj_per_mol / (_radian * _radian))  # force constant for dihedral (kJ/mol/rad^2)
            K_phiB = self._restraint_dict['force_constants']['kphiB'] / (_kj_per_mol / (_radian * _radian))  # force constant for dihedral (kJ/mol/rad^2)
            K_phiC = self._restraint_dict['force_constants']['kphiC'] / (_kj_per_mol / (_radian * _radian))  # force constant for dihedral (kJ/mol/rad^2)

            # Convert all the units to float before this calculation as BSS cannot handle root
            arg = ((8.0 * np.pi ** 2 * V) /
                   (r0 ** 2 * np.sin(thA) * np.sin(thB)) *
                   (((K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC) ** 0.5) /
                    ((2.0 * np.pi * K * T) ** 3)))

            dG = - K * T * np.log(arg)
            # Attach unit
            dG *= _kj_per_mol
            return dG
