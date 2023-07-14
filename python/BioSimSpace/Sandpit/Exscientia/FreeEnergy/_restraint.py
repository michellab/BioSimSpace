######################################################################
# BioSimSpace: Making biomolecular simulation a breeze!
#
# Copyright: 2017-2023
#
# Authors: Lester Hedges <lester.hedges@gmail.com>
#
# BioSimSpace is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
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
from math import e
from scipy import integrate as _integrate
import numpy as _np
import warnings as _warnings

from sire.legacy.Units import angstrom3 as _Sire_angstrom3
from sire.legacy.Units import k_boltz as _k_boltz
from sire.legacy.Units import meter3 as _Sire_meter3
from sire.legacy.Units import mole as _Sire_mole

from .._SireWrappers import Atom as _Atom
from .._SireWrappers import System as _System
from ..Types import Angle as _Angle
from ..Types import Length as _Length
from ..Types import Temperature as _Temperature
from ..Units.Angle import degree as _degree
from ..Units.Angle import radian as _radian
from ..Units.Area import angstrom2 as _angstrom2
from ..Units.Energy import kcal_per_mol as _kcal_per_mol
from ..Units.Energy import kj_per_mol as _kj_per_mol
from ..Units.Length import angstrom as _angstrom
from ..Units.Length import nanometer as _nanometer
from ..Units.Temperature import kelvin as _kelvin


class Restraint:
    """The Restraint class which holds the restraint information for the ABFE
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
    """

    def __init__(self, system, restraint_dict, temperature, restraint_type="Boresch"):
        """
        Constructor.

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
                "temperature must be of type 'BioSimSpace.Types.Temperature'"
            )
        else:
            self.T = temperature

        if restraint_type.lower() == "boresch":
            self._restraint_type = "boresch"
            # Test if the atoms are of BioSimSpace._SireWrappers.Atom
            for key in ["r3", "r2", "r1", "l1", "l2", "l3"]:
                if not isinstance(restraint_dict["anchor_points"][key], _Atom):
                    raise ValueError(
                        f"restraint_dict['anchor_points']['{key}'] "
                        "must be of type "
                        "'BioSimSpace._SireWrappers.Atom'"
                    )

            # Test if the equilibrium length of the bond r1-l1 is a length unit
            # Such as angstrom or nanometer
            if not isinstance(restraint_dict["equilibrium_values"]["r0"], _Length):
                raise ValueError(
                    "restraint_dict['equilibrium_values']['r0'] must be of type 'BioSimSpace.Types.Length'"
                )

            # Test if the equilibrium length of the angle and dihedral is a
            # angle unit such as radian or degree
            for key in ["thetaA0", "thetaB0", "phiA0", "phiB0", "phiC0"]:
                if not isinstance(restraint_dict["equilibrium_values"][key], _Angle):
                    raise ValueError(
                        f"restraint_dict['equilibrium_values']['{key}'] must be "
                        f"of type 'BioSimSpace.Types.Angle'"
                    )

            # Test if the force constant of the angle and dihedral is the correct unit
            # Such as kcal/mol/rad^2
            for key in ["kthetaA", "kthetaB", "kphiA", "kphiB", "kphiC"]:
                if restraint_dict["force_constants"][key] != 0:
                    dim = restraint_dict["force_constants"][key].dimensions()
                    if dim != (-2, 0, 2, 1, -1, 0, -2):
                        raise ValueError(
                            f"restraint_dict['force_constants']['{key}'] must be of type "
                            f"'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Angle^2'"
                        )

            # Test for unstable combinations of force constants
            non_zero_force_const = [
                i[0] for i in restraint_dict["force_constants"].items() if i[1] != 0
            ]
            if "kr" not in non_zero_force_const:
                raise ValueError('"kr" cannot be zero')
            if "kthetaA" not in non_zero_force_const:
                if "kphiA" in non_zero_force_const or "kphiB" in non_zero_force_const:
                    raise ValueError(
                        "Restraining phiA or phiB without restraining thetaA "
                        "will produce unstable Boresch restraints."
                    )
            if "kthetaB" not in non_zero_force_const:
                if "kphiB" in non_zero_force_const or "kphiC" in non_zero_force_const:
                    raise ValueError(
                        "Restraining phiB or phiC without restraining thetaB "
                        "will produce unstable Boresch restraints."
                    )

            # Test if the force constant of the bond r1-l1 is the correct unit
            # Such as kcal/mol/angstrom^2
            dim = restraint_dict["force_constants"]["kr"].dimensions()
            if dim != (0, 0, 0, 1, -1, 0, -2):
                raise ValueError(
                    "restraint_dict['force_constants']['kr'] must be of type "
                    "'BioSimSpace.Types.Energy'/'BioSimSpace.Types.Length^2'"
                )

            # Ensure restrained angles are >= 10 kT from collinear
            R = (
                _k_boltz.value() * _kcal_per_mol / _kelvin
            ).value()  # molar gas constant in kcal mol-1 K-1
            T = self.T / _kelvin  # Temperature in Kelvin

            for angle in ["thetaA", "thetaB"]:
                force_const = restraint_dict["force_constants"][f"k{angle}"] / (
                    _kcal_per_mol / (_radian * _radian)
                )
                if force_const != 0:
                    equil_val = (
                        restraint_dict["equilibrium_values"][f"{angle}0"] / _radian
                    )

                    # Convert 10 kT to angle
                    R = (
                        _k_boltz.value() * _kcal_per_mol / _kelvin
                    ).value()  # molar gas constant in kcal mol-1 K-1
                    T = self.T / _kelvin  # Temperature in Kelvin
                    min_stable_dist = _np.sqrt((20 * R * T) / force_const)
                    min_dist = min([abs(equil_val - 0), abs(equil_val - _np.pi)])

                    if min_dist < min_stable_dist:
                        _warnings.warn(
                            f"The equilibrium value of {angle} is within 10 kT of"
                            "collinearity, which may result in unstable Boresch restraints."
                            " Consider increasing the force constants or selecting equilibrium"
                            " values further from 0 or pi radians."
                        )

        else:
            raise NotImplementedError(
                f"Restraint type {type} not implemented "
                f"yet. Only Boresch restraint is supported."
            )

        self._restraint_dict = restraint_dict
        self.system = system

    @property
    def system(self):
        return self._system

    @system.setter
    def system(self, system):
        """
        Update the system object.

        Parameters
        ----------

        system : :class:`System <BioSimSpace._SireWrappers.System>`
            The molecular system.
        """
        if not isinstance(system, _System):
            raise TypeError(
                "'system' must be of type 'BioSimSpace._SireWrappers.System'"
            )
        else:
            if self._restraint_type == "boresch":
                # Check if the ligand atoms are decoupled.
                # Find the decoupled molecule, assume that only one can be
                # decoupled.
                (decoupled_mol,) = system.getDecoupledMolecules()
                for key in ["l1", "l2", "l3"]:
                    atom = self._restraint_dict["anchor_points"][key]
                    # Discussed in https://github.com/michellab/BioSimSpace/pull/337
                    if (
                        atom._sire_object.molecule().number()
                        != decoupled_mol._sire_object.number()
                    ):
                        raise ValueError(
                            f"The ligand atom {key} is not from decoupled moleucle."
                        )
                for key in ["r1", "r2", "r3"]:
                    atom = self._restraint_dict["anchor_points"][key]
                    if not atom in system:
                        raise ValueError(
                            f"The protein atom {key} is not in the system."
                        )

            # Store a copy of solvated system.
            self._system = system.copy()

    def _gromacs_boresch(self):
        """Format the Gromacs string for boresch restraint."""

        # Format the atoms into index list
        def format_index(key_list):
            formated_index = []
            for key in key_list:
                formated_index.append(
                    "{:<10}".format(
                        self._system.getIndex(
                            self._restraint_dict["anchor_points"][key]
                        )
                        + 1
                    )
                )
            return " ".join(formated_index)

        parameters_string = "{eq0:<10} {fc0:<10} {eq1:<10} {fc1:<10}"

        # Format the parameters for the bonds
        def format_bond(equilibrium_values, force_constants):
            converted_equ_val = (
                self._restraint_dict["equilibrium_values"][equilibrium_values]
                / _nanometer
            )
            converted_fc = self._restraint_dict["force_constants"][force_constants] / (
                _kj_per_mol / _nanometer**2
            )
            return parameters_string.format(
                eq0="{:.3f}".format(converted_equ_val),
                fc0="{:.2f}".format(0),
                eq1="{:.3f}".format(converted_equ_val),
                fc1="{:.2f}".format(converted_fc),
            )

        # Format the parameters for the angles and dihedrals
        def format_angle(equilibrium_values, force_constants):
            converted_equ_val = (
                self._restraint_dict["equilibrium_values"][equilibrium_values] / _degree
            )
            converted_fc = self._restraint_dict["force_constants"][force_constants] / (
                _kj_per_mol / (_radian * _radian)
            )
            return parameters_string.format(
                eq0="{:.3f}".format(converted_equ_val),
                fc0="{:.2f}".format(0),
                eq1="{:.3f}".format(converted_equ_val),
                fc1="{:.2f}".format(converted_fc),
            )

        # basic format of the Gromacs string
        master_string = "  {index} {func_type} {parameters}"

        def write_bond(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=6,
                parameters=format_bond(equilibrium_values, force_constants),
            )

        def write_angle(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=1,
                parameters=format_angle(equilibrium_values, force_constants),
            )

        def write_dihedral(key_list, equilibrium_values, force_constants):
            return master_string.format(
                index=format_index(key_list),
                func_type=2,
                parameters=format_angle(equilibrium_values, force_constants),
            )

        # Writing the string
        output = [
            "[ intermolecular_interactions ]",
        ]

        output.append("[ bonds ]")
        output.append("; ai         aj      type bA         kA         bB         kB")
        # Bonds: r1-l1 (r0, kr)
        output.append(write_bond(("r1", "l1"), "r0", "kr"))

        output.append("[ angles ]")
        output.append(
            "; ai         aj         ak      type thA        fcA        thB        fcB"
        )
        # Angles: r2-r1-l1 (thetaA0, kthetaA)
        output.append(write_angle(("r2", "r1", "l1"), "thetaA0", "kthetaA"))
        # Angles: r1-l1-l2 (thetaB0, kthetaB)
        output.append(write_angle(("r1", "l1", "l2"), "thetaB0", "kthetaB"))

        output.append("[ dihedrals ]")
        output.append(
            "; ai         aj         ak         al      type phiA       fcA        phiB       fcB"
        )
        # Dihedrals: r3-r2-r1-l1 (phiA0, kphiA)
        output.append(write_dihedral(("r3", "r2", "r1", "l1"), "phiA0", "kphiA"))
        # Dihedrals: r2-r1-l1-l2 (phiB0, kphiB)
        output.append(write_dihedral(("r2", "r1", "l1", "l2"), "phiB0", "kphiB"))
        # Dihedrals: r1-l1-l2-l3 (phiC0, kphiC)
        output.append(write_dihedral(("r1", "l1", "l2", "l3"), "phiC0", "kphiC"))

        return "\n".join(output)

    def _somd_boresch(self):
        """Format the SOMD string for the Boresch restraints."""

        # Indices
        r1 = self._system.getIndex(self._restraint_dict["anchor_points"]["r1"])
        r2 = self._system.getIndex(self._restraint_dict["anchor_points"]["r2"])
        r3 = self._system.getIndex(self._restraint_dict["anchor_points"]["r3"])
        l1 = self._system.getIndex(self._restraint_dict["anchor_points"]["l1"])
        l2 = self._system.getIndex(self._restraint_dict["anchor_points"]["l2"])
        l3 = self._system.getIndex(self._restraint_dict["anchor_points"]["l3"])
        # Equilibrium values
        r0 = self._restraint_dict["equilibrium_values"]["r0"] / _angstrom
        thetaA0 = self._restraint_dict["equilibrium_values"]["thetaA0"] / _radian
        thetaB0 = self._restraint_dict["equilibrium_values"]["thetaB0"] / _radian
        phiA0 = self._restraint_dict["equilibrium_values"]["phiA0"] / _radian
        phiB0 = self._restraint_dict["equilibrium_values"]["phiB0"] / _radian
        phiC0 = self._restraint_dict["equilibrium_values"]["phiC0"] / _radian
        # Force constants. Halve these as SOMD defines force constants as E = kx**2
        kr = (
            self._restraint_dict["force_constants"]["kr"]
            / (_kcal_per_mol / (_angstrom * _angstrom))
        ) / 2
        kthetaA = (
            self._restraint_dict["force_constants"]["kthetaA"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kthetaB = (
            self._restraint_dict["force_constants"]["kthetaB"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kphiA = (
            self._restraint_dict["force_constants"]["kphiA"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kphiB = (
            self._restraint_dict["force_constants"]["kphiB"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2
        kphiC = (
            self._restraint_dict["force_constants"]["kphiC"]
            / (_kcal_per_mol / (_radian * _radian))
        ) / 2

        restr_string = f'boresch restraints dictionary = {{"anchor_points":{{"r1":{r1}, "r2":{r2}, "r3":{r3}, "l1":{l1}, '
        restr_string += f'"l2":{l2}, "l3":{l3}}}, '
        restr_string += f'"equilibrium_values":{{"r0":{r0:.2f}, "thetaA0":{thetaA0:.2f}, "thetaB0":{thetaB0:.2f},"phiA0":{phiA0:.2f}, '
        restr_string += f'"phiB0":{phiB0:.2f}, "phiC0":{phiC0:.2f}}}, '
        restr_string += f'"force_constants":{{"kr":{kr:.2f}, "kthetaA":{kthetaA:.2f}, "kthetaB":{kthetaB:.2f}, "kphiA":{kphiA:.2f}, '
        restr_string += f'"kphiB":{kphiB:.2f}, "kphiC":{kphiC:.2f}}}}}'

        return restr_string

    def toString(self, engine=None):
        """
        The method for convert the restraint to a format that could be used
        by MD Engines.

        Parameters
        ----------

        engine : str
            The molecular dynamics engine used to generate the restraint.
            Available options currently is "GROMACS" and "SOMD". If this argument
            is omitted then BioSimSpace will choose an appropriate engine
            for you.
        """
        if engine.strip().lower() == "gromacs":
            if self._restraint_type == "boresch":
                return self._gromacs_boresch()
        elif engine.lower() == "somd":
            if self._restraint_type == "boresch":
                return self._somd_boresch()
            else:
                raise NotImplementedError(
                    f"Restraint type {self._restraint_type} not implemented "
                    f"yet. Only Boresch restraints are supported."
                )
        else:
            raise NotImplementedError(
                f"MD Engine {engine} not implemented "
                f"yet. Only Gromacs and SOMD are supported."
            )

    def getCorrection(self, method="analytical"):
        """
        Calculate the free energy of releasing the restraint
        to the standard state volume.'''

        Parameters
        ----------

        method : str
            The integration method to use for calculating the correction for
            releasing the restraint to the standard state concentration.
            "numerical" or "analytical". Note that the analytical correction
            can introduce errors when the restraints are weak, restrained
            angles are close to 0 or pi radians, or the restrained distance
            is close to 0.

        Returns
        ----------
        dG : float
            Free energy of releasing the restraint to the standard state volume,
            in kcal / mol.
        """

        if self._restraint_type == "boresch":
            # Constants. Take .value() to avoid issues with ** and log of GeneralUnit
            v0 = (
                ((_Sire_meter3 / 1000) / _Sire_mole) / _Sire_angstrom3
            ).value()  # standard state volume in A^3
            R = (
                _k_boltz.value() * _kcal_per_mol / _kelvin
            ).value()  # molar gas constant in kcal mol-1 K-1

            # Parameters
            T = self.T / _kelvin  # Temperature in Kelvin
            prefactor = (
                8 * (_np.pi**2) * v0
            )  # In A^3. Divide this to account for force constants of 0 in the
            # analytical correction

            if method == "numerical":
                # ========= Acknowledgement ===============
                # Calculation copied from restraints.py  in
                # Yank https://github.com/choderalab/yank
                # =========================================

                def numerical_distance_integrand(r, r0, kr):
                    """Integrand for harmonic distance restraint. Domain is on [0, infinity],
                    but this will be truncated to [0, 8 RT] for practicality.

                    Parameters
                    ----------
                        r : (float)
                            Distance to be integrated, in Angstrom
                        r0 : (float)
                            Equilibrium distance, in Angstrom
                        kr : (float)
                            Force constant, in kcal mol-1 A-2

                    Returns
                    ----------
                        float : Value of integrand
                    """
                    return (r**2) * _np.exp(-(kr * (r - r0) ** 2) / (2 * R * T))

                def numerical_angle_integrand(theta, theta0, ktheta):
                    """Integrand for harmonic angle restraints. Domain is on [0,pi].

                    Parameters
                    ----------
                        theta : (float)
                        Angle to be integrated, in radians
                        theta0 : (float)
                        Equilibrium angle, in radians
                        ktheta : (float)
                        Force constant, in kcal mol-1 rad-2

                    Returns
                    ----------
                        float: Value of integrand
                    """
                    return _np.sin(theta) * _np.exp(
                        -(ktheta * (theta - theta0) ** 2) / (2 * R * T)
                    )

                def numerical_dihedral_integrand(phi, phi0, kphi):
                    """Integrand for the harmonic dihedral restraints. Domain is on [-pi,pi].

                    Parameters
                    ----------
                        phi : (float)
                        Angle to be integrated, in radians
                        phi0 : (float)
                        Equilibrium angle, in radians
                        kphi : (float)
                        Force constant, in kcal mol-1 rad-2

                    Returns
                    ----------
                        float: Value of integrand
                    """
                    d_phi = abs(phi - phi0)
                    d_phi_corrected = min(
                        d_phi, 2 * _np.pi - d_phi
                    )  # correct for periodic boundaries
                    return _np.exp(-(kphi * d_phi_corrected**2) / (2 * R * T))

                # Radial
                r0 = self._restraint_dict["equilibrium_values"]["r0"] / _angstrom  # A
                kr = self._restraint_dict["force_constants"]["kr"] / (
                    _kcal_per_mol / _angstrom2
                )  # kcal mol-1 A-2
                dist_at_8RT = 4 * _np.sqrt(
                    (R * T) / kr
                )  # Dist. which gives restraint energy = 8 RT
                r_min = max(0, r0 - dist_at_8RT)
                r_max = r0 + dist_at_8RT
                integrand = lambda r: numerical_distance_integrand(r, r0, kr)
                z_r = _integrate.quad(integrand, r_min, r_max)[0]

                # Angular
                for angle in ["thetaA", "thetaB"]:
                    theta0 = (
                        self._restraint_dict["equilibrium_values"][f"{angle}0"]
                        / _radian
                    )  # rad
                    ktheta = self._restraint_dict["force_constants"][f"k{angle}"] / (
                        _kcal_per_mol / (_radian * _radian)
                    )  # kcal mol-1 rad-2
                    integrand = lambda theta: numerical_angle_integrand(
                        theta, theta0, ktheta
                    )
                    z_r *= _integrate.quad(integrand, 0, _np.pi)[0]

                # Dihedral
                for dihedral in ["phiA", "phiB", "phiC"]:
                    phi0 = (
                        self._restraint_dict["equilibrium_values"][f"{dihedral}0"]
                        / _radian
                    )  # rad
                    kphi = self._restraint_dict["force_constants"][f"k{dihedral}"] / (
                        _kcal_per_mol / (_radian * _radian)
                    )  # kcal mol-1 rad-2
                    integrand = lambda phi: numerical_dihedral_integrand(
                        phi, phi0, kphi
                    )
                    z_r *= _integrate.quad(integrand, -_np.pi, _np.pi)[0]

                # Compute dg and attach unit
                dg = -R * T * _np.log(prefactor / z_r)
                dg *= _kcal_per_mol

                return dg

            elif method == "analytical":
                # Only need three equilibrium values for the analytical correction
                r0 = (
                    self._restraint_dict["equilibrium_values"]["r0"] / _angstrom
                )  # Distance in A
                thetaA0 = (
                    self._restraint_dict["equilibrium_values"]["thetaA0"] / _radian
                )  # Angle in radians
                thetaB0 = (
                    self._restraint_dict["equilibrium_values"]["thetaB0"] / _radian
                )  # Angle in radians

                force_constants = []

                # Loop through and correct for force constants of zero,
                # which break the analytical correction. To account for this,
                # divide the prefactor accordingly. Note that setting
                # certain force constants to zero while others are non-zero
                # will result in unstable restraints, but this will be checked when
                # the restraint object is created
                for k, val in self._restraint_dict["force_constants"].items():
                    if val.value() == 0:
                        if k == "kr":
                            raise ValueError("The force constant kr must not be zero")
                        if k == "kthetaA":
                            prefactor /= 2 / _np.sin(thetaA0)
                        if k == "kthetaB":
                            prefactor /= 2 / _np.sin(thetaB0)
                        if k[:4] == "kphi":
                            prefactor /= 2 * _np.pi
                    else:
                        if k == "kr":
                            force_constants.append(val / (_kcal_per_mol / _angstrom2))
                        else:
                            force_constants.append(
                                val / (_kcal_per_mol / (_radian * _radian))
                            )

                # Calculation
                n_nonzero_k = len(force_constants)
                prod_force_constants = _np.prod(force_constants)
                numerator = prefactor * _np.sqrt(prod_force_constants)
                denominator = (
                    (r0**2)
                    * _np.sin(thetaA0)
                    * _np.sin(thetaB0)
                    * (2 * _np.pi * R * T) ** (n_nonzero_k / 2)
                )

                # Compute dg and attach unit
                dg = -R * T * _np.log(numerator / denominator)
                dg *= _kcal_per_mol

                return dg

            else:
                raise ValueError(
                    f"Correction method {method} is not supported. Please choose from 'numerical' or 'analytical'."
                )

    @property
    def correction(self):
        """Give the free energy of removing the restraint."""
        return self.getCorrection()
