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
Functionality for aligning molecules.
"""

__author__ = "Lester Hedges"
__email__ = "lester.hedges@gmail.com"

__all__ = ["generateNetwork",
           "matchAtoms",
           "viewMapping",
           "rmsdAlign",
           "flexAlign",
           "merge"]

import csv as _csv
import os as _os
import subprocess as _subprocess
import shutil as _shutil
import shlex as _shlex
import sys as _sys
import tempfile as _tempfile

from .._Utils import _try_import, _have_imported, _assert_imported

import warnings as _warnings
# Suppress numpy warnings from RDKit import.
_warnings.filterwarnings("ignore", message="numpy.dtype size changed")
_warnings.filterwarnings("ignore", message="numpy.ndarray size changed")
_warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

# Suppress duplicate to-Python converted warnings.
# Both Sire and RDKit register the same converter.
with _warnings.catch_warnings():
    _warnings.filterwarnings("ignore")
    _rdkit = _try_import("rdkit")

    if _have_imported(_rdkit):
        from rdkit import Chem as _Chem
        from rdkit.Chem import rdFMCS as _rdFMCS
        from rdkit import RDLogger as _RDLogger

        # Disable RDKit warnings.
        _RDLogger.DisableLog('rdApp.*')
    else:
        _Chem = _rdkit
        _rdFMCS = _rdkit
        _RDLogger = _rdkit

from Sire import Base as _SireBase
from Sire import Maths as _SireMaths
from Sire import Mol as _SireMol
from Sire import Units as _SireUnits

from .. import _is_notebook, _isVerbose
from .._Exceptions import AlignmentError as _AlignmentError
from .._Exceptions import MissingSoftwareError as _MissingSoftwareError
from .._SireWrappers import Molecule as _Molecule

from .. import IO as _IO
from .. import Units as _Units
from .. import _Utils

# lomap depends on RDKit and networkx
_networkx = _try_import("networkx")

if _have_imported(_rdkit) and _have_imported(_networkx):
    import lomap as _lomap
elif _have_imported(_rdkit):
    _lomap = _networkx
elif _have_imported(_networkx):
    _lomap = _rdkit
else:
    from .._Utils import _module_stub
    _lomap = _module_stub(name="rdkit, networkx")

from ._merge import merge as _merge

# Try to find the FKCOMBU program from KCOMBU: https://pdbj.org/kcombu
try:
    _fkcombu_exe = _SireBase.findExe("fkcombu").absoluteFilePath()
except:
    _fkcombu_exe = None

def generateNetwork(molecules, names=None, work_dir=None, plot_network=False,
        links_file=None, property_map={}, n_edges_forced=None):
    """Generate a perturbation network using Lead Optimisation Mappper (LOMAP).

       Parameters
       ----------

       molecules : :[class:`Molecule <BioSimSpace._SireWrappers.Molecule>`], \
                    [rdkit.Chem.rdchem.Mol]
           A list of molecules. (Both BioSimSpace and RDKit molecule objects
           are supported.)

       names : [str]
           A list of names for the molecules. If None, then the index of each
           molecule will be used.

       work_dir : str
           The working directory for the LOMAP process.

       plot_network : bool
           Whether to plot the network when running from within a notebook.
           If using a 'work_dir', then a PNG image will be located in
           'work_dir/images/network.png'.

       links_file : str
           Path to a file providing links to seed the LOMAP graph with. Each
           record in the file must contain a minimum of two entries, i.e. the
           names of the ligands of interest, e.g.
             ligA ligB
           The third column can include a pre-computed score for the ligand
           pair, e.g.
             ligA ligB score
           A value of < -1 means 'recompute', but force the link to be included.
           A negative value (>= -1) means use the absolute value as the score,
           but force the link to be included. A positive means use this value
           as the score, but treat the link as normal in the LOMAP graph
           calculation. Finally, if the fourth column contains the string
           'force', then the link is included, irrespective of its score.

       property_map : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       n_edges_forced : int
           An integer that forces the number of edges that should be used in
           the perturbation network. Must be in the range
           [1 .. (len(molecules)**2-len(molecules))/2].
           In cases where n_edges_forced > the number of edges suggested by
           LOMAP, BioSimSpace will add the top scoring n edges parsed from the
           LOMAP output file and add them to the network. Conversely if
           n_edges_forced < the number of edges suggested by LOMAP, BioSimSpace
           will remove the bottom n edges parsed from the LOMAP output file.
           This last option is discouraged as it can cause network cycle
           breakage and disconnecting of ligands/clusters from the network.

       Returns
       -------

       edges : [(int, int)]
           A tuple containing the edges of the network. Each edge is itself
           a tuple, containing the indices of the molecules that are
           connected by the edge.

       scores : [float]
           The LOMAP score for each edge. The higher the score implies that a
           perturbation between molecules along an edge is likely to be more
           accurate.
    """

    # Adapted from code by Jenke Scheen (@JenkeScheen).

    _assert_imported(_lomap)

    if not isinstance(molecules, (list, tuple)):
        raise TypeError("'molecules' must be a list of "
                        "'BioSimSpace._SireWrappers.Molecule' "
                        "or 'rdkit.Chem.rdchem.Mol' objects.")

    # Validate the molecules.
    rdkit_input = False

    # A list of BioSimSpace molecule objects.
    if all(isinstance(x, _Molecule) for x in molecules):
        pass
    # A list of RDKit molecule objects.
    elif all(isinstance(x, _Chem.rdchem.Mol) for x in molecules):
        rdkit_input = True
    else:
        raise TypeError("'molecules' must be a list of "
                        "'BioSimSpace._SireWrappers.Molecule' "
                        "or 'rdkit.Chem.rdchem.Mol' objects.")

    # Validate the names.
    if names is not None:
        if not all(isinstance(x, str) for x in names):
            raise TypeError("'names' must be a list of 'str' types.")
        if len(names) != len(molecules):
            raise ValueError("There must be one name for each molecule!")

    # Validate the working directory.
    if work_dir is not None:
        if not isinstance(work_dir, str):
            raise TypeError("'work_dir' must be of type 'str'.")

    # Validate the plotting flag.
    if not isinstance(plot_network, bool):
        raise TypeError("'plot_network' must be of type 'bool'.")

    # Validate the property map.
    if not isinstance(property_map, dict):
        raise TypeError("'property_map' must be of type 'dict'")

    # Validate the scores file.
    if links_file is not None:

        if not isinstance(links_file, str):
            raise TypeError("'links_file' must be of type 'str'")

        # Check that it exists.
        if _os.path.isfile(links_file):
            # Must have names to match against.
            if names is None:
                raise ValueError("'names' must be defined when passing 'links_file' to LOMAP.")

            # Validate the records in the links file. Must have at least two
            # entries per line, i.e. ligA ligB, third column must contain a
            # float type score if present, and fourth column can only contain
            # 'force' to insist that the link is included irrespective of its
            # score.
            with open(links_file, "r") as lf:
                for line in lf:
                    records = line.split()
                    if len(records) < 2:
                        raise ValueError(f"{links_file} should have at least two entries "
                                         f"per line. Failed line is {line.rstrip()}")
                    if len(records) > 2:
                        try:
                            float(records[2])
                        except ValueError:
                            raise ValueError(f"{links_file} contains a non-numerical value for the "
                                             f"score in the third column: {line.rstrip()}.")
                    if len(records) > 3:
                        if records[3] != "force":
                            raise ValueError(f"{links_file} can only contain 'force' in the "
                                             f"fourth column: {line.rstrip()}.")

                    # Make sure that the ligands are in the names list.
                    if not records[0] in names:
                        raise ValueError(f"Ligand '{records[0]}' not in 'names' list!")
                    if not records[1] in names:
                        raise ValueError(f"Ligand '{records[1]}' not in 'names' list!")

        else:
            raise IOError(f"The links file doesn't exist: {links_file}")

    # Validate the number of edges parameter.
    if n_edges_forced is not None:
        if not type(n_edges_forced) is int:
            raise TypeError("'n_edges_forced' must be of type 'int'")

        n_edges_fully_connected = int((len(molecules)**2 - len(molecules))/2)+1

        if not 0 < n_edges_forced < n_edges_fully_connected:
            raise ValueError(f"'n_edges_forced' must be 0 < value < {n_edges_fully_connected}.")

    # Create a temporary working directory and store the directory name.
    if work_dir is None:
        tmp_dir = _tempfile.TemporaryDirectory()
        work_dir = tmp_dir.name

    # User specified working directory.
    else:
        # Use full path.
        if work_dir[0] != "/":
            work_dir = _os.getcwd() + "/" + work_dir

        # Create the directory if it doesn't already exist.
        if not _os.path.isdir(work_dir):
            _os.makedirs(work_dir, exist_ok=True)

    # Make the LOMAP input and output directories.
    _os.makedirs(work_dir + "/inputs", exist_ok=True)
    _os.makedirs(work_dir + "/outputs", exist_ok=True)

    # Dictionary to map ligand names in the links file to the file names
    # that are used in the working directory.
    links_names = {}

    # Write all of the molecules to disk.
    if rdkit_input:
        if names is not None:
            for x, (molecule, name) in enumerate(zip(molecules, names)):
                file_name = f"{x:03d}_{name}.sdf"
                links_names[name] = file_name
                writer =  _Chem.SDWriter(f"{work_dir}/inputs/{file_name}")
                writer.write(molecule)
                writer.close()
        else:
            for x, molecule in enumerate(molecules):
                writer =  _Chem.SDWriter(work_dir + f"/inputs/{x:03d}.sdf")
                writer.write(molecule)
                writer.close()
    else:
        if names is not None:
            is_names = True
            _names = names
        else:
            is_names = False
            # Create a dummy names list so that we can handle everything in a
            # single loop.
            _names = ["txt" for x in range(len(molecules))]

        for x, (molecule, name) in enumerate(zip(molecules, _names)):
            # If the molecule came from an SDF file, then use
            # that as the format as it's generally more reliable.
            is_sdf = False
            if molecule._sire_object.hasProperty("fileformat"):
                if "SDF" in molecule._sire_object.property("fileformat").value():
                    is_sdf = True
                    if is_names:
                        _IO.saveMolecules(work_dir + f"/inputs/{x:03d}_{name}",
                            molecule, "sdf", property_map=property_map)
                        links_names[name] = f"{x:03d}_{name}.sdf"
                    else:
                        _IO.saveMolecules(work_dir + f"/inputs/{x:03d}",
                            molecule, "sdf", property_map=property_map)
                else:
                    if is_names:
                        _IO.saveMolecules(work_dir + f"/inputs/{x:03d}_{name}",
                            molecule, "mol2", property_map=property_map)
                        links_names[name] = f"{x:03d}_{name}.mol2"
                    else:
                        _IO.saveMolecules(work_dir + f"/inputs/{x:03d}",
                            molecule, "mol2", property_map=property_map)
            else:
                if is_names:
                    _IO.saveMolecules(work_dir + f"/inputs/{x:03d}_{name}",
                        molecule,  "mol2", property_map=property_map)
                    links_names[name] = f"{x:03d}_{name}.mol2"
                else:
                    _IO.saveMolecules(work_dir + f"/inputs/{x:03d}",
                        molecule,  "mol2", property_map=property_map)

    # Create a local copy of the links file in the working directory,
    # replacing the original ligand names with their actual file names.
    if links_file:
        # Read the old file and map the ligand names to their file names.
        new_lines = []
        with open(links_file, "r") as lf:
            for line in lf:
                records = line.split()
                if len(records) > 2:
                    lf_data = " ".join(records[2:])
                new_lines.append((links_names[records[0]],links_names[records[1]],lf_data))
       
        # Write the updated lomap links file.
        with open(f"{work_dir}/inputs/lomap_links_file.txt", "w") as lf:
            writer = _csv.writer(lf, delimiter=" ")
            for line in new_lines:
                writer.writerow(line)        

        # specify the file path
        lf = f"{work_dir}/inputs/lomap_links_file.txt"
        
    else:
        lf = None

    # Create the DBMolecules object.
    db_mol = _lomap.DBMolecules(f"{work_dir}/inputs",
                                name=f"{work_dir}/outputs/lomap",
                                links_file=lf,
                                output=True,
                                output_no_graph=True,
                                output_no_images=True,
                                threed=False,
                                max3d=3.0,
                                time=3,
                                parallel=10)

    # Create the similarity matrices.
    strict, loose = db_mol.build_matrices()

    # Generate the graph. (Might be able to use this later, rather than
    # reconstructing it by hand.)
    nx_graph = db_mol.build_graph()

    # Store the name to the LOMAP output file.
    lomap_file = work_dir + "/outputs/lomap_score_with_connection.txt"

    # Check that it exists.
    if not _os.path.isfile(lomap_file):
        raise _AlignmentError("LOMAP output file doesn't exist!")

    # Read the file to get the edges and scores.
    edges = []
    nodes = []
    scores = []

    edges_excluded = []

    with open(lomap_file, "r") as csv_file:
        # Load as a CSV file.
        csv_reader = _csv.reader(csv_file)

        # Loop over all rows of the CSV file.
        for row in csv_reader:
            # If the file contains all possible edges, then only take edges
            # that LOMAP indicates should be drawn.
            if row[7].strip() == "Yes":
                # Extract the nodes (molecules) connected by the edge.
                mol0 = int(row[2].rsplit(".")[0].rsplit("_")[0])
                mol1 = int(row[3].rsplit(".")[0].rsplit("_")[0])

                # Extract the score and convert to a float.
                score = float(row[4])

                # Update the lists.
                edges.append((mol0, mol1))
                nodes.append(mol0)
                nodes.append(mol1)
                scores.append(score)

            # Also collect the excluded edges in case we need to add more at a later stage.
            elif row[7].strip() == "No":
                # Extract the nodes (molecules) connected by the edge.
                mol0 = int(row[2].rsplit(".")[0].rsplit("_")[0])
                mol1 = int(row[3].rsplit(".")[0].rsplit("_")[0])

                # Extract the score and convert to a float.
                score = float(row[4])

                # Update the list while checking that the inverse edge is not already in
                # the network.
                if not (mol1, mol0) in edges:
                    edges_excluded.append((mol0, mol1, score))

    # If the user has specified a forced number of edges, adjust the network
    # to match the query. We have three situations to deal with.
    if n_edges_forced:

        # sort the list of excluded edges by LOMAP score.
        edges_excluded.sort(key=lambda x: x[2], reverse=True)

        # 1) The network already contains the specified number of edges.
        if  len(edges) == n_edges_forced:
            # Return the network as is.
            print(f"LOMAP already suggested the user-specified number of edges ({len(edges)}).")

        # 2) The network contains fewer edges than the specified number.
        elif len(edges) < n_edges_forced:
            # We need to add edges to the network to match the queried number.
            n_to_add = n_edges_forced-len(edges)
            print(f"Adding {n_to_add} edges to the LOMAP network.")

            # Get the top n excluded edges.
            for edge in edges_excluded[:n_to_add]:
                edges.append((edge[0], edge[1]))
                nodes.append(edge[0])
                nodes.append(edge[1])
                scores.append(edge[2])

        # 3) The network contains more edges than the specified number. This is not
        # recommended as this can cause breaking of network cycles or disconnecting nodes.
        elif len(edges) > n_edges_forced:
            # We need to remove edges from the network to match the queried number.
            n_to_remove = len(edges) - n_edges_forced
            print(f"Removing {n_to_remove} edges from the LOMAP network, potentially" \
                +" breaking network cycles or disconnecting ligands/ clusters.")
            lomap_network = list(zip(edges, scores))

            # Sort the network by LOMAP-score.
            lomap_network.sort(key=lambda x: x[1])

            # Create the new network by keeping the required number of top edges.
            edges = []
            nodes = []
            scores = []
            for edge in lomap_network[n_to_remove:]:
                edges.append(edge[0])
                nodes.append(edge[0][0])
                nodes.append(edge[0][1])
                scores.append(edge[1])

    # Convert nodes to a set to remove duplicates.
    nodes = set(nodes)

    # Plot the LOMAP network.
    if plot_network:
        # Conditional imports.
        _assert_imported(_rdkit)

        import matplotlib.image as _mpimg
        import matplotlib.pyplot as _plt
        _nx = _try_import("networkx")
        _pydot = _try_import("pydot")

        _assert_imported(_nx)
        _assert_imported(_pydot)

        from rdkit.Chem import AllChem as _AllChem
        from rdkit.Chem import Draw as _Draw
        from rdkit.Chem import rdmolops as _rdmolops

        # Set the DPI to make the network look nice.
        _plt.rcParams["figure.dpi"]= 150

        # Make directory for output images.
        _os.makedirs(work_dir + "/images", exist_ok=True)

        # 1) Loop over each molecule and load into RDKit.
        try:
            rdmols = []
            if names is not None:
                is_names = True
                _names = names
            else:
                is_names = False
                # Create a dummy names list so that we can handle everything in a
                # single loop.
                _names = ["txt" for x in range(len(molecules))]

            for x, name in zip(range(0, len(molecules)), _names):
                if is_sdf:
                    ext = "sdf"
                else:
                    ext = "mol2"

                if is_names:
                    file = f"{work_dir}/inputs/{x:03d}_{name}.{ext}"
                else:
                    file = f"{work_dir}/inputs/{x:03d}.{ext}"

                if is_sdf:
                    rdmol = _Chem.MolFromMolFile(file, sanitize=False, removeHs=False)
                else:
                    rdmol = _Chem.MolFromMol2File(file, sanitize=False, removeHs=False)

                # RDKit doesn't thrown an exception, rather returns None and
                # prints an error message.
                if rdmol is None:
                    raise OSError(f"RDKit was unable to read: {file}")

                # Store the molecule.
                rdmols.append(rdmol)

        except Exception as e:
            msg = "RDKit was unable to load molecule!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _AlignmentError(msg) from e
            else:
                raise _AlignmentError(str(e)) from None

        # 2) Find the MCS of the molecules to use as a template.
        try:
            # Remove hydrogens to dramatically speed up MCS algorithm.
            rdmols = [_Chem.RemoveHs(mol) for mol in rdmols]

            template = _Chem.MolFromSmarts(_rdFMCS.FindMCS(rdmols,
                                    atomCompare=_rdFMCS.AtomCompare.CompareAny,
                                    bondCompare=_rdFMCS.BondCompare.CompareAny,
                                    matchValences=False,
                                    ringMatchesRingOnly=True,
                                    completeRingsOnly=True,
                                    matchChiralTag=False).smartsString)
            _AllChem.Compute2DCoords(template)

        except Exception as e:
            msg = "Unable to compute MCS of molecules!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _AlignmentError(msg) from e
            else:
                raise _AlignmentError(msg) from None

        # 3) Load all ligands, make 2D depiction aligned to the template and save to file.
        try:
            for x, mol in enumerate(rdmols):
                mol.UpdatePropertyCache(strict=False)
                _AllChem.Compute2DCoords(mol)
                _AllChem.GenerateDepictionMatching2DStructure(mol, template)

                mol = _Chem.RemoveHs(mol)

                # Remove stereochemistry to simplify depiction in network.
                _rdmolops.RemoveStereochemistry(mol)
                _Draw.MolToFile(mol, f"{work_dir}/images/{x:03d}.png")

        except Exception as e:
            msg = "Unable to make 2D depiction of molecules!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _AlignmentError(msg) from e
            else:
                raise _AlignmentError(msg) from None

        # 4) Create the NetworkX graph.
        try:
            # Generate the graph.
            graph = _nx.Graph()

            # If ligand names aren't specified, then use the molecule index.
            if names is None:
                names = [x for x in range(1, len(molecules)+1)]

            # Create a dictionary mapping the edges to their scores.
            edge_dict = {}
            for x, (node0, node1) in enumerate(edges):
                edge_dict[(names[node0], names[node1])] = str(round(scores[x], 2))

            # Loop over the nodes and add to the graph.
            for node in nodes:
                img = f"{work_dir}/images/{node:03d}.png"
                graph.add_node(names[node], image=img, label=names[node], labelloc="t")

            # Loop over the edges and add to the graph.
            for edge in edges:
                graph.add_edge(names[edge[0]],
                               names[edge[1]],
                               label=(edge_dict[(names[edge[0]], names[edge[1]])]))

        except Exception as e:
            msg = "Unable to generate network representation!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _AlignmentError(msg) from e
            else:
                raise _AlignmentError(msg) from None

        # 5) Create and display the plot.
        try:
            # Convert to a dot graph.
            dot_graph = _nx.drawing.nx_pydot.to_pydot(graph)

            # Write to a PNG.
            network_plot = f"{work_dir}/images/network.png"
            dot_graph.write_png(network_plot)

            if _is_notebook:
                # Create a plot of the network.
                img = _mpimg.imread(network_plot)
                _plt.figure(figsize=(20, 20))
                _plt.axis("off")
                _plt.imshow(img)

        except Exception as e:
            msg = "Unable to create network plot!"
            if _isVerbose():
                msg += ": " + getattr(e, "message", repr(e))
                raise _AlignmentError(msg) from e
            else:
                raise _AlignmentError(msg) from None

    return edges, scores

def matchAtoms(molecule0,
               molecule1,
               engine="SOMD",
               scoring_function="rmsd_align",
               matches=1,
               return_scores=False,
               prematch={},
               timeout=5*_Units.Time.second,
               complete_rings_only=True,
               prune_perturbed_constraints=None,
               prune_crossing_constraints=None,
               max_scoring_matches=1000,
               property_map0={},
               property_map1={}):
    """Find mappings between atom indices in molecule0 to those in molecule1.
       Molecules are aligned using a Maximum Common Substructure (MCS) search.
       When requesting more than one match, the mappings will be sorted using
       a scoring function and returned in order of best to worst score. (Note
       that, depending on the scoring function the "best" score may have the
       lowest value.)

       Parameters
       ----------

       molecule0 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule of interest.

       molecule1 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The reference molecule.

       engine : str
           The engine with which the simulations will be run. This is used for
           setting some reasonable defaults for the other settings.

       scoring_function : str
           The scoring function used to match atoms. Available options are:
             - "rmsd"
                 Calculate the root mean squared distance between the
                 coordinates of atoms in molecule0 to those that they
                 map to in molecule1.
             - "rmsd_align"
                 Align molecule0 to molecule1 based on the mapping before
                 computing the above RMSD score.
             - "rmsd_flex_align"
                 Flexibly align molecule0 to molecule1 based on the mapping
                 before computing the above RMSD score. (Requires the
                 'fkcombu'. package: https://pdbj.org/kcombu)

       matches : int
           The maximum number of matches to return. (Sorted in order of score).

       return_scores : bool
           Whether to return a list containing the scores for each mapping.

       prematch : dict
           A dictionary of atom mappings that must be included in the match.

       timeout : BioSimSpace.Types.Time
           The timeout for the maximum common substructure search.

       complete_rings_only : bool
           Whether to only match complete rings during the MCS search. This
           option is only relevant to MCS performed using RDKit and will be
           ignored when falling back on Sire.

       prune_perturbed_constraints : bool
           Whether to remove hydrogen atoms that are perturbed to heavy atoms
           from the mapping. This is True for AMBER by default and False for
           all other engines.

       prune_crossing_constraints : bool
           Whether to remove atoms from the mapping such that there are no
           constraints between dummy and non-dummy atoms. This is True for
           AMBER by default and False for all other engines.

       max_scoring_matches : int
           The maximum number of matching MCS substructures to consider when
           computing mapping scores. Consider reducing this if you find the
           matchAtoms function to be taking an excessive amount of time. This
           option is only relevant to MCS performed using RDKit and will be
           ignored when falling back on Sire.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       Returns
       -------

       matches : dict, [dict], ([dict], list)
           The best atom mapping, a list containing a user specified number of
           the best mappings ranked by their score, or a tuple containing the
           list of best mappings and a list of the corresponding scores.

       Examples
       --------

       Find the best maximum common substructure mapping between two molecules.

       >>> import BioSimSpace as BSS
       >>> mapping = BSS.Align.matchAtoms(molecule0, molecule1)

       Find the 5 best mappings.

       >>> import BioSimSpace as BSS
       >>> mappings = BSS.Align.matchAtoms(molecule0, molecule1, matches=5)

       Find the 5 best mappings along with their ranking scores.

       >>> import BioSimSpace as BSS
       >>> mappings, scores = BSS.Align.matchAtoms(molecule0, molecule1, matches=5, return_scores=True)

       Find the 5 best mappings along with their ranking scores. Score
       by flexibly aligning molecule0 to molecule1 based on each mapping
       and computing the root mean squared displacement of the matched
       atoms.

       >>> import BioSimSpace as BSS
       >>> mappings, scores = BSS.Align.matchAtoms(molecule0, molecule1, matches=5, return_scores=True, scoring_function="rmsd_flex_align")

       Find the best mapping that contains a prematch (this is a dictionary mapping
       atom indices in molecule0 to those in molecule1).

       >>> import BioSimSpace as BSS
       >>> mapping = BSS.Align.matchAtoms(molecule0, molecule1, prematch={0 : 10, 3 : 7})
    """

    # A list of supported scoring functions.
    scoring_functions = ["RMSD", "RMSDALIGN", "RMSDFLEXALIGN"]

    # Validate input.

    if not isinstance(molecule0, _Molecule):
        raise TypeError("'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(molecule1, _Molecule):
        raise TypeError("'molecule1' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(scoring_function, str):
        raise TypeError("'scoring_function' must be of type 'str'")
    else:
        # Strip underscores and whitespace, then convert to upper case.
        _scoring_function = scoring_function.replace("_", "").upper()
        _scoring_function = _scoring_function.replace(" ", "").upper()
        if not _scoring_function in scoring_functions:
            raise ValueError("Unsupported scoring function '%s'. Options are: %s"
                % (scoring_function, scoring_functions))

    if _scoring_function == "RMSDFLEXALIGN" and _fkcombu_exe is None:
        raise _MissingSoftwareError("'rmsd_flex_align' option requires the 'fkcombu' program: "
                                    "https://pdbj.org/kcombu")

    if not type(matches) is int:
        raise TypeError("'matches' must be of type 'int'")
    else:
        if matches < 0:
            raise ValueError("'matches' must be positive!")

    if not isinstance(return_scores, bool):
        raise TypeError("'return_matches' must be of type 'bool'")

    if not isinstance(prematch, dict):
        raise TypeError("'prematch' must be of type 'dict'")
    else:
        _validate_mapping(molecule0, molecule1, prematch, "prematch")

    if not isinstance(timeout, _Units.Time._Time):
        raise TypeError("'timeout' must be of type 'BioSimSpace.Types.Time'")

    if not type(max_scoring_matches) is int:
        raise TypeError("'max_scoring_matches' must be of type 'int'")

    if max_scoring_matches <= 0:
        raise ValueError("'max_scoring_matches' must be >= 1.")

    if not isinstance(property_map0, dict):
        raise TypeError("'property_map0' must be of type 'dict'")

    if not isinstance(property_map1, dict):
        raise TypeError("'property_map1' must be of type 'dict'")

    # Set some defaults.
    is_amber = engine.upper() == "AMBER"
    if prune_perturbed_constraints is None:
        prune_perturbed_constraints = is_amber
    if prune_crossing_constraints is None:
        prune_crossing_constraints = is_amber

    # Extract the Sire molecule from each BioSimSpace molecule.
    mol0 = molecule0._getSireObject()
    mol1 = molecule1._getSireObject()

    # Convert the timeout to seconds and take the value as an integer.
    timeout = int(timeout.seconds().value())

    # Create a temporary working directory.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Use RDKkit to find the maximum common substructure.

    try:
        # Run inside a temporary directory.
        with _Utils.cd(work_dir):
            # Write both molecules to PDB files.
            _IO.saveMolecules("tmp0", molecule0, "PDB", property_map=property_map0)
            _IO.saveMolecules("tmp1", molecule1, "PDB", property_map=property_map1)

            # Load the molecules with RDKit.
            # Note that the C++ function overloading seems to be broken, so we
            # need to pass all arguments by position, rather than keyword.
            # The arguments are: "filename", "sanitize", "removeHs", "flavor"
            mols = [_Chem.MolFromPDBFile("tmp0.pdb", True, False, 0),
                    _Chem.MolFromPDBFile("tmp1.pdb", True, False, 0)]

            # Generate the MCS match.
            mcs = _rdFMCS.FindMCS(mols,
                                  atomCompare=_rdFMCS.AtomCompare.CompareAny,
                                  bondCompare=_rdFMCS.BondCompare.CompareAny,
                                  completeRingsOnly=complete_rings_only,
                                  ringMatchesRingOnly=True,
                                  matchChiralTag=False,
                                  matchValences=False,
                                  maximizeBonds=False,
                                  timeout=timeout)

            # Get the common substructure as a SMARTS string.
            mcs_smarts = _Chem.MolFromSmarts(mcs.smartsString)

    except:
        raise RuntimeError("RDKIT MCS mapping failed!")

    # Score the mappings and return them in sorted order (best to worst).
    mappings, scores = _score_rdkit_mappings(mol0, mol1, mols[0], mols[1],
        mcs_smarts, prematch, _scoring_function, max_scoring_matches,
        property_map0, property_map1)

    # Sometimes RDKit fails to generate a mapping that includes the prematch.
    # If so, then try generating a mapping using the MCS routine from Sire.
    if len(mappings) == 1 and mappings[0] == prematch:

        # Warn that we've fallen back on using Sire.
        if prematch != {}:
            _warnings.warn("RDKit mapping didn't include prematch. Using Sire MCS.")

        # Warn about unsupported options.
        if not complete_rings_only:
            _warnings.warn("Using Sire MCS. Ignoring unsupported 'complete_rings_only' option!")

        # Convert timeout to a Sire Unit.
        timeout = timeout * _SireUnits.second

        # Regular match. Include light atoms, but don't allow matches between heavy
        # and light atoms.
        m0 = mol0.evaluate().findMCSmatches(mol1, _SireMol.AtomResultMatcher(_to_sire_mapping(prematch)),
                                            timeout, True, property_map0, property_map1, 6, False)

        # Include light atoms, and allow matches between heavy and light atoms.
        # This captures mappings such as O --> H in methane to methanol.
        m1 = mol0.evaluate().findMCSmatches(mol1, _SireMol.AtomResultMatcher(_to_sire_mapping(prematch)),
                                            timeout, True, property_map0, property_map1, 0, False)

        # Take the mapping with the larger number of matches.
        if len(m1) > 0:
            if len(m0) > 0:
                if len(m1[0]) > len(m0[0]):
                    mappings = m1
                else:
                    mappings = m0
            else:
                mappings = m1
        else:
            mappings = m0

        # Score the mappings and return them in sorted order (best to worst).
        mappings, scores = _score_sire_mappings(mol0, mol1, mappings, prematch,
            _scoring_function, property_map0, property_map1)

    # Optionally post-process the MCS.
    if prune_perturbed_constraints:
        mappings = [_prune_perturbed_constraints(molecule0, molecule1, x) for x in mappings]
    if prune_crossing_constraints:
        mappings = [_prune_crossing_constraints(molecule0, molecule1, x) for x in mappings]

    if matches == 1:
        if return_scores:
            return (mappings[0], scores[0])
        else:
            return mappings[0]
    else:
        # Return a list of matches from best to worst.
        if return_scores:
            return (mappings[0:matches], scores[0:matches])
        # Return a tuple containing the list of matches from best to
        # worst along with the list of scores.
        else:
            return mappings[0:matches]

def rmsdAlign(molecule0, molecule1, mapping=None, property_map0={}, property_map1={}):
    """Align atoms in molecule0 to those in molecule1 using the mapping
       between matched atom indices. The molecule is aligned using rigid-body
       translation and rotations, with a root mean squared displacement (RMSD)
       fit to find the optimal translation vector (as opposed to merely taking
       the difference of centroids).

       Parameters
       ----------

       molecule0 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to align.

       molecule1 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The reference molecule.

       mapping : dict
           A dictionary mapping atoms in molecule0 to those in molecule1.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The aligned molecule.

       Examples
       --------

       Align molecule0 to molecule1 based on a precomputed mapping.

       >>> import BioSimSpace as BSS
       >>> molecule0 = BSS.Align.rmsdAlign(molecule0, molecule1, mapping)

       Align molecule0 to molecule1. Since no mapping is passed one will be
       autogenerated using :class:`matchAtoms <BioSimSpace.Align.matchAtoms>`
       with default options.

       >>> import BioSimSpace as BSS
       >>> molecule0 = BSS.Align.rmsdAlign(molecule0, molecule1)
    """

    if not isinstance(molecule0, _Molecule):
        raise TypeError("'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(molecule1, _Molecule):
        raise TypeError("'molecule1' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(property_map0, dict):
        raise TypeError("'property_map0' must be of type 'dict'")

    if not isinstance(property_map1, dict):
        raise TypeError("'property_map1' must be of type 'dict'")

    # The user has passed an atom mapping.
    if mapping is not None:
        if not isinstance(mapping, dict):
            raise TypeError("'mapping' must be of type 'dict'.")
        else:
            _validate_mapping(molecule0, molecule1, mapping, "mapping")

    # Get the best match atom mapping.
    else:
        mapping = matchAtoms(molecule0, molecule1, property_map0=property_map0,
                             property_map1=property_map1)

    # Extract the Sire molecule from each BioSimSpace molecule.
    mol0 = molecule0._getSireObject()
    mol1 = molecule1._getSireObject()

    # Convert the mapping to AtomIdx key:value pairs.
    sire_mapping = _to_sire_mapping(mapping)

    # Perform the alignment, mol0 to mol1.
    try:
        mol0 = mol0.move().align(mol1, _SireMol.AtomResultMatcher(sire_mapping)).molecule()
    except Exception as e:
        msg = "Failed to align molecules based on mapping: %r" % mapping
        if "Could not calculate the single value decomposition" in str(e):
            msg += ". Try minimising your molecular coordinates prior to alignment."
        if _isVerbose():
            raise _AlignmentError(msg) from e
        else:
            raise _AlignmentError(msg) from None

    # Return the aligned molecule.
    return _Molecule(mol0)

def flexAlign(molecule0, molecule1, mapping=None, fkcombu_exe=None,
        property_map0={}, property_map1={}):
    """Flexibly align atoms in molecule0 to those in molecule1 using the
       mapping between matched atom indices.

       Parameters
       ----------

       molecule0 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule to align.

       molecule1 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The reference molecule.

       mapping : dict
           A dictionary mapping atoms in molecule0 to those in molecule1.

       fkcombu_exe : str
           Path to the fkcombu executable. If None is passed, then BioSimSpace
           will attempt to find fkcombu by searching your PATH.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The aligned molecule.

       Examples
       --------

       Align molecule0 to molecule1 based on a precomputed mapping.

       >>> import BioSimSpace as BSS
       >>> molecule0 = BSS.Align.flexAlign(molecule0, molecule1, mapping)

       Align molecule0 to molecule1. Since no mapping is passed one will be
       autogenerated using :class:`matchAtoms <BioSimSpace.Align.matchAtoms>`
       with default options.

       >>> import BioSimSpace as BSS
       >>> molecule0 = BSS.Align.flexAlign(molecule0, molecule1)
    """

    # Check that we found fkcombu in the PATH.
    if fkcombu_exe is None:
        if _fkcombu_exe is None:
            raise _MissingSoftwareError("'BioSimSpace.Align.flexAlign' requires the 'fkcombu' program: "
                                        "https://pdbj.org/kcombu")
        else:
            fkcombu_exe = _fkcombu_exe
    # Check that the user supplied executable exists.
    else:
        if not _os.path.isfile(fkcombu_exe):
            raise IOError("'fkcombu' executable doesn't exist: '%s'" % fkcombu_exe)

    if not isinstance(molecule0, _Molecule):
        raise TypeError("'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(molecule1, _Molecule):
        raise TypeError("'molecule1' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(property_map0, dict):
        raise TypeError("'property_map0' must be of type 'dict'")

    if not isinstance(property_map1, dict):
        raise TypeError("'property_map1' must be of type 'dict'")

    # The user has passed an atom mapping.
    if mapping is not None:
        if not isinstance(mapping, dict):
            raise TypeError("'mapping' must be of type 'dict'.")
        else:
            _validate_mapping(molecule0, molecule1, mapping, "mapping")

    # Get the best match atom mapping.
    else:
        mapping = matchAtoms(molecule0, molecule1, property_map0=property_map0,
                             property_map1=property_map1)

    # Convert the mapping to AtomIdx key:value pairs.
    sire_mapping = _to_sire_mapping(mapping)

    # Create a temporary working directory.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Execute in the working directory.
    with _Utils.cd(work_dir):

        # Write the two molecules to PDB files.
        _IO.saveMolecules("molecule0", molecule0, "PDB", property_map=property_map0)
        _IO.saveMolecules("molecule1", molecule1, "PDB", property_map=property_map1)

        # Write the mapping to text. (Increment indices by one).
        with open("mapping.txt", "w") as file:
            for idx0, idx1 in sire_mapping.items():
                file.write("%d %d\n" % (idx0.value() + 1, idx1.value() + 1))

        # Create the fkcombu command string.
        command = "%s -T molecule0.pdb -R molecule1.pdb -alg F -iam mapping.txt -opdbT aligned.pdb" % fkcombu_exe

        # Run the command as a subprocess.
        proc = _subprocess.run(_shlex.split(command), shell=False,
            stdout=_subprocess.PIPE, stderr=_subprocess.PIPE)

        # Check that the output file exists.
        if not _os.path.isfile("aligned.pdb"):
            raise _AlignmentError("Failed to align molecules based on mapping: %r" % mapping) from None

        # Load the aligned molecule.
        aligned = _IO.readMolecules("aligned.pdb")[0]

        # Get the "coordinates" property for molecule0.
        prop = property_map0.get("coordinates", "coordinates")

        # Copy the coordinates back into the original molecule.
        molecule0._sire_object = molecule0._sire_object.edit() \
            .setProperty(prop, aligned._sire_object.property("coordinates")).commit()

    # Return the aligned molecule.
    return _Molecule(molecule0)

def merge(molecule0, molecule1, mapping=None, allow_ring_breaking=False,
        allow_ring_size_change=False, force=False,
        property_map0={}, property_map1={}):
    """Create a merged molecule from 'molecule0' and 'molecule1' based on the
       atom index 'mapping'. The merged molecule can be used in single topology
       free-energy simulations.

       Parameters
       ----------

       molecule0 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           A molecule object.

       molecule1 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           A second molecule object.

       mapping : dict
           The mapping between matching atom indices in the two molecules.
           If no mapping is provided, then atoms in molecule0 will be mapped
           to those in molecule1 using "matchAtoms", with "rmsdAlign" then
           used to align molecule0 to molecule1 based on the resulting mapping.

       allow_ring_breaking : bool
           Whether to allow the opening/closing of rings during a merge.

       allow_ring_size_change : bool
           Whether to allow changes in ring size.

       force : bool
           Whether to try to force the merge, even when the molecular
           connectivity changes not as the result of a ring transformation.
           This will likely lead to an unstable perturbation. This option
           takes precedence over 'allow_ring_breaking' and
           'allow_ring_size_change'.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       Returns
       -------

       molecule : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The merged molecule.

       Examples
       --------

       Merge molecule0 and molecule1 based on a precomputed mapping.

       >>> import BioSimSpace as BSS
       >>> merged = BSS.Align.merge(molecule0, molecule1, mapping)

       Merge molecule0 with molecule1. Since no mapping is passed one will be
       autogenerated using :class:`matchAtoms <BioSimSpace.Align.matchAtoms>`
       with default options, following which :class:`rmsdAlign <BioSimSpace.Align.rmsdAlign>`
       will be used to align molecule0 to molecule1 based on the resulting mapping.

       >>> import BioSimSpace as BSS
       >>> molecule0 = BSS.Align.merge(molecule0, molecule1)
    """

    if not isinstance(molecule0, _Molecule):
        raise TypeError("'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(molecule1, _Molecule):
        raise TypeError("'molecule1' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(property_map0, dict):
        raise TypeError("'property_map0' must be of type 'dict'")

    if not isinstance(property_map1, dict):
        raise TypeError("'property_map1' must be of type 'dict'")

    if not isinstance(allow_ring_breaking, bool):
        raise TypeError("'allow_ring_breaking' must be of type 'bool'")

    if not isinstance(allow_ring_size_change, bool):
        raise TypeError("'allow_ring_size_change' must be of type 'bool'")

    if not isinstance(force, bool):
        raise TypeError("'force' must be of type 'bool'")

    # The user has passed an atom mapping.
    if mapping is not None:
        if not isinstance(mapping, dict):
            raise TypeError("'mapping' must be of type 'dict'.")
        else:
            _validate_mapping(molecule0, molecule1, mapping, "mapping")

    # Get the best atom mapping and align molecule0 to molecule1 based on the
    # mapping.
    else:
        mapping = matchAtoms(molecule0, molecule1,
            property_map0=property_map0, property_map1=property_map1)
        molecule0 = rmsdAlign(molecule0, molecule1, mapping)

    # Convert the mapping to AtomIdx key:value pairs.
    sire_mapping = _to_sire_mapping(mapping)

    # Create and return the merged molecule.
    return _merge(molecule0, molecule1, sire_mapping, allow_ring_breaking=allow_ring_breaking,
            allow_ring_size_change=allow_ring_size_change, force=force,
            property_map0=property_map0, property_map1=property_map1)

def viewMapping(molecule0, molecule1, mapping=None, property_map0={},
        property_map1={}, style=None, orientation="horizontal", pixels=900):
    """Visualise the mapping between molecule0 and molecule1. This draws a 3D
       depiction of both molecules with the mapped atoms highlighted in green.
       Labels specify the indices of the atoms, along with the indices of
       the atoms to which they map in the other molecule.

       Parameters
       ----------

       molecule0 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The first molecule.

       molecule1 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The second molecule.

       mapping : dict
           A dictionary mapping atoms in molecule0 to those in molecule1.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       style : dict
           Drawing style. See https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html
           for some examples.

       orientation : str
           Whether to display the two molecules in a "horizontal" or "vertical"
           arrangement.

       pixels : int
           The size of the largest view dimension in pixel, i.e. either the
           "horizontal" or "vertical" size.

       Returns
       -------

       view : py3Dmol.view
           A view of the two molecules with the mapped atoms highlighted and
           labelled.
    """

    # Adapted from: https://gist.github.com/cisert/d05664d4c98ac1cf86ee70b8700e56a9

    # Only draw within a notebook.
    if not _is_notebook:
        return None

    _assert_imported(_rdkit)

    if not isinstance(molecule0, _Molecule):
        raise TypeError("'molecule0' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(molecule1, _Molecule):
        raise TypeError("'molecule1' must be of type 'BioSimSpace._SireWrappers.Molecule'")

    if not isinstance(property_map0, dict):
        raise TypeError("'property_map0' must be of type 'dict'")

    if not isinstance(property_map1, dict):
        raise TypeError("'property_map1' must be of type 'dict'")

    if style is not None:
        if not isinstance(style, dict):
            raise TypeError("'style' must be of type 'dict'")

    if not isinstance(orientation, str):
        raise TypeError("'orientation' must be of type 'str'")
    else:
        # Convert to lower case and strip whitespace.
        orientation = orientation.lower().replace(" ", "")

        if orientation not in ["horizontal", "vertical"]:
            raise ValueError("'orientation' must be equal to 'horizontal' "
                             "or 'vertical'.")

    if isinstance(pixels, float):
        pixels = int(pixels)
    if not type(pixels) is int:
        raise TypeError("'pixels' must be of type 'int'")
    if pixels <= 0:
        raise ValueError("pixels' must be > 0!")

    # The user has passed an atom mapping.
    if mapping is not None:
        if not isinstance(mapping, dict):
            raise TypeError("'mapping' must be of type 'dict'.")
        else:
            _validate_mapping(molecule0, molecule1, mapping, "mapping")

    # Get the best atom mapping and align molecule0 to molecule1 based on the
    # mapping.
    else:
        mapping = matchAtoms(molecule0, molecule1,
            property_map0=property_map0, property_map1=property_map1)
        molecule0 = rmsdAlign(molecule0, molecule1, mapping)

    # Create a temporary working directory and store the directory name.
    tmp_dir = _tempfile.TemporaryDirectory()
    work_dir = tmp_dir.name

    # Write the molecules to PDB format.
    _IO.saveMolecules(work_dir + "/molecule0",
        molecule0, "pdb", property_map=property_map0)
    _IO.saveMolecules(work_dir + "/molecule1",
        molecule1, "pdb", property_map=property_map0)

    import py3Dmol as _py3Dmol

    # Load the molecules into RDKit.
    rdmol0 = _Chem.MolFromPDBFile(work_dir + "/molecule0.pdb",
        sanitize=False, removeHs=False)
    rdmol1 = _Chem.MolFromPDBFile(work_dir + "/molecule1.pdb",
        sanitize=False, removeHs=False)

    # Set grid view properties.
    viewer0 = (0, 0)
    if orientation == "horizontal":
        viewergrid = (1, 2)
        viewer1 = (0, 1)
        width = pixels
        height = pixels / 2
    else:
        viewergrid = (2, 1)
        viewer1 = (1, 0)
        width = pixels / 2
        height = pixels

    # Create the view.
    view = _py3Dmol.view(linked=False, width=width,
            height=height, viewergrid=viewergrid)

    # Set default drawing style.
    if style is None:
        style = {"stick":{"colorscheme":"grayCarbon", "linewidth": 0.1}}

    # Add the molecules to the views.
    view.addModel(_Chem.MolToMolBlock(rdmol0), "mol0", viewer=viewer0)
    view.addModel(_Chem.MolToMolBlock(rdmol1), "mol1", viewer=viewer1)

    # Set the style.
    view.setStyle({"model": 0}, style, viewer=viewer0)
    view.setStyle({"model": 0}, style, viewer=viewer1)

    # Highlight the atoms from the mapping.
    for atom0, atom1 in mapping.items():
        p = rdmol0.GetConformer().GetAtomPosition(atom0)
        view.addSphere({"center" : {"x" : p.x, "y" : p.y, "z" : p.z},
                        "radius" : 0.5,
                        "color"  : "green", "alpha": 0.8},
                       viewer=viewer0)
        view.addLabel(f"{atom0} \u2192 {atom1}",
                      {"position" : {"x" : p.x, "y" : p.y, "z" : p.z}},
                      viewer=viewer0)
        p = rdmol1.GetConformer().GetAtomPosition(atom1)
        view.addSphere({"center" : {"x" : p.x, "y" : p.y, "z" : p.z},
                        "radius" : 0.5,
                        "color"  : "green", "alpha": 0.8},
                       viewer=viewer1)
        view.addLabel(f"{atom1} \u2192 {atom0}",
                      {"position" : {"x" : p.x, "y" : p.y, "z" : p.z}},
                      viewer=viewer1)

    # Set background colour.
    view.setBackgroundColor("white", viewer=viewer0)
    view.setBackgroundColor("white", viewer=viewer1)

    # Zoom to molecule.
    view.zoomTo(viewer=viewer0)
    view.zoomTo(viewer=viewer1)

    return view

def _score_rdkit_mappings(molecule0, molecule1, rdkit_molecule0, rdkit_molecule1,
        mcs_smarts, prematch, scoring_function, max_scoring_matches, property_map0,
        property_map1):
    """Internal function to score atom mappings based on the root mean squared
       displacement (RMSD) between mapped atoms in two molecules. Optionally,
       molecule0 can first be aligned to molecule1 based on the mapping prior
       to computing the RMSD. The function returns the mappings sorted based
       on their score from best to worst, along with a list containing the
       scores for each mapping.

       Parameters
       ----------

       molecule0 : Sire.Molecule.Molecule
           The first molecule (Sire representation).

       molecule0 : Sire.Molecule.Molecule
           The second molecule (Sire representation).

       rdkit_mol0 : RDKit.Chem.Mol
           The first molecule (RDKit representation).

       rdkit_mol1 : RDKit.Chem.Mol
           The second molecule (RDKit representation).

       mcs_smarts : RDKit.Chem.MolFromSmarts
           The smarts string representing the maximum common substructure of
           the two molecules.

       prematch : dict
           A dictionary of atom mappings that must be included in the match.

       scoring_function : str
           The RMSD scoring function.

       max_scoring_matches : int
           The maximum number of matching MCS substructures to consider when
           computing mapping scores. Consider reducing this if you find the
           matchAtoms function to be taking an excessive amount of time.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       Returns
       -------

       mapping, scores : ([dict], list)
           The ranked mappings and corresponding scores.
    """

    # Adapted from FESetup: https://github.com/CCPBioSim/fesetup

    # Make sure to re-map the coordinates property in both molecules, otherwise
    # the move and align functions from Sire will not work.
    prop0 = property_map0.get("coordinates", "coordinates")
    prop1 = property_map1.get("coordinates", "coordinates")

    if prop0 != "coordinates":
        molecule0 = molecule0.edit().setProperty("coordinates", molecule0.property(prop0)).commit()
    if prop1 != "coordinates":
        molecule1 = molecule1.edit().setProperty("coordinates", molecule1.property(prop1)).commit()

    # Get the set of matching substructures in each molecule. For some reason
    # setting uniquify to True removes valid matches, in some cases even the
    # best match! As such, we set uniquify to False and account for duplicate
    # mappings in the code below.

    matches0 = rdkit_molecule0.GetSubstructMatches(mcs_smarts, uniquify=False,
        maxMatches=max_scoring_matches, useChirality=False)
    matches1 = rdkit_molecule1.GetSubstructMatches(mcs_smarts, uniquify=False,
        maxMatches=max_scoring_matches, useChirality=False)

    # Swap the order of the matches.
    if len(matches0) < len(matches1):
        matches0, matches1 = matches1, matches0
        is_swapped = True
    else:
        is_swapped = False

    # Initialise a list to hold the mappings.
    mappings = []

    # Initialise a list of to hold the score for each mapping.
    scores = []

    # Whether there was a GSL alignment error.
    is_gsl_error = False

    # Loop over all matches from mol0.
    for x in range(len(matches0)):
        match0 = matches0[x]

        # Loop over all matches from mol1.
        for y in range(len(matches1)):
            match1 = matches1[y]

            # Initialise the mapping for this match.
            mapping = {}
            sire_mapping = {}

            # Loop over all atoms in the match.
            for i, idx0 in enumerate(match0):
                idx1 = match1[i]

                # Add to the mapping.
                if is_swapped:
                    mapping[idx1] = idx0
                    sire_mapping[_SireMol.AtomIdx(idx1)] = _SireMol.AtomIdx(idx0)
                else:
                    mapping[idx0] = idx1
                    sire_mapping[_SireMol.AtomIdx(idx0)] = _SireMol.AtomIdx(idx1)

            # This is a new mapping:
            if not mapping in mappings:
                # Check that the mapping contains the pre-match.
                is_valid = True
                for idx0, idx1 in prematch.items():
                    # Pre-match isn't found, return to top of loop.
                    if idx0 not in mapping or mapping[idx0] != idx1:
                        is_valid = False
                        break

                if is_valid:
                    # Rigidly align molecule0 to molecule1 based on the mapping.
                    if scoring_function == "RMSDALIGN":
                        try:
                            molecule0 = molecule0.move().align(molecule1, _SireMol.AtomResultMatcher(sire_mapping)).molecule()
                        except Exception as e:
                            if "Could not calculate the single value decomposition" in str(e):
                                is_gsl_error = True
                                gsl_exception = e
                            else:
                                msg = "Failed to align molecules when scoring based on mapping: %r" % mapping
                                if _isVerbose():
                                    raise _AlignmentError(msg) from e
                                else:
                                    raise _AlignmentError(msg) from None
                    # Flexibly align molecule0 to molecule1 based on the mapping.
                    elif scoring_function == "RMSDFLEXALIGN":
                        molecule0 = flexAlign(_Molecule(molecule0), _Molecule(molecule1), mapping,
                            property_map0=property_map0, property_map1=property_map1)._sire_object

                    # Append the mapping to the list.
                    mappings.append(mapping)

                    # We now compute the RMSD between the coordinates of the matched atoms
                    # in molecule0 and molecule1.

                    # Initialise lists to hold the coordinates.
                    c0 = []
                    c1 = []

                    # Loop over each atom index in the map.
                    for idx0, idx1 in sire_mapping.items():
                        # Append the coordinates of the matched atom in molecule0.
                        c0.append(molecule0.atom(idx0).property("coordinates"))
                        # Append the coordinates of atom in molecule1 to which it maps.
                        c1.append(molecule1.atom(idx1).property("coordinates"))

                    # Compute the RMSD between the two sets of coordinates.
                    scores.append(_SireMaths.getRMSD(c0, c1))

    # No mappings were found.
    if len(mappings) == 0:
        # We failed to align mappings during scoring due to convergence issues
        # during the GSL single value decomposition routine.
        if is_gsl_error:
            msg = "Failed to align molecules when scoring. " \
                  "Try minimising your molecular coordinates prior calling matchAtoms."
            if _isVerbose():
                raise _AlignmentError(msg) from gsl_exception
            else:
                raise _AlignmentError(msg) from None

        if len(prematch) == 0:
            return ([{}], [])
        else:
            return ([prematch], [])

    # Sort the scores and return the sorted keys. (Smaller RMSD is best)
    keys = sorted(range(len(scores)), key=lambda k: scores[k])

    # Sort the mappings.
    mappings = [mappings[x] for x in keys]

    # Sort the scores and convert to Angstroms.
    scores = [scores[x] * _Units.Length.angstrom for x in keys]

    # Return the sorted mappings and their scores.
    return (mappings, scores)

def _score_sire_mappings(molecule0, molecule1, sire_mappings, prematch,
        scoring_function, property_map0, property_map1):
    """Internal function to score atom mappings based on the root mean squared
       displacement (RMSD) between mapped atoms in two molecules. Optionally,
       molecule0 can first be aligned to molecule1 based on the mapping prior
       to computing the RMSD. The function returns the mappings sorted based
       on their score from best to worst, along with a list containing the
       scores for each mapping.

       Parameters
       ----------

       molecule0 : Sire.Molecule.Molecule
           The first molecule (Sire representation).

       molecule0 : Sire.Molecule.Molecule
           The second molecule (Sire representation).

       sire_mappings : [{}]
           The list of mappings generated by Sire.

       prematch : dict
           A dictionary of atom mappings that must be included in the match.

       scoring_function : str
           The RMSD scoring function.

       property_map0 : dict
           A dictionary that maps "properties" in molecule0 to their user
           defined values. This allows the user to refer to properties
           with their own naming scheme, e.g. { "charge" : "my-charge" }

       property_map1 : dict
           A dictionary that maps "properties" in molecule1 to their user
           defined values.

       Returns
       -------

       mapping, scores : ([dict], list)
           The ranked mappings and corresponding scores.
    """

    # Make sure to re-map the coordinates property in both molecules, otherwise
    # the move and align functions from Sire will not work.
    prop0 = property_map0.get("coordinates", "coordinates")
    prop1 = property_map1.get("coordinates", "coordinates")

    if prop0 != "coordinates":
        molecule0 = molecule0.edit().setProperty("coordinates", molecule0.property(prop0)).commit()
    if prop1 != "coordinates":
        molecule1 = molecule1.edit().setProperty("coordinates", molecule1.property(prop1)).commit()

    # Initialise a list to hold the mappings.
    mappings = []

    # Initialise a list of to hold the score for each mapping.
    scores = []

    # Loop over all of the mappings.
    for mapping in sire_mappings:

        # Check that the mapping contains the pre-match.
        is_valid = True
        for idx0, idx1 in prematch.items():
            # Pre-match isn't found, return to top of loop.
            if _SireMol.AtomIdx(idx0) not in mapping or mapping[_SireMol.AtomIdx(idx0)] != _SireMol.AtomIdx(idx1):
                is_valid = False
                break

        if is_valid:
            # Rigidly align molecule0 to molecule1 based on the mapping.
            if scoring_function == "RMSDALIGN":
                try:
                    molecule0 = molecule0.move().align(molecule1, _SireMol.AtomResultMatcher(mapping)).molecule()
                except Exception as e:
                    msg = "Failed to align molecules when scoring based on mapping: %r" % mapping
                    if _isVerbose():
                        raise _AlignmentError(msg) from e
                    else:
                        raise _AlignmentError(msg) from None
            # Flexibly align molecule0 to molecule1 based on the mapping.
            elif scoring_function == "RMSDFLEXALIGN":
                molecule0 = flexAlign(_Molecule(molecule0), _Molecule(molecule1), _from_sire_mapping(mapping),
                    property_map0=property_map0, property_map1=property_map1)._sire_object

            # Append the mapping to the list.
            mappings.append(_from_sire_mapping(mapping))

            # We now compute the RMSD between the coordinates of the matched atoms
            # in molecule0 and molecule1.

            # Initialise lists to hold the coordinates.
            c0 = []
            c1 = []

            # Loop over each atom index in the map.
            for idx0, idx1 in mapping.items():
                # Append the coordinates of the matched atom in molecule0.
                c0.append(molecule0.atom(idx0).property("coordinates"))
                # Append the coordinates of atom in molecule1 to which it maps.
                c1.append(molecule1.atom(idx1).property("coordinates"))

            # Compute the RMSD between the two sets of coordinates.
            scores.append(_SireMaths.getRMSD(c0, c1))

    # No mappings were found.
    if len(mappings) == 0:
        if len(prematch) == 0:
            return ([{}], [])
        else:
            return ([prematch], [])

    # Sort the scores and return the sorted keys. (Smaller RMSD is best)
    keys = sorted(range(len(scores)), key=lambda k: scores[k])

    # Sort the mappings.
    mappings = [mappings[x] for x in keys]

    # Sort the scores and convert to Angstroms.
    scores = [scores[x] * _Units.Length.angstrom for x in keys]

    # Return the sorted mappings and their scores.
    return (mappings, scores)

def _validate_mapping(molecule0, molecule1, mapping, name):
    """Internal function to validate that a mapping contains key:value pairs
       of the correct type.

       Parameters
       ----------

       molecule0 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The molecule of interest.

       molecule1 : :class:`Molecule <BioSimSpace._SireWrappers.Molecule>`
           The reference molecule.

       mapping : dict
           The mapping between matching atom indices in the two molecules.

       name : str
           The name of the mapping. (Used when raising exceptions.)
    """

    for idx0, idx1 in mapping.items():
            if type(idx0) is int and type(idx1) is int:
                pass
            elif isinstance(idx0, _SireMol.AtomIdx) and isinstance(idx1, _SireMol.AtomIdx):
                idx0 = idx0.value()
                idx1 = idx1.value()
            else:
                raise TypeError("%r dictionary key:value pairs must be of type 'int' or "
                                "'Sire.Mol.AtomIdx'" % name)
            if idx0 < 0 or idx0 >= molecule0.nAtoms() or \
               idx1 < 0 or idx1 >= molecule1.nAtoms():
                raise ValueError("%r dictionary key:value pair '%s : %s' is out of range! "
                                 "The molecules contain %d and %d atoms."
                                 % (name, idx0, idx1, molecule0.nAtoms(), molecule1.nAtoms()))

def _to_sire_mapping(mapping):
    """Internal function to convert a regular mapping to Sire AtomIdx format.

       Parameters
       ----------

       mapping : {int:int}
           The regular mapping.

       Returns
       -------

       sire_mapping : {Sire.Mol.AtomIdx:Sire.Mol.AtomIdx}
           The Sire mapping.
    """

    sire_mapping = {}

    # Convert the mapping to AtomIdx key:value pairs.
    for idx0, idx1 in mapping.items():
        # Early exit if the mapping is already the correct format.
        if isinstance(idx0, _SireMol.AtomIdx):
            return mapping
        else:
            sire_mapping[_SireMol.AtomIdx(idx0)] = _SireMol.AtomIdx(idx1)

    return sire_mapping

def _from_sire_mapping(sire_mapping):
    """Internal function to convert from a Sire mapping to regular format.

       Parameters
       ----------

       sire_mapping : {Sire.Mol.AtomIdx:Sire.Mol.AtomIdx}
           The Sire mapping.

       Returns
       -------

       mapping : {int:int}
           The regular mapping.
    """

    mapping = {}

    # Convert the mapping to int key:value pairs.
    for idx0, idx1 in sire_mapping.items():
        # Early exit if the mapping is already the correct format.
        if type(idx0) is int:
            return sire_mapping
        else:
            mapping[idx0.value()] = idx1.value()

    return mapping

def _prune_perturbed_constraints(molecule0, molecule1, mapping):
    """Prunes the maximum common substructure (MCS) so that no hydrogen
       bond constraints are perturbed.

       Parameters
       ----------
       molecule0 : BioSimSpace._SireWrappers.Molecule
           The first molecule (used at lambda = 0).
       molecule1 : BioSimSpace._SireWrappers.Molecule
           The second molecule (used at lambda = 1).
       mapping : dict(int, int)
           A maximum common substructure mapping between both molecules, as generated by e.g.
           BioSimSpace.Align.matchAtoms().

       Returns
       -------
       new_mapping : dict(int, int)
           The pruned MCS.
    """
    new_mapping = {}

    for idx0, idx1 in mapping.items():
        atom0 = molecule0.getAtoms()[idx0]
        atom1 = molecule1.getAtoms()[idx1]
        elem0 = atom0._sire_object.property("element").symbol()
        elem1 = atom1._sire_object.property("element").symbol()
        elems = {elem0, elem1}

        # Make sure we are not matching a hydrogen to a non-hydrogen
        if not ("H" in elems and len(elems) > 1):
            new_mapping[idx0] = idx1

    return new_mapping

def _prune_crossing_constraints(molecule0, molecule1, mapping):
    """Prunes the maximum common substructure (MCS) mapping so that there are no
       constrained bonds between a common core and a softcore atom.
       
       Parameters
       ----------
       molecule0 : BioSimSpace._SireWrappers.Molecule
           The first molecule (used at lambda = 0).
       molecule1 : BioSimSpace._SireWrappers.Molecule
           The second molecule (used at lambda = 1).
       mapping : dict(int, int)
           A maximum common substructure mapping between both molecules, as generated by e.g.
           BioSimSpace.Align.matchAtoms().

       Returns
       -------
       new_mapping : dict(int, int)
           The pruned mapping.
    """
    connectivity0 = _SireMol.Connectivity(
        molecule0._sire_object, _SireMol.CovalentBondHunter()
    )
    connectivity1 = _SireMol.Connectivity(
        molecule1._sire_object, _SireMol.CovalentBondHunter()
    )

    while True:
        new_mapping = {}

        for idx0, idx1 in mapping.items():
            # Get the relevant atom and whether it's a hydrogen
            atom0 = molecule0._sire_object.atom(_SireMol.AtomIdx(idx0))
            atom1 = molecule1._sire_object.atom(_SireMol.AtomIdx(idx1))
            is_H0 = "H" in atom0.property("element").symbol()
            is_H1 = "H" in atom1.property("element").symbol()

            # Get the neighbours to the atom
            neighbours0 = [
                molecule0._sire_object.atom(i)
                for i in connectivity0.connectionsTo(_SireMol.AtomIdx(idx0))
            ]
            neighbours1 = [
                molecule1._sire_object.atom(i)
                for i in connectivity1.connectionsTo(_SireMol.AtomIdx(idx1))
            ]

            # Determine whether there are any constrained bonds between the MCS and softcore part
            any_Hdummies0 = any(
                ("H" in atom.property("element").symbol() or is_H0)
                and atom.index().value() not in mapping.keys()
                for atom in neighbours0
            )
            any_Hdummies1 = any(
                ("H" in atom.property("element").symbol() or is_H1)
                and atom.index().value() not in mapping.values()
                for atom in neighbours1
            )

            if not any_Hdummies0 and not any_Hdummies1:
                new_mapping[idx0] = idx1

        # We stop iterating if the pruned mapping is the same as the input one
        if new_mapping == mapping:
            mapping = new_mapping
            break
        else:
            mapping = new_mapping

    return mapping
