import BioSimSpace as BSS

import glob
import os
import pytest
import tempfile


def test_file_cache():
    """
    Simple test to see if cached files are used when repeatedly writing
    the same system to the same file format.
    """

    # Clear the file cache.
    BSS.IO._file_cache._cache = {}

    # Load the molecular system.
    s = BSS.IO.readMolecules(["test/input/ala.crd", "test/input/ala.top"])

    # Create a temporary working directory.
    tmp_dir = tempfile.TemporaryDirectory()
    tmp_path = tmp_dir.name

    # Write the system to PDB and PRM7 format.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp", s, ["pdb", "prm7"])

    # Check that the file cache has two entries.
    assert len(BSS.IO._file_cache._cache) == 2

    # Write to PDB and GroTop format. The PDB from the cache should be reused.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp2", s, ["pdb", "grotop"])

    # Check that the file cache has three entries, i.e. only the GroTop was added.
    assert len(BSS.IO._file_cache._cache) == 3

    # The directory should now contain 4 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 4

    # Now delete one of the files on disk and re-write. This should create a new
    # entry in the cache.
    os.remove(f"{tmp_path}/tmp.pdb")
    BSS.IO.saveMolecules(f"{tmp_path}/tmp3", s, ["pdb", "grotop"])

    # Check that the file cache still has three entries.
    assert len(BSS.IO._file_cache._cache) == 3

    # The directory should now contain 5 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 5

    # Now "corrupt" a file on disk so that its MD5 checksum is no longer
    # valid.
    with open(f"{tmp_path}/tmp2.pdb", "w") as f:
        pass

    # Write back to PDB and Gro87 format. The PDB file is now invalid, so
    # a new one will be written and added to the cache.
    BSS.IO.saveMolecules(f"{tmp_path}/tmp4", s, ["pdb", "gro87"])

    # Check that the file cache still has three entries.
    assert len(BSS.IO._file_cache._cache) == 4

    # The directory should now contain 7 files.
    assert len(glob.glob(f"{tmp_path}/*")) == 7
