# BioSimSpace.IO

This package contains functionality for reading and writing a large number of
common molecular file formats. Unlike many other molecular parsers, BioSimSpace
does not attempt to infer the format from the file extension. Instead, it reads
the files in parallel using its blazing fast suite of parsers in order to
automatically determine the correct format.

To read a system from file:

```python
import BioSimSpace as BSS

# Load a molecular system from a set of files.
system = BSS.IO.readMolecules(["ala.crd", "ala.top"])
```

It is possible to mix and match formats as long as the molecular information in the
different files is consistent. BioSimSpace will report an error whenever it is
unable to parse a particular file, or whenever it encounters conflicting information.

It is possible to query the system to see the file formats from which it was constructed:

```python
print(system.fileFormat())
'PRM7,RST7'
```

(Note that the file formats need not match the extensions of the original files.)

BioSimSpace also makes it easy to write a molecular system:

```python
# Save the system using the original file format, i.e. "filename.prm7, filename.rst7"
BSS.IO.saveMolecules("filename", system, system.fileFormat())
```

Where possible, BioSimSpace can also interconvert between formats. For example,
to save a [PDB](https://www.rcsb.org) representation of the system:

```python
BSS.IO.saveMolecules("filename", system, "PDB")
```

To see what file formats are available within BioSimSpace:

```python
BSS.IO.fileFormats()
```

To get a description of a given fileformat:

```python
BSS.IO.formatInfo("PDB")
'Protein Data Bank (PDB) format files.'
```
