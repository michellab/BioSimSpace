def decouple(molecule, property_map0={}, property_map1={}, intramol=True):
    """Make the molecule as being decoupled, where the interactions with the
    rest of the environment are removed, or annihilate, where the interactions
    within the molecule are removed as well (choose this mode with
    intramol=False).

        Parameters
        ----------

        molecule0 : BioSimSpace._SireWrappers.Molecule
            The molecule to be decoupled or annihilated.

        property_map0 : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values at the start of the transformation. This allows
            user to selectively turn on or off the charge or the vdw
            interactions e.g. { "charge" : True, "vdw" : True}

        property_map1 : dict
            A dictionary that maps "properties" in this molecule to their
            user defined values at the end of the transformation. This allows
            user to selectively turn on or off the charge or the vdw
            interactions e.g. { "charge" : False, "vdw" : False}

        intramol : bool
            Whether to remove the intra-molecule forces, thereby transforming


        Returns
        -------

        merged : Sire.Mol.Molecule
            The merged molecule.
    """
