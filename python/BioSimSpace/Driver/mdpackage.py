"""
@package biosimspace
@author  Lester Hedges
@brief   Functionality for driving molecular dynamics simulations.
"""

from Sire.Base import findExe

import BioSimSpace.Process as Process
import BioSimSpace.Protocol as Protocol

# A dictionary mapping MD packages to their executable names and GPU support.
#                PACKAGE        EXE       GPU
_md_packages = { "AMBER"   : { "sander" : False, "pmemd" : False, "pmemd.cuda" : True },
                 "GROMACS" : { "gmx" : True },
                 "NAMD"    : { "namd2" : False }
               }

def findMDPackage(supports_gpu=False):
    """Find a molecular dynamics package on the system and return
       a handle to it as a MDPackage object.

       Keyword arguments:

       supports_gpu -- Whether GPU support is required.
    """

    if type(supports_gpu) is not bool:
        raise TypeError("'supports_gpu' keyword must be of type 'bool'")

    return MDPackage()

class MDPackage():
    """A class for driving molecular dynamics simulations."""

    def __init__(self):
        """Constructor."""
