""" The ADIOS2 high-level API module
License:
  Distributed under the OSI-approved Apache License, Version 2.0.  See
  accompanying file Copyright.txt for details.
"""

import adios2.bindings

from adios2.adios import *
from adios2.attribute import *
from adios2.engine import *
from adios2.io import *
from adios2.operator import *
from adios2.stream import *
from adios2.variable import *
from adios2.file_reader import *
from adios2.bindings import LocalValueDim

__license__ = "Apache-2.0"
__version__ = adios2.bindings.__version__
is_built_with_mpi = adios2.bindings.is_built_with_mpi
