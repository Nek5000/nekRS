.. _boundary_initial_conditions:

Boundary conditions
===============================

Boundary conditions for each mesh boundary should normally be set in the 
:ref:`parameter_file`, using the ``boundaryTypeMap`` parameter. This is used 
within the ``VELOCITY``, ``TEMPERATURE`` or ``SCALARXX`` sections to set the 
boundary conditions of the respective components of the case.

The potential values are summarised in the command line help function for 
nekRS, ``parHelp.txt`` and below.

.. literalinclude:: ../../parHelp.txt
   :language: none
   :lines: 1-3, 144-166

For each of the values that have a function name present, this must be created
within the udf file when used (E.G setting a ``codedFixedValue`` boundary will 
require the ``velocityDirichletConditions`` to be present).

To setup some cases (E.G. periodic boundary conditions or other 
scenarios where conditions have been applied to the mesh), you may need to use the 
``boundaryIDMap`` parameter of the ``MESH`` section to apply the ``boundaryTypeMap``
options to the correct boundary IDs. NekRS assumed that the boundaryIDs start at 1.

Types
-----

Periodic
""""""""

Periodic boundaries are not considered to be boundary conditional as part of the 
``boundaryTypeMap`` and must be set to ``none`` for the relevant boundaries. 
Instead it must be set as part of the mesh generation (see :ref:`periodic_mesh`)

Flow - constFlowRate, boundary force???????????

Plugins
-------

Velocity Recycling
""""""""""""""""""

