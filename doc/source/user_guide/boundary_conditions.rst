.. _boundary_conditions:

Boundary conditions
===================

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
within the ``.udf`` file when used (E.G setting a ``codedFixedValue`` boundary will 
require the ``velocityDirichletConditions`` to be present in the OKL section).

.. tip::

   To setup some cases (E.G. periodic boundary conditions), you may need to use the 
   ``boundaryIDMap`` parameter of the ``MESH`` section to apply the ``boundaryTypeMap``
   options to the correct boundary IDs (NekRS assumed that the boundaryIDs start at 1).

Zero
----

codedFixedValue
---------------

Interpolation
-------------

zero<X/Y/Z/N>Value / zeroGradient
---------------------------------

zero<X/Y/Z/N>Value VELOCITY only.

zero<X/Y/Z/N>Value / codedFixedGradient
---------------------------------------

VELOCITY only

zero<X/Y/Z/N>Value / fixedGradient
----------------------------------

codedFixedGradient
------------------

SCALAR only

None / Periodic
---------------

None is typically used when a periodic boundary condition has been set as part
of the mesh and it does not need to be considered as part of the standard processing
of boundary conditions.

.. tip:: 

   Each of the different meshing methods from :ref:`meshing` have a different
   way in which the periodic boundary will be defined.

   .. tabs::

      .. tab:: Scripting

      .. tab:: Cubit

      .. tab:: Gmsh
      
      .. tab:: CGNS

.. NekRS supports periodic boundary conditions. To set up a periodic case, first
.. you need to run ``exo2nek`` to establish the pairings between the periodic sidesets.
.. All this information will be prompted on the screen by ``exo2nek``;
.. You will provide the sideset IDs of the periodic boundaries, a search tolerance
.. for identifying paired sides, and a translation vector that points from one of the
.. paired sidesets to the other. For example, if you want to have one periodic surface
.. that is a :math:`z`-plane at :math:`z=-1.0` that is paired to another :math:`z`-plane
.. at :math:`z=1.0`, the translation vector would be :math:`(0.0, 0.0, 2.0)`.

.. _periodic_boundary:

Turbulent Inflow
----------------

Velocity Recycling
""""""""""""""""""

