.. _boundary_conditions:

Boundary conditions
===================

Boundary conditions for each mesh boundary should normally be set in the 
:ref:`parameter_file`, using the ``boundaryTypeMap`` parameter. This is used 
within the ``VELOCITY``, ``TEMPERATURE`` or ``SCALARXX`` sections to set the 
boundary conditions of the respective components of the case.

Available Types
---------------

The potential values are summarised in the command line help function for 
nekRS, ``parHelp.txt`` and below.

.. literalinclude:: ../../parHelp.txt
   :language: none
   :lines: 1-3, 149-179

.. tip::

   To setup some cases, you may need to use the ``boundaryIDMap`` parameter of 
   the ``MESH`` section to apply the ``boundaryTypeMap`` options to the correct
   boundary IDs of the mesh (By default NekRS assumes that the boundaryIDs start
   at 1). Below is an example where the four boundary conditions are
   applied to the boundary IDs 389 (``codeFixedValue``), 231 (``zeroGradient``),
   4 (``zeroValue``) and 23 (``zeroValue``).

   .. code-block::

      [VELOCITY]
      boundaryTypeMap = codedFixedValue, zeroGradient, zeroValue, zeroValue

      [MESH]
      boundaryIDMap = 389, 231, 4, 23

User Defined Value/Gradient
"""""""""""""""""""""""""""

If a boundary condition requires a value setting rather than just a type, a 
suitable function will need to be provided within the ``.udf`` file. The name of
the function used should be in the form of ``codedFixed`` + ``Value/Gradient`` 
+ ``Velocity/Scalar/Mesh``, E.G. ``codedFixedValueVelocity``, 
``codedFixedGradientScalar`` and ``codedFixedValueMesh``.

All of these functions are passed the ``bcData`` struct which has the following
parameters available:

+-----------------------------------------------------+------------------------------+------------------------------+
| Name                                                | Type                         | Description                  |
+=====================================================+==============================+==============================+
| ``idM``                                             | ``int``                      |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``fieldOffset``                                     | ``int``                      |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``id``                                              | ``int``                      |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``time``                                            | ``double``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``x``, ``y``, ``z``                                 | ``dfloat``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``nx``, ``ny``, ``nz``                              | ``dfloat``                   | normals                      |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``1x``, ``t1y``, ``t1z``, ``t2x``, ``t2y``, ``t2z`` | ``dfloat``                   | tangential directions        |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``tr1``, ``tr2``                                    | ``dfloat``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``u``, ``v``, ``w``                                 | ``dfloat``                   | velocity                     |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``p``                                               | ``dfloat``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``uinterp``, ``vinterp``, ``winterp``               | ``dfloat``                   | interpolated velocity values |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``scalarId``                                        | ``int``                      |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``s``                                               | ``dfloat``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``flux``                                            | ``dfloat``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``sinterp``                                         | ``dfloat``                   | interpolated scalar value    |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``meshu``, ``meshv``, ``meshw``                     | ``dfloat``                   |                              |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``trans``, ``diff``                                 | ``dfloat``                   | properties                   |
+-----------------------------------------------------+------------------------------+------------------------------+
| ``usrwrk``                                          | ``@globalPtr const dfloat*`` |                              |
+-----------------------------------------------------+------------------------------+------------------------------+

Internal / Periodic
"""""""""""""""""""

None is used when a internal boundary condition is required or a periodic 
boundary condition has been set as part of the mesh and it does not need to be
considered as part of the standard processing of boundary conditions.

Perodicity is linked to the mesh connectivity and is handled by the meshing tool.

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

