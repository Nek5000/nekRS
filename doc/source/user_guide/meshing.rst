.. _meshing:

Meshing
=======

The first step in setting up a case is to have a mesh of the geometry to correctly
constrain the simulation and detail the locations at which boundary conditions will
later be applied. Simple meshes can be created  using the ``genbox`` tool, but the
most general and flexible approach is to use commercial meshing software such as
Cubit or Gmsh. NekRS requires the subsequent mesh to be in the ``.re2`` binary 
format which can be generated (through conversions scripts) from common mesh formats
such as CGNS (CFD General Notation System), EXO (Exodus II) or GMSH.

The ``genbox``, ``cgns2nek``, ``exo2nek`` and ``gmsh2nek`` tools are installed 
by following the instructions in the :ref:`scripts` section. We recommend looking
at the The README's for the `cgns2nek <https://github.com/Nek5000/Nek5000/blob/master/tools/cgns2nek/README.md>`__,
`exo2nek <https://github.com/Nek5000/Nek5000/blob/master/tools/exo2nek/README.md>`__ or
`gmsh2nek <https://github.com/Nek5000/Nek5000/blob/master/tools/gmsh2nek/README.md>`__ tools
if using them, with further information also available in the 
`Nek5000 documentation <http://nek5000.github.io/NekDoc/tools.html>`__.

.. _cubit_mesh:

Mesh with Cubit
---------------

TODO Example of creating Exodus Mesh

.. _gmsh_mesh:

Meshing with Gmsh
-----------------

TODO Example of creating a Gmsh mesh

All your mesh should be hexahedral elements. Before exporting from Gmsh, you will need to set the mesh order to 2.
The Gmsh mesh format must also be version 2, ASCII/binary format. If your Gmsh version
shows a pop-up box when exporting the mesh, do *not* click "Save all elements"
or "Save parametric elements".

.. _cgns_mesh:

Creating a mesh with ?????? - CGNS 
----------------------------------

TODO Example of creating a cgns mesh

.. _cht_mesh:

Conjugate Heat Transfer
-----------------------

