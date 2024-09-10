.. _meshing:

Meshing
=======

.. _converting_mesh:

Converting a Mesh to .re2 Format
--------------------------------

The most general and flexible approach for creating a mesh is to use commercial meshing software
such as Cubit or Gmsh. After creating the mesh, it must be converted to the ``.re2`` binary format.
The following two sections describe how to convert Exodus and Gmsh meshes into ``.re2`` binary format
with scripts that ship with the Nek5000 dependency. First build these scripts following
the instructions in the :ref:`Building the Nek5000 Tool Scripts <scripts>` section.

Converting an Exodus mesh
"""""""""""""""""""""""""

To convert from an Exodus format mesh
(for this case, named ``my_mesh.exo``) to the ``.re2`` format, use the ``exo2nek`` script:

.. code-block::

  user$ exo2nek

Then, follow the on-screen prompts associated with the ``exo2nek`` script.
``exo2nek`` will convert all elements in the Exodus mesh (TET6, WEDGE6, HEX8, HEX20) to HEX20 elements and dump into ``.re2`` format.

Converting a Gmsh mesh
""""""""""""""""""""""

To convert from a Gmsh format mesh (for this case, named ``my_mesh.msh``) to the
``.re2`` format, use the ``gmsh2nek`` script:

.. code-block::

  user$ gmsh2nek

  Enter mesh dimension: 3
  Input (.msh) file name: my_mesh


All your mesh should be hexahedral elements. Before exporting from Gmsh, you will need to set the mesh order to 2.
The Gmsh mesh format must also be version 2, ASCII/binary format. If your Gmsh version
shows a pop-up box when exporting the mesh, do *not* click "Save all elements"
or "Save parametric elements".

.. _periodic_mesh:

Setting up Mesh for Periodic Boundary Condition
"""""""""""""""""""""""""""""""""""""""""""""""

NekRS supports periodic boundary conditions. To set up a periodic case, first
you need to run ``exo2nek`` to establish the pairings between the periodic sidesets.
All this information will be prompted on the screen by ``exo2nek``;
You will provide the sideset IDs of the periodic boundaries, a search tolerance
for identifying paired sides, and a translation vector that points from one of the
paired sidesets to the other. For example, if you want to have one periodic surface
that is a :math:`z`-plane at :math:`z=-1.0` that is paired to another :math:`z`-plane
at :math:`z=1.0`, the translation vector would be :math:`(0.0, 0.0, 2.0)`.

.. _cht_mesh:

Creating a Mesh for Conjugate Heat Transfer
-------------------------------------------

Mesh generation for conjugate heat transfer requires an additional pre-processing
step before performing other steps of the mesh generation process such as those
described in the :ref:`Converting a Mesh to .re2 Format <converting_mesh>` section.
The nekRS approach for conjugate heat transfer is still dependent on legacy limitations
from Nek5000. Unfortunately, you cannot
simply use a standard commercial meshing tool and define fluid and solid
regions according to block IDs - you must individually create the mesh for the fluid and
the solid, and then merge them with the ``pretex`` script.