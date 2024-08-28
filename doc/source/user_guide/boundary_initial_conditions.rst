.. _boundary_initial_conditions:

Boundary and Initial conditions
===============================

.. _setting_ICs:

Setting Initial Conditions with ``UDF_Setup``
---------------------------------------------

This section provides an example for setting initial conditions with the
``UDF_Setup`` user-defined function that was introduced on the :ref:`Input Files <input>` page.
The following code snippet sets initial conditions for all three components of
velocity, the pressure, and two passive scalars. You may not necessarily have all of these
variables in your model - this example is just intended to cover all possibilities.

For this example, the initial conditions are
:math:`V_x=sin(x)cos(y)cos(z)`, :math:`V_y=-cos(x)sin(y)cos(z)`, and :math:`V_z=0`
for the three components of velocity;
:math:`P=101325` for the pressure; and :math:`\phi_0=573` and :math:`\phi_1=100+z` for the
two passive scalars indicated generically as :math:`\phi_0` and :math:`\phi_1`.

.. note::

  If present, the temperature variable is represented internally in nekRS as a passive
  scalar, since the form of the equation is the same as those solver for other passive
  scalars such as chemical concentration.

Because these initial conditions will
be a function of space, we must first obtain the mesh information, for which we
use the ``nrs->mesh`` pointer. All solution fields are stored in nekRS in terms of the
quadrature points (also referred to as the :term:`GLL` points). So, we will apply
the initial conditions by looping over all of these quadrature points, which for
the current :term:`MPI` process is equal to ``mesh->Np``, or the number of quadrature
points per element, and ``mesh->Nelements``, the number of elements on this process.

Next, we can get the :math:`x`, :math:`y`, and :math:`z` coordinates for the current
quadrature point with the ``x``, ``y``, and ``z`` pointers on the ``mesh`` object.
Finally, we programmatically set initial conditions for the solution fields. ``nrs->U``
is a single array that holds all three components of velocity; the ``nrs->fieldOffset``
variable is used to shift between components in this array. ``nrs->P`` represents the
pressure. Finally, ``nrs->S`` is a single array that holds all of the passive scalars.
Similar to the offset performed to index into the velocity array, the
``nrs->cds->fieldOffset`` variable is used to shift between components in the ``nrs->S``
array.

.. code-block:: cpp

   void UDF_Setup(nrs_t* nrs)
   {
    mesh_t* mesh = nrs->mesh;
    int num_quadrature_points = mesh->Np * mesh->Nelements;

    for (int n = 0; n < num_quadrature_points; n++) {
      float x = mesh->x[n];
      float y = mesh->y[n];
      float z = mesh->z[n];

      nrs->U[n + 0 * nrs->fieldOffset] = sin(x) * cos(y) * cos(z);
      nrs->U[n + 1 * nrs->fieldOffset] = -cos(x) * sin(y) * cos(z);
      nrs->U[n + 2 * nrs->fieldOffset] = 0;

      nrs->P[n] = 101325.0;

      nrs->S[n + 0 * nrs->cds->fieldOffset] = 573.0;
      nrs->S[n + 1 * nrs->cds->fieldOffset] = 100.0 + z;
    }
   }


Periodic Boundary Conditions
----------------------------

NekRS supports periodic boundary conditions. To set up a periodic case, first
you need to run ``exo2nek`` to establish the pairings between the periodic sidesets.
All this information will be prompted on the screen by ``exo2nek``;
You will provide the sideset IDs of the periodic boundaries, a search tolerance
for identifying paired sides, and a translation vector that points from one of the
paired sidesets to the other. For example, if you want to have one periodic surface
that is a :math:`z`-plane at :math:`z=-1.0` that is paired to another :math:`z`-plane
at :math:`z=1.0`, the translation vector would be :math:`(0.0, 0.0, 2.0)`.

After generating the mesh, you then need to modify the sideset IDs inside the
``usrdat2`` function. Any boundary that is now periodic, you need to set
``boundaryID(ifc,iel)`` to 0. For all non-periodic boundaries, you need to
"renormalize" those boundaries to "begin counting" from 1. For example, consider
an original (non-periodic) mesh with sidesets 1, 2, 3, and 4. You run ``exo2nek``
and set up sidesets 2 and 3 as periodic. Then, in the code snippet below, you
would reset sidesets 2 and 3 in ``boundaryID`` to zero. For the remaining two
boundaries (originally 1 and 4), you need to renormalized those to boundaries
1 and 2 (because NekRS wants the boundaries to be ordered sequentially beginning
from 1).

.. code-block::

      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'
      integer e,f

      n = lx1*ly1*lz1*nelv
      nxz = nx1*nz1
      nface = 2*ldim

      do iel=1,nelt
      do ifc=1,2*ndim
         if (boundaryID(ifc,iel).eq. 1) then
           boundaryID(ifc,iel) = 1
         else if (boundaryID(ifc,iel).eq. 2) then
           boundaryID(ifc,iel) = 0
         else if (boundaryID(ifc,iel) .eq. 3) then
           boundaryID(ifc,iel) = 0
         else if (boundaryID(ifc,iel) .eq. 4) then
           boundaryID(ifc,iel) = 2
         endif
      enddo
      enddo

      return
      end

Then, in the other case files, you do not need any boundary conditions for the periodic
boundaries - for instance, in the ``<case>.par`` file for this example, the boundary conditions
set in ``boundaryTypeMap`` would only display the boundary conditions for the non-periodic
boundaries (and similarly in the ``<case>.oudf`` file). Finally, in order to enforce periodic
flow with a constant flow rate, specify the ``constFlowRate`` parameter in the ``<case>.par``
file, such as

.. code-block::

    [GENERAL]
      constFlowRate = meanVelocity=1.0 + direction=Z

Stamping Initial Conditions
---------------------------

For many periodic flows, you can save significant computing time by solving the flow equations
on a shorter-height mesh, and then "stamping" that solution onto a full-height mesh (where you
might then be solving for passive scalar transport). NekRS allows you to "stamp" a partial-height
solution onto a full-height mesh using the ``gfldr`` utility. To do so, you simply need to call
the ``gfldr`` function in a loop inside of ``userchk()``. Below, ``nd`` represents the number
of times you want to stamp a short-height solution to obtain the full-height case and ``delta``
represents the height of one short-height domain. So, the example below would represent
a previous solution (``short.fld``) on a short-height domain of height 62.42, that you want to stamp five times
onto a new mesh that has a height of 312.1.

.. code-block::

      subroutine userchk()
      include 'SIZE'
      include 'TOTAL'


      ntot = lx1*ly1*lz1*nelv

      do nd = 0,5

        delta = 62.421731741003335

        do i = 1,ntot
         zm1(i,1,1,1) = zm1(i,1,1,1) - delta*nd
        enddo

        call gfldr('short.fld')

        do i = 1,ntot
         zm1(i,1,1,1) = zm1(i,1,1,1) + delta*nd
        enddo

      enddo

      return
      end

Velocity Recycling Plugin
-------------------------

