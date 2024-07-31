.. _tutorial_rans:

------------
RANS Channel
------------

This tutorial describes the essential setup details for a wall-resolved RANS simulation, illustrated through a 3D channel case. 
The :math:`k-\tau` RANS model is employed for this tutorial, which is the recommended RANS model in *nekRS*. 
For more information see :ref:`here <intro_ktau>`. 

..........................
Before You Begin
..........................

It is highly recommended that new users familiarize themselves with the basic *nekRS* simulation
setup files and procedures outlined in the :ref:`fdlf` and :ref:`perhill` tutorials before proceeding.

..............................
Mesh and Boundary Conditions
..............................
 
The mesh is generated with ``genbox`` using the following input file

.. literalinclude:: ktauTutorial/input.box
   :language: none 
   
It creates an infinite 3D half-channel of non-dimensional width :math:`1`. 
The streamwise (:math:`x`) direction has 5  elements with periodic (``P``) boundary conditions.
The wall-normal (:math:`y`) direction has 12  geometrically spaced elements with a symmetry (``SYM``) boundary condition specified at the bottom face and a wall (``W``) boundary on the top face. 
The spanwise (:math:`z`) direction has 4 elements with periodic (``P``) boundary conditions.
The elements expand from the bottom of the domain to the top with a ratio of :math:`0.35`.

.. In addition to velocity and temperature, RANS simulation requires two additional fields for the turbulent scalars, in this cases turbulent kinetic energy, :math:`k`, and specific dissipation period, :math:`\tau`. The boundary conditions for these scalars are assigned in the :ref:`.par`.

.. Note::
  Boundary conditions for the turbulent scalars are specified in the :ref:`.udf` file.

..............................
Control Parameters (.par file)
..............................

Details of the structure of the parameter file can be found :ref:`here <case_files_par>`. 
For RANS simulations it is critical to include two additional scalars which correspond to the :math:`k` and 
:math:`\tau` fields respectively. 
Here, they are included as the ``[SCALAR01]`` and ``[SCALAR02]`` cards.
In addition, it is essential to also include the ``[PROBLEMTYPE]`` card and enable ``variableViscosity``.

For this particular tutorial, the simulation is :ref:`non-dimensionalized <intro_ns_nondim>` and flow properties are ``density=1.0`` and ``viscosity=-43500``. 
It is strongly recommended to run RANS simulations in non-dimensional form. 
By assigning a negative value to the viscosity, *nekRS* will treat it as a Reynolds number.
This setup corresponds to :math:`Re = 43,500`, which is well above the critical Reynolds number for a channel.
``density`` and ``diffusivity``  for ``SCALAR01`` and ``SCALAR02`` should be assigned identical values as ``density`` and ``viscosity`` for velocity field for consistency.

.. literalinclude:: ktauTutorial/channel.par
   :language: ini 
   :emphasize-lines: 33-46

The Temperature field is also  solved in this tutorial, but can be turned off  by uncommenting  the ``Solver = None`` line on the ``.par`` file.

..............................
User Routines (.usr file)
..............................

This section describes the essential code snippets required in the ``.usr`` file for RANS simulations.
Other details of all the subroutines can be found :ref:`here <case_files_usr>`. 

Foremost, it is essential to include the following header at the beginning of the ``.usr`` file.

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 2
      
Files in the above relative locations in the *nekRS* repo load the essential RANS subroutines.

.. _rans_init:

RANS initialization is done through the ``rans_init`` subroutine call from ``usrdat2``. 
The required code snippet is shown below.

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 162-198
   
The ``rans_init`` subroutine sets up the necessary solver variables to use RANS models. 
This includes loading the model coefficients, setting up the character boundary condition (``cbc``) array for the turbulent scalars, and calculating the regularization [Tombo2018]_ for the :math:`k-\omega` models.
The ``ifld_tke`` and ``ifld_tau`` variables specify the field index location of the transport variables of the two-equation RANS model. 
The specific RANS model used is identified by the ``m_id`` variable. 
All available RANS models are annotated in the above code. 
The recommended :math:`k-\tau` model is invoked with ``m_id=4``.
Selecting ``ifcoeffs=.false.`` sets the standard RANS coefficients outlined in [Wilcox2008]_. 
Advanced users have the option of specifying their own coefficients which have to populated in the ``coeffs`` array, with ``ifcoeffs=.true.``. 
The final parameter to be aware of is the wall distance function code ``w_id``. 
The recommended value for it is ``w_id=2`` which provides a smoother distance function, populated in the ``wd`` array. 
A cheaper option is available through ``w_id=1`` which provides "taxicab" distance and can search across periodic boundaries. 
The user also has the option of specifying their own wall distance array by setting ``w_id=0`` which will require ``wd`` array to
be populated with user computed wall distance before the ``rans_init`` call. 

Diffusion coefficients for all fields in RANS simulation runs must be modified to include eddy viscosity.
This is done by the following inclusions in the ``uservp`` subroutine

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 5-39
	
As above, eddy viscosity, ``mu_t``, is added to the momentum diffusion coefficient and :math:`k` and :math:`\tau` diffusion coefficients are modified as described by Eq. :eq:`ktau` in the :ref:`RANS theory section<intro_ktau>`. 

.. Note:: 
  The transport and diffusion coefficients (``utrans`` and ``udiff``) for the RANS scalars are set with the :ref:`field coefficient array <tab:cpfld>` vaules (``cpfld``) for field 1, velocity. This is done to ensure consistency in the problem setup.

Source terms in the :math:`k-\tau` transport equations, Eq. :eq:`ktau`, are added with the following inputs in the ``userq`` subroutine

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 61-88
	
Implicit source terms contributions are specified to the ``avol`` variable, while remaining terms are 
in the ``qvol`` variable. 

For the wall-resolved :math:`k-\tau` RANS model, the wall boundary conditions for :math:`k` and :math:`\tau` are both zero. 
These are set in ``userbc`` as

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 90-119
   :emphasize-lines: 21

The velocity boundary condition is checked in ``cb1`` for a wall on the highlighted line.
For this periodic case, this is not strictly necessary. 
However, for any inlet/outlet case it is needed to distinguish between the Dirichlet BC values for :math:`k` and :math:`\tau` on the wall versus an inlet, where ``cb1`` would have a value of ``v``.
	
Initial conditions are specified in the ``useric`` routine. 
For RANS simulations, positive non-zero initial values for :math:`k` and :math:`\tau` are recommended. 
The following is used for the channel simulation,

.. literalinclude:: rans/wallResolved/chan_WR.usr
   :language: fortran
   :lines: 121-143
	
..............................
SIZE file
..............................

The ``SIZE`` file can be copied from the available template included at the ``nekRS/core/SIZE.template`` and changing the basic parameters as necessary for this case.
The user needs to ensure that the auxiliary fields specified in the SIZE file is at minimum ``ldimt=3`` for RANS simulations. 
This includes the reserved field for temperature and the two additional fields for the turbulent scalars.
In order to use the stress formulation, additional memory must be allocated by setting, ``lx1m=lx1``
These two changes are highlighted below.
Other details on the contents of the ``SIZE`` file can be found :ref:`here<case_files_SIZE>`.

.. literalinclude:: rans/wallResolved/SIZE
   :language: fortran
   :lines: 11-40
   :emphasize-lines: 10,24

..............................
Compilation
..............................

All required case files for RANS wall-resolved channel simulation can be downloaded using the links below:

 * :download:`chan_WR.usr <rans/wallResolved/chan_WR.usr>`
 * :download:`chan_WR.par <rans/wallResolved/chan_WR.par>`
 * :download:`chan_WR.re2 <rans/wallResolved/chan_WR.re2>`
 * :download:`chan_WR.ma2 <rans/wallResolved/chan_WR.ma2>`
 * :download:`SIZE <rans/wallResolved/SIZE>`
 
..............................
Results
..............................

For reference, the results obtained from the :math:`k-\tau` RANS wall-resolved simulation are shown below and
compared with results from low-Re :math:`k-\tau` (``m_id=5``), regularized :math:`k-\omega` (``m_id=0``) and low-Re regularized 
:math:`k-\omega` (``m_id=1``) models.

.. _fig:streamwise_vel:

.. figure:: rans/U1.png
   :align: center
   :figclass: align-center

   Normalized stream-wise velocity from different RANS models.
   
.. _fig:chan_tke:

.. figure:: rans/k1.png
   :align: center
   :figclass: align-center

   Normalized TKE from different RANS models.
   
RANS models may be simply switched by using the appropriate ``m_id`` in ``usrdat2``. 
The low-Re versions of the model should be employed only for capturing the near wall TKE peak. 
It may require additional resolution in the near-wall region and has only minimal impact on the velocity profile.
