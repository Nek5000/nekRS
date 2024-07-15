.. _running:

Running nekRS
=============

This page gives information on how to run nekRS has been installed 
(see :ref:`installing` ) and appropriate input files have been generated 
(see :ref:`input`).

Basic Running
-------------

Native Running
""""""""""""""

The most basic way of running nekRS is in the directory with all the appropriate 
files with the scenario name:

.. code-block::

    nekrs --setup <scenario>.par

Running with MPI
""""""""""""""""

The majority of users will want to run with MPI enabled to utilise multiple CPU 
cores. If nek has been installed with this enabled this can be done with:

.. code-block::

    mpirun -np <number of MPI tasks> nekrs --setup <scenario>.par

There are also helper scripts that can be used to run nekRS with MPI and
optionally in the background. See :ref:`nekrs_scripts` for more information.

.. code-block::

    nrsmpi <scenario>.par <number of MPI tasks>
    nrsbmpi <scenario>.par <number of MPI tasks>

Cluster Running
"""""""""""""""

TODO

Command Line arguments
""""""""""""""""""""""

Below are the command line arguments that can be used to further modifiy how 
nekRS is run.

+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
|    Parameter     | Short option |                       Options                       |                                           Description                                            | Required |
+==================+==============+=====================================================+==================================================================================================+==========+
| ``--help``       | ``-h``       | None or ``par``                                     | Print help, either summary of command line argument                                              | No       |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
| ``--setup``      | ``-s``       | None, ``par`` or ``sess file``                      | Specifies the location of files needed to initialise the simulation                              | Yes      |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
| ``--build-only`` | ``-b``       | None or ``#procs``                                  | Initialise the simulation and run :ref:`just_in_time_compilation` only.                          | No       |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
| ``--cimode``     | ``-c``       | None or ``<id>`` (N.B. If set must be >=0)          | Runs specific CI tests if available in the chosen simulation (see :ref:`contributing`)           | No       |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
| ``--debug``      | ``-d``       | None                                                | Run in debug mode, IE print values of many of the variables while running                        | No       |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
| ``--backend``    | ``-t``       | ``CPU``, ``CUDA``, ``HIP``, ``DPCPP`` or ``OPENCL`` | Manually set the backend device for running                                                      | No       |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+
| ``--device-id``  | ``-i``       | ``id`` or ``LOCAL-RANK``                            | Manually set OCCA device ID (I.E. for machines with multiple GPUs) or use ``LOCAL-RANK`` for CPU | No       |
+------------------+--------------+-----------------------------------------------------+--------------------------------------------------------------------------------------------------+----------+

.. _nekrs_scripts:

nekRS helper scripts
--------------------

A number of scripts ship with nekRS itself and are located in the 
``$NEKRS_HOME/bin`` directory (see :ref:`nekrs_home`). A brief summary of these 
scripts and their usage is as follows.

* ``nrsmpi <casename> <processes>``: run nekRS in parallel with ``<processes>`` parallel
  processes for the case files that are prefixed with ``casename``.
* ``nrsbmpi <casename> <processes>``: same as ``nrsmpi``, except that nekRS runs
  in the background
* ``nrspre <casename> <target GPUs>``: precompile nekRS case (see
  :ref:`just_in_time_compilation`)
* ``nrsqsub_lassen <casename> <nodes> <wall time>``: submission script for
  `Lassen <https://computing.llnl.gov/computers/lassen>`_, a supercomputer
  at Lawrence Livermore National Laboratory. A number of other settings are specified
  within the script itself.
* ``nrsqsub_summit <casename> <nodes> <wall time>``: submission script for
  `Summit <https://www.olcf.ornl.gov/summit/>`_, a supercomputer
  at Oak Ridge National Laboratory. A number of other settings are specified within the
  script itself.
* ``nrsvis <casename>``: postprocess ``fld``-type nekRS output files into a form
  readable by Paraview or Visit.