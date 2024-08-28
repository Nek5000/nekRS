.. _postprocessing:

Postprocessing
==============

.. _writing_output:

Writing an Output File - using UDF_ExecuteStep
----------------------------------------------

nekRS will automatically write output files according to the ``writeControl`` criterion
set in the ``.par`` file. However, it may be desirable to have finer-grained control of
output writing, such as if you want the solution at a specific time step, but that
time step is not an integer multiple of ``writeInterval``. In this case, you can force
the output file writing to occur by calling the ``outfld(double time, double outputTime)``
function in the ``nekrs`` namespace. This function performs the following actions:

1. Copy the nekRS solution from the nekRS device arrays directly to the backend
   Nek5000 arrays.
2. Write an output file.

Note that this function is slightly different from the ``nek_ocopyFrom`` function described
in the :ref:`Copying Device to Host <copy_device_to_host>` section. This function is
solely intended for writing output, so no effort is expended in copying the device
solution into the nekRS host arrays - that step is bypassed, and the device solution is
copied straight into the Nek5000 backend arrays. The ``nek_ocopyFrom`` routine should really
only be used if you require access to the nekRS solution arrays on the host, while the
``outfld`` routine should be used strictly for writing output files.

By default, nekRS will only write the velocity, pressure, and temperature to an output file.
However, you may have problem-specific fields that you want to view, such as :math:`y^+`.
To write other fields to files, nekRS re-uses the
functions that are used to write the velocity, pressure, and temperature
to write other fields. Note that this imposes limitations on both the dimensionality of fields that
can be output, as well as how they are named in the output files.
For example, suppose you would like to write three fields to a file:

  * ``o_yPlus``, a device array that holds :math:`y^+` values, and
  * ``o_Uavg``, a device array that holds a time-averaged velocity field, and
  * ``o_rst``, a device array that holds the one component of the Reynolds stress tensor.

To write these three fields to an output file, use the ``writeFld`` function as follows.
The ``writeFld`` function takes eight arguments, and has a signature
``void writeFld(const char* suf, dfloat t, int coords, int FP64, void* o_u, void* o_p, void* o_s, int NSf)``.
In this example, the first parameter, ``"usr"``, is a three-character
prefix that will determine how the new output file is written. While the velocity, pressure,
and temperatures are written to files named ``case0.f<time_step>``, where ``case`` is the case
name and ``<time_step>`` is a six-digit number indicating the time step, any additional fields
we will write are written to separate files. So for this example, we will write three fields
to files named ``usrcase0.f<time_step>``. The next three parameters simply indicate the time
step that is being written, whether coordinates are written, and if the results should be written
in double precision. Next, the three fields that are to be output are provided. The order is very
important - the first of these fields must be of length ``nrs->fieldOffset * nrs->NVfields``
because it represents a component vector field (this is how velocity is written in the usual output
file). The second of these fields must be of length ``nrs->fieldOffset``, because it represents
a non-component field (this is how pressure is written in the usual output file). Finally,
the third of these fields must be of length ``nrs->cds->fieldOffset * Nscalar``, because it
represents a passive scalar field (this is how the passive scalars are written in the usual
output file).

.. code-block:: cpp

   void UDF_ExecuteStep(nrs_t* nrs, dfloat time, int tstep)
   {
     // get o_yPlus, o_Uavg, and o_rst in the scope of this function

     bool coords = true;
     bool FP64 = true;
     int Nscalar = nrs->cds->Nscalar;
     writeFld("usr", time, coords, FP64, &o_Uavg, &o_rst, &o_yPlus, Nscalar);
   }

.. warning::

  ``writeFld`` can only write data of type ``dfloat``. So, if you want to write an
  integer field to a field, you must first convert that data to ``dfloat``.

nekRS's output system does not have any means by which to understand *what* these fields
represent. Therefore, the names of these fields in the output file will be ``velocity``,
``pressure``, and ``temperature``, even if those names have no relationship to what is
being output. Therefore, for this example, the ``usrcase0.f<time_step>`` files will
contain the following:

* ``o_Uavg`` is written to a field named ``velocity``
* ``o_rst`` is written to a field named ``pressure``
* ``o_yPlus`` is written to a field named ``temperature``

nekRS's output system requires additional maneuvering if you wish to output
more than one of each of each of these three categories of fields. For instance, suppose
you want to output three different fields, ``o_field1``, ``o_field2``, and ``o_field3``,
each of size ``nrs->fieldOffset``. Because only one input argument to ``writeFld`` can have
these dimensions, three separate output files need to be written, and in *each* of these
files, our field of interest is named ``pressure``. To fill the other two field arguments
of the ``writeFld`` function, a void pointer is passed in to indicate that neither of
the other two fields are written.

.. code-block:: cpp

   void UDF_ExecuteStep(nrs_t* nrs, dfloat time, int tstep)
   {
     // get o_field1, o_field2, o_field3 in the scope of this function

     bool coords = true;
     bool FP64 = true;
     int Nscalar = nrs->cds->Nscalar;
     occa::memory o_null;
     writeFld("fl1", time, coords, FP64, &o_null, &o_field1, &o_null, Nscalar);
     writeFld("fl2", time, coords, FP64, &o_null, &o_field2, &o_null, Nscalar);
     writeFld("fl3", time, coords, FP64, &o_null, &o_field3, &o_null, Nscalar);
   }

This will write three output files, which contain the following.

* ``fl1case0.f<time_step>`` contains ``o_field1``, but named ``pressure``
* ``fl2case0.f<time_step>`` contains ``o_field2``, but named ``pressure``
* ``fl3case0.f<time_step>`` contains ``o_field3``, but named ``pressure``

.. _vis_output:

Visualizing Output Files
------------------------

nekRS output files all have the form ``<case0>.fld<n>``, where ``<case>`` is the case
name and ``<n>`` is a five-digit number indicating the number of the output file (each output
file represents a single time step that is output according to the settings for
``writeControl`` and ``writeInterval`` in the ``.par`` file). These output files are in a custom
binary format that requires an additional postprocessing step in order to visualize in Paraview.
In the directory where the case files are located, run the ``visnek`` script:

.. code-block::

  user$ visnek case

which will create a ``case.nek5000`` file that is viewable in Paraview. See
:ref:`Building the Nek5000 Tool Scripts <scripts>` for instructions on compiling the ``visnek`` program.

Turbulence Statistics
"""""""""""""""""""""