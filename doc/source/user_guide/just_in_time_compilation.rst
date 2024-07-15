.. _just_in_time_compilation:

Just-in-time Compilation
========================

nekRS uses just-in-time compilation to build the functions in the ``.udf`` and ``.oudf`` 
case files, as well as for compiling certain fixed-size arrays based on the order of
the polynomial approximation or other problem settings.

For most cases, no special actions need to be taken by the user for this
process to work correctly. However, a high-level understanding of the just-in-time
compilation is useful to know what steps need to be taken to fully clear the cached
build files, as well as how to perform the pre-compilation separately from a full run
to obtain more accurate runtime measurements.

When nekRS performs just-in-time compilation, object files are created in the
``.cache`` directory within the current case directory. To completely clear the
cached state of the case, simply delete the ``.cache`` directory:

.. code-block::

  user$ rm -rf .cache/

.. tip::

   If you experience strange behavior when running your case during the precompilation
   step (such as failures to build in COMMON blocks or other parts of the code that you
   are not touching in the ``.udf`` and ``.oudf`` files), try deleting the ``.cache``
   directory and trying again. It is not uncommon for the precompilation process to miss
   the need to build new versions of object files if you are making frequent changes to
   the nekRS source. This is also sometimes encountered if you are using multiple nekRS
   versions in different projects (such as standalone nekRS or nekRS wrapped within
   a multiphysics coupling application such as :term:`ENRICO`), but don't have your
   environment completely self-consistent.

The precompilation process usually takes on the order of one minute. Depending on
the use case, it may be advantageous to force the precompilation separately from the run itself.
To precompile the case, use the ``nrspre`` script. See the
:ref:`Scripts That Ship with nekRS <nekrs_scripts>` section for where to find this script.

As an example, to precompile a case with name ``my_case`` for a later run with less than
or equal to 4 GPUs, the ``nrspre`` script would be used as follows. After the precompilation,
run as normal with the ``nrsmpi`` script, which then skips the precompilation step since
the case has already been compiled.

.. code-block::

  user$ nrspre my_case 4
  user$ nrsmpi my_case 32
