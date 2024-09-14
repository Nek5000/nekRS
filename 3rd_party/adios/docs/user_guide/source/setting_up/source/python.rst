****************************
Enabling the Python bindings
****************************

To enable the Python bindings in ADIOS2, based on `PyBind11 <http://pybind11.readthedocs.io/en/stable/>`_, make sure to follow these guidelines:

- **Minimum requirements:**

    * Python 2.7 and above.
    * `numpy`
    * `mpi4py`

- **Running:** If CMake enables Python compilation, an ``adios2.so`` library containing the Python module is generated in the build directory under ``lib/pythonX.X/site-packages/``

    * make sure your ``PYTHONPATH`` environment variable contains the path to ``adios2.so``.

    * make sure the Python interpreter is compatible with the version used for compilation via ``python --version``.

    * Run the Python tests with ``ctest -R Python``

    * Run `helloBPWriter.py <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpWriter/helloBPWriter.py>`_ and `helloBPTimeWriter.py <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/bpTimeWriter/helloBPTimeWriter.py>`_ via

    .. code-block:: bash

            $ mpirun -n 4 python helloBPWriter.py
            $ python helloBPWriter.py
