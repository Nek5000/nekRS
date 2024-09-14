*******************
Python Example Code
*******************

The Python APIs follows closely Python style directives. They rely on numpy and, optionally, on ``mpi4py``, if the underlying ADIOS2 library is compiled with MPI.

For online examples on MyBinder :

- `Python-MPI Notebooks <https://mybinder.org/v2/gh/ornladios/ADIOS2-Jupyter.git/python-mpi>`_

- `Python-noMPI Notebooks <https://mybinder.org/v2/gh/ornladios/ADIOS2-Jupyter.git/python-nompi>`_


Examples in the ADIOS2 repository
---------------------------------

- Simple file-based examples
    - examples/hello/helloWorld/hello-world.py
    - examples/hello/bpReader/bpReaderHeatMap2D.py
    - examples/hello/bpWriter/bpWriter.py

- Staging examples using staging engines SST and DataMan
    - examples/hello/sstWriter/sstWriter.py
    - examples/hello/sstReader/sstReader.py
    - examples/hello/datamanWriter/dataManWriter.py
    - examples/hello/datamanReader/dataManReader.py

Python Write example
--------------------

.. code-block:: python
   
   from mpi4py import MPI
   import numpy as np
   from adios2 import Stream
   
   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()
   size = comm.Get_size()
   
   nx = 10
   shape = [size * nx]
   start = [rank * nx]
   count = [nx]
   
   temperature = np.zeros(nx, dtype=np.double)
   pressure = np.ones(nx, dtype=np.double)
   delta_time = 0.01
   physical_time = 0.0
   nsteps = 5
   
   with Stream("cfd.bp", "w", comm) as s:
      # NSteps from application
      for _ in s.steps(nsteps):
         if rank == 0 and s.current_step() == 0:
            # write a Python integer
            s.write("nproc", size)
         
         # write a Python floating point value
         s.write("physical_time", physical_time)
         # temperature and pressure are numpy arrays
         s.write("temperature", temperature, shape, start, count)
         s.write_attribute("temperature/unit", "K")
         s.write("pressure", pressure, shape, start, count)
         s.write_attribute("pressure/unit", "Pa")
         physical_time += delta_time

.. code-block:: bash

    $ mpirun -n 4 python3 ./adios2-doc-write.py
    $ bpls -la cfd.bp 
      int64_t  nproc             scalar = 4
      double   physical_time     5*scalar = 0 / 0.04
      double   pressure          5*{40} = 1 / 1
      string   pressure/unit     attr   = "Pa"
      double   temperature       5*{40} = 0 / 0
      string   temperature/unit  attr   = "K"


Python Read "step-by-step" example
----------------------------------

.. code-block:: python
   
    import numpy as np
    from adios2 import Stream

    with Stream("cfd.bp", "r") as s:
        # steps comes from the stream
        for _ in s.steps():

            # track current step
            print(f"Current step is {s.current_step()}")

            # inspect variables in current step
            for name, info in s.available_variables().items():
                print("variable_name: " + name, end=" ")
                for key, value in info.items():
                    print("\t" + key + ": " + value, end=" ")
                print()

            if s.current_step() == 0:
                nproc = s.read("nproc")
                print(f"nproc is {nproc} of type {type(nproc)}")

            # read variables return a numpy array with corresponding selection
            physical_time = s.read("physical_time")
            print(f"physical_time is {physical_time} of type {type(physical_time)}")
            temperature = s.read("temperature")
            temp_unit = s.read_attribute("temperature/unit")
            print(f"temperature array size is {temperature.size} of shape {temperature.shape}")
            print(f"temperature unit is {temp_unit} of type {type(temp_unit)}")
            pressure = s.read("pressure")
            press_unit = s.read_attribute("pressure/unit")
            print(f"pressure unit is {press_unit} of type {type(press_unit)}")
            print()

.. code-block:: bash

    $ python3 adios2-doc-read.py
    Current step is 0
    variable_name: nproc    AvailableStepsCount: 1  Max: 4  Min: 4  Shape:          SingleValue: true       Type: int64_t
    variable_name: physical_time    AvailableStepsCount: 1  Max: 0  Min: 0  Shape:          SingleValue: true       Type: double
    variable_name: pressure         AvailableStepsCount: 1  Max: 1  Min: 1  Shape: 40       SingleValue: false      Type: double
    variable_name: temperature      AvailableStepsCount: 1  Max: 0  Min: 0  Shape: 40       SingleValue: false      Type: double
    nproc is 4 of type <class 'numpy.ndarray'>
    physical_time is 0.0 of type <class 'numpy.ndarray'>
    temperature array size is 40 of shape (40,)
    temperature unit is K of type <class 'str'>
    pressure unit is Pa of type <class 'str'>

    Current step is 1
    variable_name: physical_time    AvailableStepsCount: 1  Max: 0.01       Min: 0.01       Shape:          SingleValue: true   Type: double
    variable_name: pressure         AvailableStepsCount: 1  Max: 1  Min: 1  Shape: 40       SingleValue: false      Type: double
    variable_name: temperature      AvailableStepsCount: 1  Max: 0  Min: 0  Shape: 40       SingleValue: false      Type: double
    physical_time is 0.01 of type <class 'numpy.ndarray'>
    temperature array size is 40 of shape (40,)
    temperature unit is K of type <class 'str'>
    pressure unit is Pa of type <class 'str'>

    ...


Python Read Random Access example
----------------------------------

.. code-block:: python

    import numpy as np
    from adios2 import FileReader

    with FileReader("cfd.bp") as s:
        # inspect variables
        vars = s.available_variables()
        for name, info in vars.items():
            print("variable_name: " + name, end=" ")
            for key, value in info.items():
                print("\t" + key + ": " + value, end=" ")
            print()
        print()

        nproc = s.read("nproc")
        print(f"nproc is {nproc} of type {type(nproc)} with ndim {nproc.ndim}")
        
        # read variables return a numpy array with corresponding selection
        steps = int(vars['physical_time']['AvailableStepsCount'])
        physical_time = s.read("physical_time", step_selection=[0, steps])
        print(
            f"physical_time is {physical_time} of type {type(physical_time)} with "
            f"ndim {physical_time.ndim} shape = {physical_time.shape}"
        )

        steps = int(vars['temperature']['AvailableStepsCount'])
        temperature = s.read("temperature", step_selection=[0, steps])
        temp_unit = s.read_attribute("temperature/unit")
        print(f"temperature array size is {temperature.size} of shape {temperature.shape}")
        print(f"temperature unit is {temp_unit} of type {type(temp_unit)}")

        steps = int(vars['pressure']['AvailableStepsCount'])
        pressure = s.read("pressure", step_selection=[0, steps])
        press_unit = s.read_attribute("pressure/unit")
        print()

.. code-block:: bash

    $ python3 adios2-doc-read-filereader.py
    variable_name: nproc    AvailableStepsCount: 1  Max: 4  Min: 4  Shape:          SingleValue: true       Type: int64_t
    variable_name: physical_time    AvailableStepsCount: 5  Max: 0.04       Min: 0  Shape:          SingleValue: true       Type: double
    variable_name: pressure         AvailableStepsCount: 5  Max: 1  Min: 1  Shape: 40       SingleValue: false      Type: double
    variable_name: temperature      AvailableStepsCount: 5  Max: 0  Min: 0  Shape: 40       SingleValue: false      Type: double

    nproc is 4 of type <class 'numpy.ndarray'> with ndim 0
    physical_time is [0.   0.01 0.02 0.03 0.04] of type <class 'numpy.ndarray'> with ndim 1 shape = (5,)
    temperature array size is 200 of shape (200,)
    temperature unit is K of type <class 'str'>

