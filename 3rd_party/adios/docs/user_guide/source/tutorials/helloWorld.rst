Hello World
===========

.. _sec:tutorials_basics_hello_world:

Like in any language, the first program you write is the "Hello World" program.

In this tutorial, we will see how to write "Hello World from ADIOS2" and read it back with ADIOS2.
So let's dig in!

Start editing the skeleton file `ADIOS2/examples/hello/helloWorld/hello-world_tutorialSkeleton.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/helloWorld/hello-world_tutorialSkeleton.cpp>`_.


1. We create an ADIOS instance, and define the greeting message in our main function as follows:

.. code-block:: c++

  int main()
  {
    adios2::ADIOS adios();
    const std::string greeting("Hello World from ADIOS2");
    ...
    return 0;
  }

2. Then we create a writer function in which we pass the adios instance, and the greeting message as follows:

.. code-block:: c++

  void writer(adios2::ADIOS& adios, const std::string& greeting)
  {
    ...
  }

3. In this writer function, we define an IO object, a string variable for the message as follows:

.. code-block:: c++

  adios2::IO io = adios.DeclareIO("hello-world-writer");
  adios2::Variable<std::string> varGreeting = io.DefineVariable<std::string>("Greeting");

.. note::

  Using the IO object, we can define the engine type that we want to utilize using the *io.SetEngine()* function.
  If *SetEngine()* is not used, the default engine type is *BPFile* which is an alias for the latest version of the BP
  engine of the ADIOS2 library. See :ref:`Available Engines` and :ref:`Supported Engines` for more information.
  It's important to note that the file extension of an output file, although it's not a good practice, it can differ
  from the engine type, e.g. write a foo.h5 file with the BPFile engine. When reading foo.h5 you should explicitly
  specify the engine type as BPFile to read it properly.

4. Then we open a file with the name "hello-world-cpp.bp" and write the greeting message to it as follows:

.. code-block:: c++

  adios2::Engine writer = io.Open("hello-world-cpp.bp", adios2::Mode::Write);
  writer.BeginStep();
  writer.Put(varGreeting, greeting);
  writer.EndStep();
  writer.Close();

.. note::

  The ``BeginStep`` and ``EndStep`` calls are optional when **writing** one step, but they are required
  for multiple steps, so it is a good practice to always use them.

5. Now we create a reader function in which we pass the adios instance, and get the greeting message back as follows:

.. code-block:: c++

  std::string reader(adios2::ADIOS& adios)
  {
    ...
  }

6. In this reader function, we define an IO object and inquire a string variable for the message as follows:

.. code-block:: c++

  adios2::IO io = adios.DeclareIO("hello-world-reader");
  reader.BeginStep();
  adios2::Variable<std::string> varGreeting = io.InquireVariable<std::string>("Greeting");

7. Then we open the file with the name "hello-world-cpp.bp", read the greeting message from it and return it as follows:

.. code-block:: c++

  adios2::Engine reader = io.Open("hello-world-cpp.bp", adios2::Mode::Read);
  std::string greeting;
  reader.BeginStep();
  reader.Get(varGreeting, greeting);
  reader.EndStep();
  reader.Close();
  return greeting;

.. note::

  In Mode::Read, the ``BeginStep`` and ``EndStep`` calls are required when **reading** one step and multiple steps. We will see in
  another tutorial how to read multiple steps. It's important to note that the ``BeginStep`` should be called **before**
  all ``Inquire*`` / ``Available*`` function calls.

8. Finally, we call the writer and reader functions in our main function as follows:

.. code-block:: c++

  int main()
  {
    adios2::ADIOS adios();
    const std::string greeting("Hello World from ADIOS2");
    writer(adios, greeting);
    std::string message = reader(adios);
    std::cout << message << std::endl;
    return 0;
  }

9. The final code should look as follows (excluding try/catch and the optional usage of MPI), and it was derived from
   the example `ADIOS2/examples/hello/helloWorld/hello-world.cpp <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/helloWorld/hello-world.cpp>`_.

.. literalinclude:: ../../../../examples/hello/helloWorld/hello-world.cpp
   :language: cpp

10. You can compile and run it as follows:

.. code-block:: bash

   cd Path-To-ADIOS2/examples/hello/helloWorld
   mkdir build
   cd build
   cmake -DADIOS2_DIR=Path-To-ADIOS2/build/ ..
   cmake --build .
   ./adios2_hello_helloWorld

11. You can check the content of the output file "hello-world-cpp.bp" using *bpls* as follows:

.. code-block:: bash

   Path-To-ADIOS2/build/bin/bpls ./hello-world-cpp.bp

      string   Greeting  scalar

12. The Python version of this tutorial can be found at `ADIOS2/examples/hello/helloWorld/hello-world.py <https://github.com/ornladios/ADIOS2/blob/master/examples/hello/helloWorld/hello-world.py>`_.
    and it looks as follows:

.. literalinclude:: ../../../../examples/hello/helloWorld/hello-world.py
   :language: python
