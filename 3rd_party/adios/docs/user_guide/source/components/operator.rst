********
Operator
********

.. _sec:basics_interface_components_operator:

The Operator abstraction allows ADIOS2 to act upon the user application data, either from a ``adios2::Variable`` or a set of Variables in an ``adios2::IO`` object.
Current supported operations are:

1. Data compression/decompression, lossy and lossless.
2. Callback functions (C++11 bindings only) supported by specific engines

ADIOS2 enables the use of third-party libraries to execute these tasks.

Operators can be attached onto a variable in two modes: private or shared. In most situations, it is recommended to add an operator as a private one, which means it is owned by a certain variable. A simple example code is as follows.

.. code-block:: c++

    #include <vector>
    #include <adios2.h>
    int main(int argc, char *argv[])
    {
        std::vector<double> myData = {
            0.0001, 1.0001, 2.0001, 3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001,
            1.0001, 2.0001, 3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001,
            2.0001, 3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001,
            3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001,
            4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001,
            5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001,
            6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001,
            7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001, 2.0001,
            8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001, 2.0001, 1.0001,
            9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001, 2.0001, 1.0001, 0.0001,
        };
        adios2::ADIOS adios;
        auto io = adios.DeclareIO("TestIO");
        auto varDouble = io.DefineVariable<double>("varDouble", {10,10}, {0,0}, {10,10}, adios2::ConstantDims);

        // add operator
        varDouble.AddOperation("mgard",{{"accuracy","0.01"}});
        // end add operator

        auto engine = io.Open("hello.bp", adios2::Mode::Write);
        engine.Put<double>(varDouble, myData.data());
        engine.Close();
        return 0;
    }

For users who need to attach a single operator onto multiple variables, a shared operator can be defined using the adios2::ADIOS object, and then attached to multiple variables using the reference to the operator object. Note that in this mode, all variables sharing this operator will also share the same configuration map. It should be only used when a number of variables need *exactly* the same operation. In real world use cases this is rarely seen, so please use this mode with caution.

.. code-block:: c++

    #include <vector>
    #include <adios2.h>
    int main(int argc, char *argv[])
    {
        std::vector<double> myData = {
            0.0001, 1.0001, 2.0001, 3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001,
            1.0001, 2.0001, 3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001,
            2.0001, 3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001,
            3.0001, 4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001,
            4.0001, 5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001,
            5.0001, 6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001,
            6.0001, 7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001,
            7.0001, 8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001, 2.0001,
            8.0001, 9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001, 2.0001, 1.0001,
            9.0001, 8.0001, 7.0001, 6.0001, 5.0001, 4.0001, 3.0001, 2.0001, 1.0001, 0.0001,
        };
        adios2::ADIOS adios;
        auto io = adios.DeclareIO("TestIO");
        auto varDouble = io.DefineVariable<double>("varDouble", {10,10}, {0,0}, {10,10}, adios2::ConstantDims);

        // define operator
        auto op = adios.DefineOperator("SharedCompressor","mgard",{{"accuracy","0.01"}});
        // add operator
        varDouble.AddOperation(op);
        // end add operator

        auto engine = io.Open("hello.bp", adios2::Mode::Write);
        engine.Put<double>(varDouble, myData.data());
        engine.Close();
        return 0;
    }

.. warning::

   Make sure your ADIOS2 library installation used for writing and reading was linked with a compatible version of a third-party dependency when working with operators.
   ADIOS2 will issue an exception if an operator library dependency is missing.
