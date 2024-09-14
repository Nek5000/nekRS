**********
C bindings
**********

The C bindings are specifically designed for C applications and those using an older C++ standard (98 and 03). If you are using a C++11 or more recent standard, please use the C++11 bindings.

The C bindings are based on opaque pointers to the components described in the :ref:`Components Overview` section. 

C API handlers:

.. code-block:: c
   
   adios2_adios*
   adios2_io*
   adios2_variable*
   adios2_attribute*
   adios2_engine*
   adios2_operator* 


Every ADIOS2 function that generates a new ``adios2_*`` unique handler returns the latter explicitly.
Therefore, checks can be applied to know if the resulting handler is NULL.
Other functions used to manipulate these valid handlers will return a value of type ``enum adios2_error`` explicitly.
These possible errors are based on the `C++ standardized exceptions <https://en.cppreference.com/w/cpp/error/exception>`_ .
Each error will issue a more detailed description in the standard error output: ``stderr``. 

``adios2_error`` possible values:

.. code-block:: C

   typedef enum {
       /** success */
       adios2_error_none = 0,

       /**
        * user input error
        */
       adios2_error_invalid_argument = 1,
   
       /** low-level system error, e.g. system IO error */
       adios2_error_system_error = 2,
   
       /** runtime errors other than system errors, e.g. memory overflow */
       adios2_error_runtime_error = 3,
   
       /** any other error exception */
       adios2_error_exception = 4
   } adios2_error; 


Usage:

.. code-block:: C

   adios2_variable* var = adios2_define_variable(io, ...)
   if(var == NULL )
   {
       // ... something went wrong with adios2_define_variable
       // ... check stderr
   }
   else
   {
       adios2_type type;
       adios2_error errio = adios2_variable_type(&type, var)
       if(errio){
         // ... something went wrong with adios2_variable_type
         if( errio == adios2_error_invalid_argument)
         {
             // ... user input error
             // ... check stderr
         }
       }
   }


.. note::
    
    Use ``#include "adios2_c.h"`` for the C bindings, ``adios2.h`` is the C++11 header.


``adios2_adios`` handler functions
----------------------------------

.. doxygenfile:: adios2_c_adios.h
   :project: C
   :path: ../../bindings/C/adios2/c/

``adios2_io`` handler functions
-------------------------------

.. doxygenfile:: adios2_c_io.h
   :project: C
   :path: ../../bindings/C/adios2/c/

``adios2_variable`` handler functions
-------------------------------------

.. doxygenfile:: adios2_c_variable.h
   :project: C
   :path: ../../bindings/C/adios2/c/

``adios2_attribute`` handler functions
--------------------------------------

.. doxygenfile:: adios2_c_attribute.h
   :project: C
   :path: ../../bindings/C/adios2/c/

``adios2_engine`` handler functions
-----------------------------------

.. doxygenfile:: adios2_c_engine.h
   :project: C
   :path: ../../bindings/C/adios2/c/

``adios2_operator`` handler functions
-------------------------------------

.. doxygenfile:: adios2_c_operator.h
   :project: C
   :path: ../../bindings/C/adios2/c/
