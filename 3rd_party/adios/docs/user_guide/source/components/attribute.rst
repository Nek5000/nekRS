*********
Attribute
*********

Attributes are extra information associated with a particular IO component.
They can be thought of as a very simplified ``Variable``, but with the goal of adding extra metadata.
The most common use is the addition of human-readable metadata (*e.g.* ``"experiment name"``, ``"date and time"``, ``"04,27,2017"``, or a schema).

Currently, ADIOS2 supports single values and arrays of primitive types (excluding ``complex<T>``) for the template type in the ``IO::DefineAttribute<T>`` and ``IO::InquireAttribute<T>`` function (in C++).

The data types supported for ADIOS2 ``Attributes`` are

.. code-block:: c++

   std::string
   char
   signed char  
   unsigned char  
   short  
   unsigned short  
   int  
   unsigned int  
   long int  
   long long int  
   unsigned long int  
   unsigned long long int  
   float  
   double  
   long double 
   
The returned object (``DefineAttribute`` or ``InquireAttribute``) only serves the purpose to inspect the current ``Attribute<T>`` information within code.

.. note:

   Attributes are not forcibly associated to a particular variable in ADIOS2.
   Developers are free to create associations through their own naming conventions.
