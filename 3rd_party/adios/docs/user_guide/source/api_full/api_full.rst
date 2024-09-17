#########
Full APIs
#########

.. note::

   Application developers who desire fine-grained control of IO tasks should use the Full APIs.
   In simple cases (e.g. reading a file for analysis, interactive Python, or saving some data for a small project or tests) please refer to the :ref:`High-Level APIs`.


Currently ADIOS2 support bindings for the following languages and their minimum standards:

+----------+----------+-------------------------+
| Language | Standard | Interface               |
+----------+----------+-------------------------+
|          | 11/newer | ``#include adios2.h``   |
| C++      |          |                         |
|          | older    | use C bindings          |
+----------+----------+-------------------------+
| C        | 99       | ``#include adios2_c.h`` |
+----------+----------+-------------------------+
| Fortran  | 90       | ``use adios2``          |
+----------+----------+-------------------------+
|          | 2.7      | ``import adios2``       |
| Python   |          |                         |
|          | 3        | ``import adios2``       |
+----------+----------+-------------------------+

.. tip::

   Prefer the C++11 bindings if your application C++ compiler supports the 2011 (or later) standard.
   For code using previous C++ standards (98 or 03) use the C bindings for ABI compatibility.


.. caution::

   Statically linked libraries (``*.a``) might result in conflicting ABIs between an older C++ project, the C bindings, and the adios native C++11 library. Test to make sure it works for your platform.


The current interaction flow for each language binding API with the ADIOS2 library is specified as follows

.. blockdiag ::
   
   diagram {
      default_fontsize = 22;
      default_shape = roundedbox;
      default_linecolor = blue;
   
      "ADIOS2" -> "C++11", C, Python;
      C -> Fortran, Matlab; 
      
      "ADIOS2"[width = 200, height = 60]
      
      group{
         label = "C++11 Library"
         color = orange
        
         "ADIOS2";
      }
      
      group{
         label = "Language bindings APIs"
         color = yellow
         
         "C++11", C, Python, Fortran, Matlab;
      }
   }

The following sections provide a summary of the API calls on each language and links to Write and Read examples to put it all together.

.. include:: cxx11.rst
.. include:: fortran.rst
.. include:: c.rst
