*******************
Components Overview
*******************

.. note::

   If you are doing simple tasks where performance is a non-critical aspect please go to the :ref:`High-Level APIs` section for a quick start.
   If you are an HPC application developer or you want to use ADIOS2 functionality in full please read this chapter.


The simple way to understand the big picture for the ADIOS2 unified user interface components is to map each class to the actual definition of the ADIOS acronym.

+------------+-----------+---------------------------+
| Component  | Acronym   | Function                  |
+------------+-----------+---------------------------+
|            |           | Set MPI comm domain       |
|            |           |                           |
| **ADIOS**  | ADaptable | Set runtime settings      |
|            |           |                           |
|            |           | Own other components      |
+------------+-----------+---------------------------+
|            |           | Set engine                |
|            |           |                           |
| **IO**     | I/O       | Set variables/attributes  |
|            |           |                           |
|            |           | Set compile-time settings |
+------------+-----------+---------------------------+
|            |           | Execute heavy IO tasks    |
| **Engine** | System    |                           |
|            |           | Manage system resources   |
+------------+-----------+---------------------------+


ADIOS2's public APIs are based on the natural choice for each supported language to represent each ADIOS2 components and its interaction with application datatypes. Thus,


============== ========================== ==================================
 **Language**      **Component API**       **Application Data**
============== ========================== ==================================
 C++(11/newer)  objects/member functions    pointers/references/std::vector
 C              handler/functions           pointers
 Fortran        handler/subroutines         arrays up to 6D
 Python         objects/member functions    numpy arrays.
============== ========================== ==================================

The following section provides a common overview to all languages based on the C++11 APIs. For each specific language go to the :ref:`Full APIs` section, but it's highly recommended to read this section as components map 1-to-1 in other languages.

The following figure depicts the components hierarchy from the application's point of view.

.. image:: https://i.imgur.com/y7bkQQt.png

* **ADIOS**: the ADIOS component is the starting point between an application and the ADIOS2 library. Applications provide:
    1. the scope of the ADIOS object through the MPI communicator,
    2. an optional runtime configuration file (in XML format) to allow changing settings without recompiling.

    The ADIOS component serves as a factory of adaptable IO components. Each IO must have a unique name within the scope of the ADIOS class object that created them with the DeclareIO function.

* **IO**: the IO component is the bridge between the application specific settings, transports. It also serves as a factory of:
    1. Variables
    2. Attributes
    3. Engines

* **Variable**: Variables are the link between self-describing representation in the ADIOS2 library and data from applications. Variables are identified by unique names in the scope of the particular IO that created them. When the Engine API functions are called, a Variable must be provided along with the application data.

* **Attribute**: Attributes add extra information to the overall variables dataset defined in the IO class. They can be single or array values.

* **Engine**: Engines define the actual system executing the heavy IO tasks at Open, BeginStep, Put, Get, EndStep and Close. Due to polymorphism, new IO system solutions can be developed quickly reusing internal components and reusing the same API. If IO.SetEngine is not called, the default engine is the binary-pack bp file reader and writer: **BPFile**.

* **Operator**: These define possible operations to be applied on adios2-managed data, for example, compression. This higher level abstraction is needed to provide support for callbacks, transforms, analytics, data models, etc. Any required task will be executed within the Engine. One or many operators can be associated with any of the adios2 objects or a group of them.
