###############
High-Level APIs
###############

The high-level APIs are designed for simple tasks for which performance is not critical. Unlike the :ref:`Full APIs`, the high-level APIs only require a single object handler resembling a C++ fstream or a Python file I/O idiom. The high-level APIs are recommended to both first-time and advanced users; the low-level APIs being recommended only when performance testing identifies a bottleneck or when more control is needed.

Typical scenarios for using the simple high-level APIs are:

* Reading a file to perform data analysis with libraries (matplotlib, scipy, etc.)
* Interactive: few calls make interactive usage easier. 
* Saving data to files is small or personal projects
* Online frameworks: *e.g.* Jupyter notebooks, see python-mpi examples running on `MyBinder <https://mybinder.org/v2/gh/ornladios/ADIOS2-Jupyter.git/python-mpi>`_ 

The designed functionality syntax is closely related to the native language IO bindings for formatted text files *e.g.* C++ ``fstream`` ``getline``, and Python file IO.
The main function calls are: ``open`` (or constructor in C++), ``write``, ``read`` and ``close`` (or destructor in C++).
In addition, ADIOS2 borrows the corresponding language native syntax for advancing lines to advance the step in write mode, and for a "step-by-step" streaming basis in read mode.
See each language section in this chapter for a write/read example.

.. note::

   The simplified APIs are based on language native file IO interface. Hence ``write`` and ``read`` calls are always synchronized and variables data memory is ready to use immediately after these calls.


Currently ADIOS2 support bindings for the following languages and their minimum standards:

+----------+----------+-----------------------+-------------+
| Language | Standard | Interface             | Based on    |
+----------+----------+-----------------------+-------------+
| C++      | 11/newer | ``#include adios2.h`` | ``fstream`` |
+----------+----------+-----------------------+-------------+
| Matlab   |          |                       |             |
+----------+----------+-----------------------+-------------+

The following sections provide a summary of the API calls on each language and links to Write and Read examples to put it all together.

.. include:: cxx11.rst
.. include:: matlab.rst

