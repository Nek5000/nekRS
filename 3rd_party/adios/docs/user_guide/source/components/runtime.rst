***************************
Runtime Configuration Files
***************************

ADIOS2 supports passing an optional runtime configuration file to the :ref:`ADIOS` component constructor (``adios2_init`` in C, Fortran).

This file contains key-value pairs equivalent to the compile time ``IO::SetParameters`` (``adios2_set_parameter`` in C, Fortran), and ``IO::AddTransport`` (``adios2_set_transport_parameter`` in C, Fortran). 

Each ``Engine`` and ``Operator`` must provide a set of available parameters as described in the :ref:`Supported Engines` section.
Prior to version v2.6.0 only XML is supported; v2.6.0 and later support both XML and YAML.

.. warning::

   Configuration files must have the corresponding format extension ``.xml``, ``.yaml``: ``config.xml``, ``config.yaml``, etc.

XML
---

.. code-block:: xml

   <?xml version="1.0"?>
   <adios-config>
     <io name="IONAME_1">  

       <engine type="ENGINE_TYPE"> 
            
         <!-- Equivalent to IO::SetParameters--> 
         <parameter key="KEY_1" value="VALUE_1"/>
         <parameter key="KEY_2" value="VALUE_2"/>
         <!-- ... -->
         <parameter key="KEY_N" value="VALUE_N"/> 
        
       </engine>

       <!-- Equivalent to IO::AddTransport -->
       <transport type="TRANSPORT_TYPE">
         <!-- Equivalent to IO::SetParameters--> 
         <parameter key="KEY_1" value="VALUE_1"/>
         <parameter key="KEY_2" value="VALUE_2"/>
         <!-- ... -->
         <parameter key="KEY_N" value="VALUE_N"/>
       </transport>
     </io>
         
     <io name="IONAME_2">  
       <!-- ... -->
     </io>
   </adios-config>
            
           
YAML
----

Starting with v2.6.0, ADIOS supports YAML configuration files.
The syntax follows strict use of the YAML node keywords mapping to the ADIOS2 components hierarchy.
If a keyword is unknown ADIOS2 simply ignores it.
For an example file refer to `adios2 config file example in our repo. <https://github.com/ornladios/ADIOS2/tree/master/testing/adios2/yaml/proto.yaml>`_

.. code-block:: yaml
   
   ---
   # adios2 config.yaml
   # IO YAML Sequence (-) Nodes to allow for multiple IO nodes
   # IO name referred in code with DeclareIO is mandatory
   
   - IO: "IOName"   
     
     Engine:
        # If Type is missing or commented out, default Engine is picked up
        Type: "BP5"
        # optional engine parameters
        key1: value1
        key2: value2
        key3: value3
     
     Variables:

         # Variable Name is Mandatory
       - Variable: "VariableName1"
         Operations:
             # Operation Type is mandatory (zfp, sz, etc.)
           - Type: operatorType
             key1: value1
             key2: value2
       
       - Variable: "VariableName2"
         Operations:
             # Operations sequence of maps
           - {Type: operatorType, key1: value1}
           - {Type: z-checker, key1: value1, key2: value2}
      
     Transports: 
         # Transport sequence of maps 
       - {Type: file, Library: fstream}
       - {Type: rdma, Library: ibverbs}
        
     ...


.. caution::
   
   YAML is case sensitive, make sure the node identifiers follow strictly the keywords: IO, Engine, Variables, Variable, Operations, Transports, Type.
   
.. tip::
   
   Run a YAML validator or use a YAML editor to make sure the provided file is YAML compatible.
