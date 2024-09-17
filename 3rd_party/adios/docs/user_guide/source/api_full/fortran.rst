****************
Fortran bindings
****************

.. role:: f90(code)
   :language: fortran
   :class: highlight

The Fortran API is a collection of subroutine calls. The first argument is usually a Fortran type (struct) to an ADIOS2 component, while the last argument is an error integer flag, :f90:`integer ierr`. ``ierr==0`` represents successful execution whereas a non-zero value represents an error or a different state. ADIOS2 Fortran bindings provide a list of possible errors coming from the C++ standardized error exception library:

.. code-block:: fortran

    ! error possible values for ierr
    integer, parameter :: adios2_error_none = 0
    integer, parameter :: adios2_error_invalid_argument = 1,
    integer, parameter :: adios2_error_system_error = 2,
    integer, parameter :: adios2_error_runtime_error = 3,
    integer, parameter :: adios2_error_exception = 4

Click here for a `Fortran write and read example`_ to illustrate the use of the APIs calls. This test will compile under your build/bin/ directory.

.. _`Fortran write and read example`: https://github.com/ornladios/ADIOS2/blob/master/testing/adios2/bindings/fortran/TestBPWriteReadHeatMap3D.F90

The following subsections describe the overall components and subroutines in the Fortran bindings API.

ADIOS2 typed handlers
---------------------

ADIOS2 Fortran bindings handlers are mapped 1-to-1 to the ADIOS2 components described in the :ref:`Components Overview` section. For convenience, each type handler contains descriptive components used for read-only inspection.
 
.. code-block:: fortran

   type(adios2_adios)     :: adios
   type(adios2_io)        :: io
   type(adios2_variable)  :: variable
   type(adios2_attribute) :: attribute
   type(adios2_engine)    :: engine
   
   !Read-only components for inspection and ( = defaults)
   
   type adios2_adios
        logical :: valid = .false.
    end type

    type adios2_io
        logical :: valid = .false.
        character(len=15):: engine_type = 'BPFile'
    end type

    type adios2_variable
        logical :: valid = .false.
        character(len=4095):: name = ''
        integer :: type = -1
        integer :: ndims = -1
    end type

    type adios2_attribute
        logical :: valid = .false.
        character(len=4095):: name = ''
        integer :: type = -1
        integer :: length = 0
    end type

    type adios2_engine
        logical :: valid = .false.
        character(len=63):: name = ''
        character(len=15):: type = ''
        integer :: mode = adios2_mode_undefined
    end type

    type adios2_operator
        logical :: valid = .false.
        character(len=63):: name = ''
        character(len=63):: type = ''
    end type
   

.. caution::

   Use the type read-only components for information purposes only.
   Changing their values directly, *e.g.* ``variable%name = new_name`` does not have any effect inside the ADIOS2 library 
   

:ref:`ADIOS` subroutines
------------------------

* :f90:`subroutine adios2_init` starting point for the ADIOS2 library 

   .. code-block:: fortran

      ! MPI versions
      subroutine adios2_init(adios, comm, ierr)
      subroutine adios2_init(adios, config_file, comm, ierr)
      
      ! Non-MPI serial versions
      subroutine adios2_init(adios, ierr)
      subroutine adios2_init(adios, config_file, ierr) 
   
      ! WHERE:
      
      ! ADIOS2 handler to allocate
      type(adios2_adios), intent(out):: adios 
      
      ! MPI Communicator
      integer, intent(in):: comm 
      
      ! Optional runtime configuration file (*.xml), see Runtime Configuration Files
      character*(*), intent(in) :: config_file
      
      ! error code
      integer, intent(out) :: ierr
      

* :f90:`subroutine adios2_declare_io` spawn/create an IO component

   .. code-block:: fortran

      subroutine adios2_declare_io(io, adios, io_name, ierr)
      
      ! WHERE:
      
      ! output ADIOS2 IO handler
      type(adios2_io), intent(out):: io
      
      ! ADIOS2 component from adios2_init spawning io tasks 
      type(adios2_adios), intent(in):: adios
      
      ! unique name associated with this IO component inside ADIOS2
      character*(*), intent(in):: io_name

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_at_io` retrieve an existing io component. Useful when the original IO handler goes out of scope

   .. code-block:: fortran

      subroutine adios2_at_io(io, adios, io_name, ierr)
      
      ! WHERE:
      
      ! output IO handler
      type(adios2_io), intent(out):: io
      
      ! ADIOS2 component from adios2_init that owns IO tasks 
      type(adios2_adios), intent(in):: adios
      
      ! unique name associated with an existing IO component (created with adios2_declare_io)
      character*(*), intent(in):: io_name

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_define_operator` define an ADIOS2 data compression/reduction operator

   .. code-block:: fortran

      subroutine adios2_define_operator(op, adios, op_name, op_type, ierr)

      ! WHERE

      ! Operator handler
      type(adios2_operator), intent(out) :: op

      ! ADIOS2 handler
      type(adios2_adios), intent(in) :: adios

      ! Operator name
      character*(*), intent(in)  :: op_name
      
      ! Operator type
      character*(*), intent(in)  :: op_type

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_inquire_operator` inquire an ADIOS2 data compression/reduction operator

   .. code-block:: fortran

      subroutine adios2_inquire_operator(op, adios, op_name, ierr)

      ! WHERE

      ! Operator handler
      type(adios2_operator), intent(out) :: op

      ! ADIOS2 handler
      type(adios2_adios), intent(in) :: adios

      ! Operator name
      character*(*), intent(in)  :: op_name

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_flush_all` flush all current engines in all IO objects

   .. code-block:: fortran

      subroutine adios2_flush_all(adios, ierr)
      
      ! WHERE:
      
      ! ADIOS2 component from adios2_init owning IO objects and engines 
      type(adios2_adios), intent(in):: adios

      ! error code
      integer, intent(out) :: ierr
   

* :f90:`subroutine adios2_remove_io` DANGER ZONE: remove an IO object. This will effectively eliminate any parameter from the config xml file

   .. code-block:: fortran

      subroutine adios2_remove_io(result, adios, name, ierr)

      ! WHERE

      ! Returns True if IO was found, False otherwise
      logical, intent(out):: result

      ! ADIOS2 handler
      type(adios2_adios), intent(in) :: adios

      ! IO input name
      character*(*), intent(in):: name

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_remove_all_ios` DANGER ZONE: remove all IO objects created for this ADIOS2 handler. This will effectively eliminate any parameter from the config xml file as well.

   .. code-block:: fortran

      subroutine adios2_remove_all_ios(adios, ierr)

      ! WHERE

      ! ADIOS2 handler
      type(adios2_adios), intent(in) :: adios

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_finalize` final point for the ADIOS2 component

   .. code-block:: fortran

      subroutine adios2_finalize(adios, ierr)
      
      ! WHERE:
      
      ! ADIOS2 handler to be deallocated 
      type(adios2_adios), intent(in):: adios

      ! error code
      integer, intent(out) :: ierr


.. caution::
   
   Make sure that for every call to ``adios2_init`` there is a call to ``adios2_finalize`` for the same ADIOS2 handler. Not doing so will result in memory leaks. 


* :f90:`subroutine adios2_enter_computation_block` inform ADIOS2 about entering communication-free computation block in main thread. Useful when using Async IO.

   .. code-block:: fortran

      subroutine adios2_enter_computation_block(adios, ierr)

      ! WHERE

      ! ADIOS2 handler
      type(adios2_adios), intent(in) :: adios

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_exit_computation_block` inform ADIOS2 about exiting communication-free computation block in main thread. Useful when using Async IO.

   .. code-block:: fortran

      subroutine adios2_exit_computation_block(adios, ierr)

      ! WHERE

      ! ADIOS2 handler
      type(adios2_adios), intent(in) :: adios

      ! error code
      integer, intent(out) :: ierr

      
:ref:`IO` subroutines
---------------------
   
* :f90:`subroutine adios2_set_engine` set the engine type, see :ref:`Supported Engines` for a list of available engines
   
   .. code-block:: fortran
      
      subroutine adios2_set_engine(io, engine_type, ierr)
      
      ! WHERE:
      
      ! IO component
      type(adios2_io), intent(in):: io
      
      ! engine_type: BP (default), HDF5, DataMan, SST, SSC
      character*(*), intent(in):: engine_type

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_in_config_file` checks if an IO object exists in a config file passed to ADIOS2.

   .. code-block:: fortran

      subroutine adios2_in_config_file(result, io, ierr)

      ! WHERE

      ! Output result to indicate whether IO exists
      logical, intent(out):: result

      ! IO handler
      type(adios2_io), intent(in):: io

      ! error code
      integer, intent(out):: ierr


* :f90:`subroutine adios2_set_parameter` set IO key/value pair parameter in an IO object, see :ref:`Supported Engines` for a list of available parameters for each engine type
   
   .. code-block:: fortran
      
      subroutine adios2_set_parameter(io, key, value, ierr)
      
      ! WHERE:
      
      ! IO component owning the attribute
      type(adios2_io), intent(in):: io
      
      ! key in the key/value pair parameter
      character*(*), intent(in):: key
      
      ! value in the key/value pair parameter
      character*(*), intent(in):: value

      ! error code
      integer, intent(out) :: ierr
      

* :f90:`subroutine adios2_set_parameters` set a map of key/value parameters in an IO object. Replaces any existing parameters. Otherwise use set_parameter for adding single parameters.

   .. code-block:: fortran

      subroutine adios2_set_parameters(io, parameters, ierr)

      ! WHERE

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! Comma-separated parameter list. E.g. "Threads=2, CollectiveMetadata=OFF"
      character*(*), intent(in) :: parameters

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_get_parameter` get parameter value from IO object for a given parameter name

   .. code-block:: fortran

      subroutine adios2_get_parameter(value, io, key, ierr)

      ! WHERE

      ! parameter value
      character(len=:), allocatable, intent(out) :: value

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! parameter key to look for in the IO object
      character*(*), intent(in) :: key

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_clear_parameters` clear all parameters from the IO object

   .. code-block:: fortran

      subroutine adios2_clear_parameters(io, ierr)

      ! WHERE

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_add_transport` add a transport to current IO. Must be supported by the currently used engine.

   .. code-block:: fortran

      subroutine adios2_add_transport(transport_index, io, type, ierr)

      ! WHERE

      ! returns a transport_index handler
      integer, intent(out):: transport_index

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! transport type. must be supported by the engine. CAN’T use the keywords “Transport” or “transport”
      character*(*), intent(in) :: type

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_set_transport_parameter` set a parameter for a transport. Overwrites existing parameter with the same key.

   .. code-block:: fortran

      subroutine adios2_set_transport_parameter(io, transport_index, key, value, ierr)

      ! WHERE

      ! IO handler
      type(adios2_io), intent(in):: io

      ! transport_index handler
      integer, intent(in):: transport_index

      ! transport key
      character*(*), intent(in) :: key

      ! transport value
      character*(*), intent(in) :: value

      ! error code
      integer, intent(out):: ierr


* :f90:`subroutine adios2_available_variables` get a list of available variables

   .. code-block:: fortran

      subroutine adios2_available_variables(io, namestruct, ierr)

      ! WHERE

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! name struct handler 
      type(adios2_namestruct), intent(out) :: namestruct

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_retrieve_names` retrieve variable names from namestruct obtained from ``adios2_available_variables``. namelist must be pre-allocated.

   .. code-block:: fortran

      subroutine adios2_retrieve_names(namestruct, namelist, ierr)

      ! WHERE

      ! namestruct obtained from adios2_available_variables
      type(adios2_namestruct), intent(inout) :: namestruct

      ! namelist that will contain variable names
      character(*), dimension(*), intent(inout) :: namelist

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_available_attributes` get list of attributes in the IO object

   .. code-block:: fortran

      subroutine adios2_available_attributes(io, namestruct, ierr)

      ! WHERE

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! list of available attributes
      type(adios2_namestruct), intent(out) :: namestruct

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_flush_all_engines` flush all existing engines opened by this IO object
   
   .. code-block:: fortran
   
      subroutine adios2_flush_all_engines(io, ierr)
        
      ! WHERE:
      
      ! IO in which search and flush for all engines is performed
      type(adios2_io), intent(in) :: io 

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_remove_variable` remove an existing variable from an IO object
   
   .. code-block:: fortran
   
      subroutine adios2_remove_variable(io, name, result, ierr)
        
      ! WHERE:
      
      ! IO in which search and removal for variable is performed
      type(adios2_io), intent(in) :: io
      
      ! unique key name to search for variable 
      character*(*), intent(in) :: name
      
      ! true: variable removed, false: variable not found, not removed
      logical, intent(out) :: result

      ! error code
      integer, intent(out) :: ierr
      

* :f90:`subroutine adios2_remove_all_variables` remove all existing variables from an IO object
   
   .. code-block:: fortran
   
      subroutine adios2_remove_all_variables(io, ierr)
        
      ! WHERE:
      
      ! IO in which search and removal for all variables is performed
      type(adios2_io), intent(in) :: io

      ! error code
      integer, intent(out) :: ierr
         
      
* :f90:`subroutine adios2_remove_attribute` remove existing attribute by its unique name
   
   .. code-block:: fortran
   
      subroutine adios2_remove_attribute(io, name, result, ierr)
        
      ! WHERE:
      
      ! IO in which search and removal for attribute is performed
      type(adios2_io), intent(in) :: io
      
      ! unique key name to search for attribute 
      character*(*), intent(in) :: name
      
      ! true: attribute removed, false: attribute not found, not removed
      logical, intent(out) :: result

      ! error code
      integer, intent(out) :: ierr
         
      
* :f90:`subroutine adios2_remove_all_attributes` remove all existing attributes
   
   .. code-block:: fortran
   
      subroutine adios2_remove_all_attributes(io, ierr)
        
      ! WHERE:
      
      ! IO in which search and removal for all attributes is performed
      type(adios2_io), intent(in) :: io

      ! error code
      integer, intent(out) :: ierr


:ref:`Variable` subroutines
---------------------------
     
* :f90:`subroutine adios2_define_variable` define/create a new variable

   .. code-block:: fortran

      ! Global array variables
      subroutine adios2_define_variable(variable, io, variable_name, adios2_type, &
                                        ndims, shape_dims, start_dims, count_dims, & 
                                        adios2_constant_dims, ierr) 
      ! Global single value variables
      subroutine adios2_define_variable(variable, io, variable_name, adios2_type, ierr)
      
      ! WHERE:
      
      ! handler to newly defined variable
      type(adios2_variable), intent(out):: variable
      
      ! IO component owning the variable
      type(adios2_io), intent(in):: io
      
      ! unique variable identifier within io
      character*(*), intent(in):: variable_name
      
      ! defines variable type from adios2 parameters, see next 
      integer, intent(in):: adios2_type 
      
      ! number of dimensions
      integer, value, intent(in):: ndims
      
      ! variable shape, global size, dimensions
      ! to create local variables optional pass adios2_null_dims 
      integer(kind=8), dimension(:), intent(in):: shape_dims
      
      ! variable start, local offset, dimensions
      ! to create local variables optional pass adios2_null_dims 
      integer(kind=8), dimension(:), intent(in):: start_dims
      
      ! variable count, local size, dimensions
      integer(kind=8), dimension(:), intent(in):: count_dims

      ! error code
      integer, intent(out) :: ierr
      
      ! .true. : constant dimensions, shape, start and count won't change 
      !          (mesh sizes, number of nodes)
      !          adios2_constant_dims = .true. use for code clarity
      ! .false. : variable dimensions, shape, start and count could change
      !           (number of particles)
      !           adios2_variable_dims = .false. use for code clarity
      logical, value, intent(in):: adios2_constant_dims
      

* available :f90:`adios2_type` parameters in :f90:`subroutine adios2_define_variable` 
   
   .. code-block:: fortran
      
      integer, parameter :: adios2_type_character = 0
      integer, parameter :: adios2_type_real = 2
      integer, parameter :: adios2_type_dp = 3
      integer, parameter :: adios2_type_complex = 4
      integer, parameter :: adios2_type_complex_dp = 5
      
      integer, parameter :: adios2_type_integer1 = 6
      integer, parameter :: adios2_type_integer2 = 7
      integer, parameter :: adios2_type_integer4 = 8
      integer, parameter :: adios2_type_integer8 = 9
      
      integer, parameter :: adios2_type_string = 10
      integer, parameter :: adios2_type_string_array = 11
  

.. tip::

   Always prefer using ``adios2_type_xxx`` parameters explicitly rather than raw numbers. 
   *e.g.* use ``adios2_type_dp`` instead of ``3``
  
  
* :f90:`subroutine adios2_inquire_variable` inquire and get a variable. See `variable%valid` to check if variable exists.
   
   .. code-block:: fortran
   
      subroutine adios2_inquire_variable(variable, io, name, ierr)
        
      ! WHERE:
      
      ! output variable handler:
      ! variable%valid = .true. points to valid found variable
      ! variable%valid = .false. variable not found
      type(adios2_variable), intent(out) :: variable
      
      ! IO in which search for variable is performed
      type(adios2_io), intent(in) :: io
      
      ! unique key name to search for variable 
      character*(*), intent(in) :: name

      ! error code
      integer, intent(out) :: ierr
      

* :f90:`subroutine adios2_set_shape` set new ``shape_dims`` for a variable if its dims are marked as varying in the define call ``adios2_define_variable``
   
   .. code-block:: fortran
   
      subroutine adios2_set_shape(variable, ndims, shape_dims, ierr)
      
      ! WHERE
      
      ! variable handler
      type(adios2_variable), intent(in) :: variable
      
      ! number of dimensions in shape_dims
      integer, intent(in) :: ndims
      
      ! new shape_dims
      integer(kind=8), dimension(:), intent(in):: shape_dims

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_set_selection` selects part of a variable through start_dims and count_dims
   
   .. code-block:: fortran
   
      subroutine adios2_set_selection(variable, ndims, start_dims, count_dims, ierr)
      
      ! WHERE
      
      ! variable handler
      type(adios2_variable), intent(in) :: variable
      
      ! number of dimensions in start_dims and count_dims
      integer, intent(in) :: ndims
      
      ! new start_dims
      integer(kind=8), dimension(:), intent(in):: start_dims
      
      ! new count_dims
      integer(kind=8), dimension(:), intent(in):: count_dims

      ! error code
      integer, intent(out) :: ierr
      

* :f90:`subroutine adios2_set_step_selection` set a list of steps by specifying the starting step and the step count
   
   .. code-block:: fortran
   
      subroutine adios2_set_step_selection(variable, step_start, step_count, ierr)
      
      ! WHERE
      
      ! variable handler
      type(adios2_variable), intent(in) :: variable
      
      ! new step_start 
      integer(kind=8), intent(in):: step_start
      
      ! new step_count (or number of steps to read from step_start)
      integer(kind=8), intent(in):: step_count

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_variable_max` get the maximum value in the variable array
  
   .. code-block:: fortran

      subroutine adios2_variable_max(maximum, variable, ierr)

      ! WHERE

      ! scalar variable that will contain the maximum value
      Generic Fortran types, intent(out) :: maximum

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_variable_min` get the minimum value in the variable array
  
   .. code-block:: fortran

      subroutine adios2_variable_min(minimum, variable, ierr)

      ! WHERE

      ! scalar variable that will contain the minimum value
      Generic Fortran types, intent(out) :: minimum

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_add_operation` add an operation to a variable

   .. code-block:: fortran

      subroutine adios2_add_operation(operation_index, variable, op, key, value, ierr)

      ! WHERE

      ! reference to the operator handle that will be created
      integer, intent(out):: operation_index

      ! variable handler
      type(adios2_variable), intent(in):: variable

      ! Operator handler
      type(adios2_operator), intent(in):: op

      ! Operator key
      character*(*), intent(in):: key

      ! Operator value
      character*(*), intent(in):: value

      ! error code
      integer, intent(out):: ierr


* :f90:`subroutine adios2_set_operation_parameter` set a parameter for a operator. Replaces value if parameter already exists.

   .. code-block:: fortran

      subroutine adios2_set_operation_parameter(variable, operation_index, key, value, ierr)

      ! WHERE

      ! variable handler
      type(adios2_variable), intent(in):: variable

      ! Operation index handler
      integer, intent(in):: operation_index

      ! parameter key
      character*(*), intent(in):: key

      ! parameter value
      character*(*), intent(in):: value

      ! error code
      integer, intent(out):: ierr


* :f90:`subroutine adios2_variable_name` retrieve variable name

   .. code-block:: fortran

      subroutine adios2_variable_name(name, variable, ierr)

      ! WHERE

      ! variable name
      character(len=:), allocatable, intent(out) :: name

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_variable_type` retrieve variable datatype

   .. code-block:: fortran

      subroutine adios2_variable_type(type, variable, ierr)

      ! WHERE

      ! variable type
      integer, intent(out) :: type

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_variable_ndims` retrieve number of dimensions for a variable

   .. code-block:: fortran

      subroutine adios2_variable_ndims(ndims, variable, ierr)

      ! WHERE

      ! No. of dimensions
      integer, intent(out) :: ndims

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_variable_shape` retrieve the shape of a variable

   .. code-block:: fortran

      subroutine adios2_variable_shape(shape_dims, ndims, variable, ierr)

      ! WHERE

      ! array that contains the shape
      integer(kind=8), dimension(:), allocatable, intent(out) :: shape_dims

      ! no. of dimensions
      integer, intent(out) :: ndims

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_variable_steps` retrieve the number of available steps

   .. code-block:: fortran

      subroutine adios2_variable_steps(steps, variable, ierr)

      ! WHERE

      ! no. of steps
      integer(kind=8), intent(out) :: steps

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_set_block_selection` Read mode only. Required for reading local variables. For global arrays it will set the appropriate Start and Count selection for the global array coordinates.

   .. code-block:: fortran

      subroutine adios2_set_block_selection(variable, block_id, ierr)

      ! WHERE

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! variable block index defined at write time. Blocks can be inspected with `bpls -D variableName`
      integer(kind=8), intent(in) :: block_id

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_set_memory_selection` set the local start (offset) point to the memory pointer passed at adios2_put and the memory local dimensions (count). Used for non-contiguous memory writes and reads (e.g. multidimensional ghost-cells). Currently Get only works for formats based on BP.

   .. code-block:: fortran

      subroutine adios2_set_memory_selection(variable, ndims, memory_start_dims, memory_count_dims, ierr)

      ! WHERE

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! no. of dimensions of the variable
      integer, intent(in) :: ndims

      ! memory start offsets
      integer(kind=8), dimension(:), intent(in) :: memory_start_dims

      ! no. of elements in each dimension
      integer(kind=8), dimension(:), intent(in) :: memory_count_dims

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_set_step_selection` set a step selection modifying current step_start, step_count. step_count is the number of steps from step_start

   .. code-block:: fortran

      subroutine adios2_set_step_selection(variable, step_start, step_count, ierr)

      ! WHERE

      ! variable handler
      type(adios2_variable), intent(in) :: variable

      ! starting step
      integer(kind=8), intent(in) :: step_start

      ! no. of steps from start
      integer(kind=8), intent(in) :: step_count

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_remove_operations` remove all current operations associated with the variable. Provides the posibility to apply operators on a block basis.

   .. code-block:: fortran

      subroutine adios2_remove_operations(variable, ierr)

      ! WHERE

      ! variable handler
      type(adios2_variable), intent(in):: variable

      ! error code
      integer, intent(out):: ierr


:ref:`Engine` subroutines
-------------------------

* :f90:`subroutine adios2_open` opens an engine to execute IO tasks 
   
   .. code-block:: fortran
   
      ! MPI version: duplicates communicator from adios2_init
      ! Non-MPI serial version  
      subroutine adios2_open(engine, io, name, adios2_mode, ierr)
      
      ! MPI version only to pass a communicator other than the one from adios_init 
      subroutine adios2_open(engine, io, name, adios2_mode, comm, ierr)
      
      ! WHERE:
      
      ! handler to newly opened adios2 engine
      type(adios2_engine), intent(out) :: engine
      
      ! IO that spawns an engine based on its configuration
      type(adios2_io), intent(in) :: io
      
      ! unique engine identifier within io, file name for default BPFile engine 
      character*(*), intent(in) :: name
      
      ! Optional MPI communicator, only in MPI library
      integer, intent(in) :: comm

      ! error code
      integer, intent(out) :: ierr
      
      ! open mode parameter: 
      !                      adios2_mode_write,
      !                      adios2_mode_append,
      !                      adios2_mode_read,  
      integer, intent(in):: adios2_mode


* :f90:`subroutine adios2_begin_step` begin a new step or progress to the next step. Starts from 0
   
   .. code-block:: fortran
   
      subroutine adios2_begin_step(engine, adios2_step_mode, timeout_seconds, status, ierr)
      ! Default Timeout = -1.    (block until step available)
      subroutine adios2_begin_step(engine, adios2_step_mode, ierr)
      ! Default step_mode for read and write
      subroutine adios2_begin_step(engine, ierr)
      
      ! WHERE
      
      ! engine handler
      type(adios2_engine), intent(in) :: engine
      
      ! step_mode parameter:
      !     adios2_step_mode_read (read mode default)
      !     adios2_step_mode_append (write mode default)
      integer, intent(in):: adios2_step_mode
      
      ! optional 
      ! engine timeout (if supported), in seconds
      real, intent(in):: timeout_seconds

      ! status of the stream from adios2_step_status_* parameters
      integer, intent(out):: status

      ! error code
      integer, intent(out) :: ierr
   
      
* :f90:`subroutine adios2_current_step` extracts current step number
   
   .. code-block:: fortran
   
      subroutine adios2_current_step(current_step, engine, ierr)
      
      ! WHERE:

      ! engine handler  
      type(adios2_engine), intent(in) :: engine
      
      ! populated with current_step value
      integer(kind=8), intent(out) :: current_step

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_steps` inspect total number of available steps, use for file engines in read mode only
   
   .. code-block:: fortran
   
      subroutine adios2_steps(steps, engine, ierr)
      
      ! WHERE:

      ! engine handler  
      type(adios2_engine), intent(in) :: engine
      
      ! populated with steps value
      integer(kind=8), intent(out) :: steps

      ! error code
      integer, intent(out) :: ierr

      
* :f90:`subroutine adios2_end_step` end current step and execute transport IO (flush or read). 
   
   .. code-block:: fortran
   
      subroutine adios2_end_step(engine, ierr)
      
      ! WHERE:

      ! engine handler  
      type(adios2_engine), intent(in) :: engine

      ! error code
      integer, intent(out) :: ierr
   

* :f90:`subroutine adios2_put` put variable data and metadata into adios2 for IO operations. Default is deferred mode. For optional sync mode, see :ref:`Put: modes and memory contracts`. Variable and data types must match.
   
   .. code-block:: fortran
   
      subroutine adios2_put(engine, variable, data, adios2_mode, ierr)
      
      ! Default adios2_mode_deferred
      subroutine adios2_put(engine, variable, data, ierr)
      
      ! WHERE:
      
      ! engine handler  
      type(adios2_engine), intent(in) :: engine
      
      ! variable handler containing metadata information  
      type(adios2_variable), intent(in) :: variable
      
      ! Fortran bindings supports data types from adios2_type in variables, 
      ! up to 6 dimensions 
      ! Generic Fortran type from adios2_type
      Generic Fortran types, intent(in):: data 
      Generic Fortran types, dimension(:), intent(in):: data
      Generic Fortran types, dimension(:,:), intent(in):: data
      Generic Fortran types, dimension(:,:,:), intent(in):: data
      Generic Fortran types, dimension(:,:,:,:), intent(in):: data
      Generic Fortran types, dimension(:,:,:,:,:), intent(in):: data
      Generic Fortran types, dimension(:,:,:,:,:,:), intent(in):: data
      
      ! mode:
      ! adios2_mode_deferred: won't execute until adios2_end_step, adios2_perform_puts or adios2_close
      ! adios2_mode_sync: special case, put data immediately, can be reused after this call
      integer, intent(in):: adios2_mode

      ! error code
      integer, intent(out) :: ierr
      
      
* :f90:`subroutine adios2_perform_puts` execute deferred calls to ``adios2_put``
      
   .. code-block:: fortran
   
      subroutine adios2_perform_puts(engine, ierr)
      
      ! WHERE:
      
      ! engine handler  
      type(adios2_engine), intent(in) :: engine

      ! error code
      integer, intent(out) :: ierr
      
      
* :f90:`subroutine adios2_get` get variable data into ADIOS2 for IO operations. Default is deferred mode. For optional sync mode, see :ref:`Get: modes and memory contracts`. Variable and data types must match, variable can be obtained from ``adios2_inquire_variable``. Memory for data must be pre-allocated.

   .. code-block:: fortran
   
      subroutine adios2_get(engine, variable, data, adios2_mode, ierr)
      
      ! Default adios2_mode_deferred
      subroutine adios2_get(engine, variable, data, ierr)
      
      ! WHERE:
      
      ! engine handler  
      type(adios2_engine), intent(in) :: engine
      
      ! variable handler containing metadata information  
      type(adios2_variable), intent(in) :: variable
      
      ! Fortran bindings supports data types from adios2_type in variables, 
      ! up to 6 dimensions. Must be pre-allocated 
      ! Generic Fortran type from adios2_type
      Generic Fortran types, intent(out):: data 
      Generic Fortran types, dimension(:), intent(out):: data
      Generic Fortran types, dimension(:,:), intent(out):: data
      Generic Fortran types, dimension(:,:,:), intent(out):: data
      Generic Fortran types, dimension(:,:,:,:), intent(out):: data
      Generic Fortran types, dimension(:,:,:,:,:), intent(out):: data
      Generic Fortran types, dimension(:,:,:,:,:,:), intent(out):: data
      
      ! mode:
      ! adios2_mode_deferred: won't execute until adios2_end_step, adios2_perform_gets or adios2_close
      ! adios2_mode_sync: special case, get data immediately, can be reused after this call
      integer, intent(in):: adios2_mode

      ! error code
      integer, intent(out) :: ierr
      
      
* :f90:`subroutine adios2_perform_gets` execute deferred calls to ``adios2_get``
      
   .. code-block:: fortran
   
      subroutine adios2_perform_gets(engine, ierr)
      
      ! WHERE:
      
      ! engine handler  
      type(adios2_engine), intent(in) :: engine

      ! error code
      integer, intent(out) :: ierr
      
      
* :f90:`subroutine adios2_close` close engine. May re-open.
      
   .. code-block:: fortran
   
      subroutine adios2_close(engine, ierr)
      
      ! WHERE:
      
      ! engine handler  
      type(adios2_engine), intent(in) :: engine

      ! error code
      integer, intent(out) :: ierr
      

* :f90:`subroutine adios2_io_engine_type` get current engine type

   .. code-block:: fortran

      subroutine adios2_io_engine_type(type, io, ierr)

      ! WHERE

      ! engine type (BP, SST, SSC, HDF5, DataMan)
      character(len=:), allocatable, intent(out) :: type

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_lock_writer_definitions` promise that no more definitions or changes to defined variables will occur. Useful information if called before the first EndStep() of an output Engine, as it will know that the definitions are complete and constant for the entire lifetime of the output and may optimize metadata handling.

   .. code-block:: fortran

      subroutine adios2_lock_writer_definitions(engine, ierr)

      ! WHERE

      ! adios2 engine handler
        type(adios2_engine), intent(in) :: engine

      ! error code
        integer, intent(out) :: ierr


* :f90:`subroutine adios2_lock_reader_selections` promise that the reader data selections of are fixed and will not change in future timesteps. This information, provided before the end_step() representing a fixed read pattern, may be utilized by the input Engine to optimize data flow.

   .. code-block:: fortran

      subroutine adios2_lock_reader_selections(engine, ierr)

      ! WHERE

      ! adios2 engine handler
        type(adios2_engine), intent(in) :: engine

      ! error code
        integer, intent(out) :: ierr


:ref:`Operator` subroutines
---------------------------

* :f90:`subroutine adios2_operator_type` get current Operator type

   .. code-block:: fortran

      subroutine adios2_operator_type(type, op, ierr)

      ! WHERE

      ! Operator type name. See list of supported operator types.
      character(len=:), allocatable, intent(out) :: type
      
      ! Operator handler
      type(adios2_operator), intent(in) :: op

      ! error code
      integer, intent(out) :: ierr


:ref:`Attribute` subroutines
----------------------------

* :f90:`subroutine adios2_define_attribute` define/create a new user attribute
   
   .. code-block:: fortran

      ! Single value attributes
      subroutine adios2_define_attribute(attribute, io, attribute_name, data, ierr)
                                         
      ! 1D array attributes
      subroutine adios2_define_attribute(attribute, io, attribute_name, data, elements, ierr)
         
      ! WHERE:
      
      ! handler to newly defined attribute
      type(adios2_attribute), intent(out):: attribute 
      
      ! IO component owning the attribute
      type(adios2_io), intent(in):: io
      
      ! unique attribute identifier within io
      character*(*), intent(in):: attribute_name
      
      ! overloaded subroutine allows for multiple attribute data types
      ! they can be single values or 1D arrays
      Generic Fortran types, intent(in):: data 
      Generic Fortran types, dimension(:), intent(in):: data
                                        
      ! number of elements if passing a 1D array in data argument
      integer, intent(in):: elements

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_inquire_attribute` inquire for existing attribute by its unique name
   
   .. code-block:: fortran
   
      subroutine adios2_inquire_attribute(attribute, io, name, ierr)
        
      ! WHERE:
      
      ! output attribute handler:
      ! attribute%valid = .true. points to valid found attribute
      ! attribute%valid = .false. attribute not found
      type(adios2_attribute), intent(out) :: attribute
      
      ! IO in which search for attribute is performed
      type(adios2_io), intent(in) :: io
      
      ! unique key name to search for attribute 
      character*(*), intent(in) :: name

      ! error code
      integer, intent(out) :: ierr

..  caution::

   Use the ``adios2_remove_*`` subroutines with extreme CAUTION.
   They create outdated dangling information in the ``adios2_type`` handlers.
   If you don't need them, don't use them. 


* :f90:`subroutine adios2_attribute_data` retrieve attribute data

   .. code-block:: fortran

      subroutine adios2_attribute_data(data, attribute, ierr)

      ! WHERE

      ! data handler
      character*(*), intent(out):: data
      real, intent(out):: data
      real(kind=8), intent(out):: data
      integer(kind=1), intent(out):: data
      integer(kind=2), intent(out):: data
      integer(kind=4), intent(out):: data
      integer(kind=8), intent(out):: data
      character*(*), dimension(:), intent(out):: data
      real, dimension(:), intent(out):: data
      real(kind=8), dimension(:), intent(out):: data
      integer(kind=2), dimension(:), intent(out):: data
      integer(kind=4), dimension(:), intent(out):: data
      integer(kind=8), dimension(:), intent(out):: data


      ! attribute
      type(adios2_attribute), intent(in):: attribute

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_attribute_name` inspect attribute name

   .. code-block:: fortran

      subroutine adios2_attribute_name(name, attribute, ierr)

      ! WHERE

      ! name to be output
      character(len=:), allocatable, intent(out) :: name

      ! attribute handler
      type(adios2_attribute), intent(in) :: attribute

      ! error code
      integer, intent(out) :: ierr


* :f90:`subroutine adios2_inquire_variable_attribute` retrieve a handler to a previously defined attribute associated to a variable

   .. code-block:: fortran

      subroutine adios2_inquire_variable_attribute(attribute, io, attribute_name, variable_name, separator, ierr)

      ! WHERE

      ! attribute handler
      type(adios2_attribute), intent(out) :: attribute

      ! IO handler
      type(adios2_io), intent(in) :: io

      ! attribute name
      character*(*), intent(in) :: attribute_name

      ! variable name
      character*(*), intent(in) :: variable_name

      ! hierarchy separator (e.g. “/” in variable/attribute )
      character*(*), intent(in) :: separator

      ! error code
      integer, intent(out) :: ierr
