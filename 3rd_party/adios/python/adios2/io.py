"""License:
  Distributed under the OSI-approved Apache License, Version 2.0.  See
  accompanying file Copyright.txt for details.
"""

import numpy as np
from adios2.attribute import Attribute
from adios2.variable import Variable
from adios2.engine import Engine


class IO:
    """High level representation of the IO class in the adios2.bindings"""

    def __init__(self, impl, name, adiosobj):
        self.impl = impl
        self._name = name
        self._adios = adiosobj

    @property
    def impl(self):
        """Bindings implementation of the class"""
        return self._impl

    @impl.setter
    def impl(self, implementation):
        self._impl = implementation

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.flush_all()

    def __eq__(self, other):
        if isinstance(other, IO):
            return self._name == other._name
        return False

    def adios(self):
        """Adios instance associated to this IO"""
        return self._adios

    def define_attribute(
        self,
        name,
        content=None,
        variable_name="",
        separator="/",
    ):
        """
        Define an Attribute

        Parameters
            name
                attribute name

            content
                attribute numpy array data

            variable_name
                if attribute is associated with a variable

            separator
                concatenation string between variable_name and attribute
                e.g. variable_name + separator + name ("var/attr")
                Not used if variable_name is empty
        """

        if isinstance(content, (np.ndarray, str, list)):
            return Attribute(self.impl, name, content, variable_name, separator)

        if not hasattr(content, "__len__") and isinstance(content, (int, float, complex)):
            return Attribute(self.impl, name, content, variable_name, separator)

        # Anything else, try to pass as a numpy array (and get an error here if that fails)
        return Attribute(self.impl, name, np.asarray(content), variable_name, separator)

    def inquire_attribute(self, name, variable_name="", separator="/"):
        """
        Inquire an Attribute

        Parameters
            name
                attribute name

            variable_name
                if attribute is associated with a variable

            separator
                concatenation string between variable_name and attribute
                e.g. variable_name + separator + name ("var/attr")
                Not used if variable_name is empty
        """
        attr = None
        attrimpl = self.impl.InquireAttribute(name, variable_name, separator)
        if attrimpl:
            attr = Attribute.__new__(Attribute)
            attr.impl = attrimpl
        return attr

    def available_attributes(self, varname="", separator="/"):
        """
        Returns a 2-level dictionary with attribute information.
        Read mode only.

        Parameters
            variable_name
                If varname is set, attributes assigned to that variable are returned.
                The keys returned are attribute names with the prefix of varname + separator
                removed.

            separator
                concatenation string between variable_name and attribute
                e.g. varname + separator + name ("var/attr")
                Not used if varname is empty

        Returns
            attributes dictionary
                key
                    attribute name
                value
                    attribute information dictionary
        """
        return self.impl.AvailableAttributes(varname, separator)

    def remove_attribute(self, name):
        """
        Remove an Attribute

        Parameters
            name
                attribute name
        """
        self.impl.RemoveAttribute(name)

    def remove_all_attributes(self):
        """
        Remove all define attributes
        """
        self.impl.RemoveAllAttributes()

    def define_variable(
        self,
        name,
        content=None,
        shape=[],
        start=[],
        count=[],
        is_constant_dims=False,
    ):
        """
        writes a variable

        Parameters
            name
                variable name

            content
                variable data values

            shape
                variable global MPI dimensions.

            start
                variable offset for current MPI rank.

            count
                variable dimension for current MPI rank.

            isConstantDims
                Whether dimensions are constant
        """
        var_impl = None
        if isinstance(content, list) and len(content) > 0:
            var_impl = self.impl.DefineVariable(
                name, content[0], shape, start, count, is_constant_dims
            )

        elif content is None:
            var_impl = self.impl.DefineVariable(name, "", shape, start, count, is_constant_dims)

        else:
            var_impl = self.impl.DefineVariable(
                name, content, shape, start, count, is_constant_dims
            )

        return Variable(var_impl)

    def inquire_variable(self, name):
        """
        Inquire a variable

        Parameters
            name
                variable name
        Returns
            The variable if it is defined, otherwise None
        """
        var = None
        var_impl = self.impl.InquireVariable(name)
        if var_impl:
            var = Variable.__new__(Variable)
            var.impl = var_impl
        return var

    def available_variables(self):
        """

        Returns a 2-level dictionary with variable information.
        Read mode only.

        Parameters
            keys
               list of variable information keys to be extracted (case insensitive)
               keys=['AvailableStepsCount','Type','Max','Min','SingleValue','Shape']
               keys=['Name'] returns only the variable names as 1st-level keys
               leave empty to return all possible keys

        Returns
            variables dictionary
                key
                    variable name
                value
                    variable information dictionary
        """
        return self.impl.AvailableVariables()

    def remove_variable(self, name):
        """
        Remove a variable

        Parameters
            name
                Variable name
        """
        self.impl.RemoveVariable(name)

    def remove_all_variables(self):
        """
        Remove all variables in the IO instance
        """
        self.impl.RemoveAllVariables()

    def open(self, name, mode, comm=None):
        """
        Open an engine

        Parameters
            name
                Engine name

            mode
                engine mode

            comm
                MPI communicator, optional
        """

        if comm:
            return Engine(self.impl.Open(name, mode, comm))

        return Engine(self.impl.Open(name, mode))

    def set_engine(self, name):
        """
        Set engine for this IO instance

        Args:
            name (str): name of engine
        """
        self.impl.SetEngine(name)

    def engine_type(self):
        """Return engine type"""
        return self.impl.EngineType()

    def add_transport(self, kind, parameters={}):
        """
        Adds a transport and its parameters to current IO. Must be
        supported by current engine type.

        Parameters
            kind
                must be a supported transport type for current engine.

            parameters
                acceptable parameters for a particular transport
                CAN'T use the keywords "Transport" or "transport" in key

        Returns
            transport_index
                handler to added transport
        """
        return self.impl.AddTransport(kind, parameters)

    def parameters(self):
        """
        Return parameter associated to this io instance

        Return:
            dict: str->str parameters
        """
        return self.impl.Parameters()

    def set_parameter(self, key, value):
        """
        Set a parameter for this IO instance

        Args:
            key (str):
            value (str):
        """
        self.impl.SetParameter(key, value)

    def set_parameters(self, parameters):
        """
        Set parameters for this IO instance

        Args:
            parameters (dict):
        """
        self.impl.SetParameters(parameters)

    def flush_all(self):
        """Flush all engines attached to this IO instance"""
        self.impl.FlushAll()
