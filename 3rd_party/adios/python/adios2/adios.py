"""License:
  Distributed under the OSI-approved Apache License, Version 2.0.  See
  accompanying file Copyright.txt for details.
"""

from adios2 import bindings
from adios2.io import IO
from adios2.operator import Operator


class Adios:
    """High level representation of the ADIOS class in the adios2.bindings"""

    def __init__(self, config_file=None, comm=None):
        if comm and not bindings.is_built_with_mpi:
            raise RuntimeError("Cannot use MPI since ADIOS2 was built without MPI support")

        if config_file:
            if comm:
                self.impl = bindings.ADIOS(config_file, comm)
            else:
                self.impl = bindings.ADIOS(config_file)
        else:
            if comm:
                self.impl = bindings.ADIOS(comm)
            else:
                self.impl = bindings.ADIOS()

    @property
    def impl(self):
        """Bindings implementation of the class"""
        return self._impl

    @impl.setter
    def impl(self, implementation):
        self._impl = implementation

    def declare_io(self, name):
        """
        Declare IO instance

        Args:
            name (str): IO instance name
        """
        return IO(self.impl.DeclareIO(name), name, self)

    def at_io(self, name):
        """
        Inquire IO instance

        Args:
            name (str): IO instance name
        Returns:
            IO: IO instance
        """
        io = None
        io_instance = self.impl.AtIO(name)
        if io_instance:
            io = IO(io_instance, name, self)
        return io

    def remove_io(self, name):
        """
        Remove IO instance

        Args:
            name (str): IO instance name
        """
        self.impl.RemoveIO(name)

    def remove_all_ios(self):
        """Remove all IO instance"""
        self.impl.RemoveAllIOs()

    def define_operator(self, name, kind, parameters={}):
        """
        Add an operation (operator).

        Args:
            name (str): name of the operator.
            kind (str): name type of the operation.
            params (dict): parameters as a form of a dict for the operation.
        """
        return Operator(self.impl.DefineOperator(name, kind, parameters), name)

    def inquire_operator(self, name):
        """
        Add an operation (operator).

        Args:
            name (str): name of the operator.
        Return:
            Operator: requested operator
        """
        operator = None
        operator_instance = self.impl.InquireOperator(name)
        if operator_instance:
            operator = Operator(self.impl.InquireOperator(name), name)

        return operator

    def flush_all(self):
        """
        Flush all IO instances
        """
        self.impl.FlushAll()
