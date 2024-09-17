"""License:
  Distributed under the OSI-approved Apache License, Version 2.0.  See
  accompanying file Copyright.txt for details.
"""


class Operator:
    """High level representation of the Attribute class in the adios2.bindings"""

    def __init__(self, implementation, name):
        self.impl = implementation
        self._name = name

    @property
    def impl(self):
        """Bindings implementation of the class"""
        return self._impl

    @impl.setter
    def impl(self, implementation):
        self._impl = implementation

    def __eq__(self, other):
        if isinstance(other, Operator):
            return self._name == other._name
        return False

    def get_parameters(self):
        """
        Get parameters associated to this Operator

        Returns:
            dict: parameters
        """
        return self.impl.Parameters()

    def set_parameter(self, key, value):
        """
        Set parameter associated to this Operator

        Args:
            str: key
            str: value
        """
        self.impl.SetParameter(key, value)
