"""License:
  Distributed under the OSI-approved Apache License, Version 2.0.  See
  accompanying file Copyright.txt for details.
"""


class Variable:
    """High level representation of the Attribute class in the adios2.bindings"""

    def __init__(self, implementation):
        self.impl = implementation

    @property
    def impl(self):
        """Bindings implementation of the class"""
        return self._impl

    @impl.setter
    def impl(self, implementation):
        self._impl = implementation

    def __eq__(self, other):
        if isinstance(other, Variable):
            return self.name() == other.name()
        return False

    def block_id(self):
        """
        BlockID of this variable.

        Returns:
            int: BlockID of this variable.
        """
        return self.impl.BlockID()

    def count(self):
        """
        Current selected count for this variable.

        Returns:
            int: Current selected count.
        """
        return self.impl.Count()

    def selection_size(self):
        """
        Current selection size selected for this variable.

        Returns:
            int: Current selection size selected.
        """
        return self.impl.SelectionSize()

    def set_block_selection(self, block_id):
        """
        Set BlockID for this variable.

        Args:
            block_id (int): Selected BlockID.
        """
        self.impl.SetBlockSelection(block_id)

    def set_selection(self, selection):
        """
        Set selection for this variable.

        Args:
            selection (list): list of the shape [[start], [count]], note that start and
            count can contain more than one element.
        """
        self.impl.SetSelection(selection)

    def set_shape(self, shape):
        """
        Set Shape (dimensions) for this variable

        Args:
            shape (list): desired shape (dimensions).
        """
        self.impl.SetShape(shape)

    def set_step_selection(self, step_selection):
        """
        Set current step selection (For RRA or ReadRandomAccess)

        Args:
            step_selection (list): On the form of [start, count].
        """
        self.impl.SetStepSelection(step_selection)

    def shape(self, step=None):
        """
        Get the shape assigned to the given step for this variable.

        Args:
            step (int): Desired step. Only in ReadRandomAccess mode

        Returns:
            list: shape of the specified step in the form of [start, count].
        """
        if step is None:
            return self.impl.Shape()
        return self.impl.Shape(step)

    def shape_id(self):
        """
        Get the ShapeID assigned to this variable.

        Returns:
            int: ShapeID assigned to this variable.
        """
        return self.impl.ShapeID()

    def type(self):
        """
        Type of the Variable

        Returns:
            str: Type of the Variable.
        """
        return self.impl.Type()

    def single_value(self):
        """
        Check if this variable is a single value.

        Returns:
            bool: true if this variable is a single value.
        """
        return bool(self.impl.SingleValue())

    def sizeof(self):
        """
        Size in bytes of the contents of the variable.

        Returns:
            int: size in bytes of the contents.
        """
        return self.impl.Sizeof()

    def start(self):
        """
        The current selected start of the variable.

        Returns:
            int: The current selected start of the variable.
        """
        return self.impl.Start()

    def steps(self):
        """
        The number of steps of the variable. This is always 1 in a stream.
        In ReadRandomAccess mode, this function returns the total number
        of steps available, which can be used when selecting steps for read.

        Returns:
            int: The avaialble steps of the variable.
        """
        return self.impl.Steps()

    def steps_start(self):
        """
        The avaliable start step, this is needed variables can start
        at any time step. This is for ReadRandomAccess.

        Returns:
            int: the starting step of for this Variable.
        """
        return self.impl.StepsStart()

    def name(self):
        """
        Name of the Variable

        Returns:
            str: Name of the Variable.
        """
        return self.impl.Name()

    def add_operation_string(self, name, params={}):
        """
        Add an operation (operator) as a string

        Args:
            name (str): name of the operation.
            params (dict): parameters as a form of a dict for the operation.
        """
        return self.impl.AddOperation(name, params)

    def add_operation(self, operation, params={}):
        """
        Add an operation (operator).

        Args:
            name (Operator): name of the operation.
            params (dict): parameters as a form of a dict for the operation.
        """
        return self.impl.AddOperation(operation.impl, params)

    def operations(self):
        """
        Current operations (operators) assigned to this Variable.

        Returns:
            list(Operators): operators assigned.
        """
        return self.impl.Operations()

    def remove_operations(self):
        """
        Remove operations (operators) assigned to this Variable.
        """
        self.impl.RemoveOperations()
