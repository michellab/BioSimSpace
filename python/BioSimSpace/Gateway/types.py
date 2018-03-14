"""
@package biosimspace
@author  Lester Hedges
@brief   A collection of variable types.
"""

class Boolean():
    def __init__(self, value):

    self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value):
        if type(value) is bool:
            self._value = value
        else:
            raise TypeError("Argument is not of type 'bool'")
