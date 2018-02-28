
import os
from collections import OrderedDict as ODict

def parameter_lookup_dict(dct):
    """ Make a dict to lookup where a given parameter is stored. """
    lt = ODict()
    for group in dct:
        for key in dct[group]:
            if key in lt:
                lt[key].append(group)
            else:
                lt[key] = [group]
    # Finally add the key itself at the end
    for key in lt:
        lt[key].append(key)
    return lt

from functools import reduce
import operator

def getFromDict(dataDict, mapList):
    return reduce(operator.getitem, mapList, dataDict)

def setInDict(dataDict, mapList, value):
    getFromDict(dataDict, mapList[:-1])[mapList[-1]] = value

class ParamSet:

    def __init__(self, parser = None):
        if parser is None:
            raise RuntimeError("Need to initialize with a parameter file parser!")
        self.parser = parser
        self.dct = parser.get_param_dict()
        self.make_parameter_lookup_dict()

    def make_parameter_lookup_dict(self):
        """ Construct a dict to map names of parameters to a list of
        keys to find them in the dict """
        self.parameter_lookup_dict = parameter_lookup_dict(self.dct)

    def set_param(self, param, value):
        """ Set the parameter param to its new value. """
        # Parsers hold values inside of lists to handle
        # multiple values per param.
        if isinstance(value, str):
            value = [value]
        # make a string out of anything else if it doesn't have a length
        try:
            len(value)
        except TypeError:
            value = ["".format(value)]
        keyList = self.parameter_lookup_dict[param]
        setInDict(self.dct, keyList, value)

    def get_param(self, param):
        """ Get the parameter param. """
        keyList = self.parameter_lookup_dict[param]
        value = getFromDict(self.dct, keyList)
        # Parsers hold values inside of lists to handle
        # multiple values per param.
        if len(value) == 1:
            return value[0]
        else:
            return value
