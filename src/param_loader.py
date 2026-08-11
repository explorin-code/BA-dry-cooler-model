"""
param_loader.py
================
Builds a dataclass instance from parameters.py constants via an explicit
{field_name: PARAMETERS_CONSTANT_NAME} mapping. Fields left out of the
mapping fall back to the dataclass's own default; if it has none,
dataclasses raises for us.
"""

import src.parameters as parameters


def from_parameters(cls, mapping: dict):
    kwargs = {}
    for field, const_name in mapping.items():
        value = getattr(parameters, const_name)
        if value is not None:
            kwargs[field] = value
    return cls(**kwargs)