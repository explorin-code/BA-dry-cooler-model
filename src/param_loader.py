"""
param_loader.py
================
Builds a dataclass instance from parameters.py constants via an explicit
{field_name: PARAMETERS_CONSTANT_NAME} mapping. Fields left out of the
mapping fall back to the dataclass's own default; if it has none,
dataclasses raises for us. extra_kwargs (e.g. an explicitly-passed geo)
are merged in as-is, taking priority over the parameters.py mapping.
"""

import src.parameters as parameters


def from_parameters(cls, mapping: dict, **extra_kwargs):
    kwargs = {}
    for field, const_name in mapping.items():
        value = getattr(parameters, const_name)
        if value is not None:
            kwargs[field] = value
    kwargs.update({k: v for k, v in extra_kwargs.items() if v is not None})
    return cls(**kwargs)