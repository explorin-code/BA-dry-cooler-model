"""
param_loader.py
================
Builds a dataclass instance from parameters.py constants via an explicit
{field_name: PARAMETERS_CONSTANT_NAME} mapping; a mapping value can also be
(CONSTANT_NAME, MM) / (CONSTANT_NAME, PERCENT) to convert an input given in mm /
% to SI metres / a fraction. Fields left out of the
mapping fall back to the dataclass's own default; if it has none,
dataclasses raises for us. extra_kwargs (e.g. an explicitly-passed geo)
are merged in as-is, taking priority over the parameters.py mapping.
"""

import src.parameters as parameters

MM = 1000.0          # divisor mm -> m (division keeps e.g. 9 mm -> exactly 0.009 m)
PERCENT = 100.0      # divisor % -> fraction


def from_parameters(cls, mapping: dict, **extra_kwargs):
    kwargs = {}
    for field, source in mapping.items():
        const_name, divisor = source if isinstance(source, tuple) else (source, None)
        value = getattr(parameters, const_name)
        if value is not None:
            kwargs[field] = value / divisor if divisor else value
    kwargs.update({k: v for k, v in extra_kwargs.items() if v is not None})
    return cls(**kwargs)