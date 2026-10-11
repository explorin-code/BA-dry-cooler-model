"""
quiet.py
=========
Two ways to silence the model while a lot of solves run.

  silence_all()   benchmarks: no log output, no warnings, no stdout -- output would distort the timings
  quiet_solvers() studies of thousands of points: the solvers' per-point log lines and the range warnings go away
                  (the flags are in the result rows instead); progress messages of the study itself still show
"""

import io
import logging
import warnings
from contextlib import contextmanager, redirect_stdout


@contextmanager
def silence_all():
    logging.disable(logging.CRITICAL)
    try:
        with warnings.catch_warnings(), redirect_stdout(io.StringIO()):
            warnings.simplefilter("ignore")
            yield
    finally:
        logging.disable(logging.NOTSET)


@contextmanager
def quiet_solvers():
    src_logger = logging.getLogger('src')
    old_level = src_logger.level
    src_logger.setLevel(logging.ERROR)
    try:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='Wang correlation')
            warnings.filterwarnings('ignore', message='NTU element field')
            yield
    finally:
        src_logger.setLevel(old_level)
