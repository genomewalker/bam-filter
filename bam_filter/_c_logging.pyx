# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
# -*- coding: utf-8 -*-
# Thin Cython wrapper exposing bf_set_verbosity/get_verbosity and nogil logging for Python

cimport cython
from libc.stdlib cimport NULL
from bam_filter._c_logging cimport (
    bf_nogil_log_level,
    bf_set_verbosity,
    bf_get_verbosity,
    bf_should_log,
)

@cython.boundscheck(False)
@cython.wraparound(False)
def set_verbosity(int level):
    bf_set_verbosity(level)

def get_verbosity():
    return bf_get_verbosity()

def should_log(int level):
    return bf_should_log(level) != 0

@cython.boundscheck(False)
@cython.wraparound(False)
def nogil_log(tag: str, msg: str, int level=0):
    """Call nogil C logger with optional level gating."""
    if not bf_should_log(level):
        return
    cdef bytes msg_b = msg.encode("utf-8")
    if tag is None or tag == "":
        bf_nogil_log_level(level, NULL, msg_b)
    else:
        cdef bytes tag_b = tag.encode("utf-8")
        bf_nogil_log_level(level, tag_b, msg_b)
