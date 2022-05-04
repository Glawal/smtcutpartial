#! /usr/bin/env python3

# $+HEADER$
#
# Copyright 2019-2021 Christoph Lueders
#
# This file is part of the SMTcut project: <http://wrogn.com/smtcut>
#
# SMTcut is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SMTcut is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SMTcut.  If not, see <http://www.gnu.org/licenses/>.
#
# $-HEADER$

# is imported by: preproc, main, smtcut, tropicalize

try:
    if False:
        from time import perf_counter as mytime
    else:
        from time import process_time as mytime
except ImportError:
    from time import clock as mytime 


def combos(bb):
    """
    Theoretical number of combinations from bag bb.
    """
    from functools import reduce
    from operator import mul
    return reduce(mul, (len(b) for b in bb), 1)


def nice_combos(bb):
    """
    Theoretical number of combinations from bag bb.
    Print as int for small numbers, as power of 10 for larger ones.
    """
    from math import log
    prod = combos(bb)
    if prod < 10000:
        return str(prod)
    return f"10^{log(prod,10):.1f}"


def get_isodatetime():
    """
    Print date in ISO format like this: "2021-02-04T16:51:24+01:00"
    
    >>> import re
    >>> re.sub("[0-9]", "x", get_isodatetime())
    'xxxx-xx-xxTxx:xx:xx+xx:xx'
    """
    import datetime, time
    # Calculate the offset taking into account daylight saving time
    utc_offset_sec = time.altzone if time.localtime().tm_isdst else time.timezone
    utc_offset = datetime.timedelta(seconds=-utc_offset_sec)
    return datetime.datetime.now().replace(tzinfo=datetime.timezone(offset=utc_offset), microsecond=0).isoformat()


# ------ get used memory -- Windows version

import ctypes
from ctypes import wintypes

GetCurrentProcess = ctypes.windll.kernel32.GetCurrentProcess
GetCurrentProcess.argtypes = []
GetCurrentProcess.restype = wintypes.HANDLE

SIZE_T = ctypes.c_size_t

class PROCESS_MEMORY_COUNTERS_EX(ctypes.Structure):
    _fields_ = [
        ('cb', wintypes.DWORD),
        ('PageFaultCount', wintypes.DWORD),
        ('PeakWorkingSetSize', SIZE_T),
        ('WorkingSetSize', SIZE_T),
        ('QuotaPeakPagedPoolUsage', SIZE_T),
        ('QuotaPagedPoolUsage', SIZE_T),
        ('QuotaPeakNonPagedPoolUsage', SIZE_T),
        ('QuotaNonPagedPoolUsage', SIZE_T),
        ('PagefileUsage', SIZE_T),
        ('PeakPagefileUsage', SIZE_T),
        ('PrivateUsage', SIZE_T),
    ]

GetProcessMemoryInfo = ctypes.windll.psapi.GetProcessMemoryInfo
GetProcessMemoryInfo.argtypes = [
    wintypes.HANDLE,
    ctypes.POINTER(PROCESS_MEMORY_COUNTERS_EX),
    wintypes.DWORD,
]
GetProcessMemoryInfo.restype = wintypes.BOOL

def get_memory_info(process=None):
    """Return Win32 process memory counters structure as a dict."""
    if process is None:
        process = GetCurrentProcess()
    counters = PROCESS_MEMORY_COUNTERS_EX()
    ret = GetProcessMemoryInfo(process, ctypes.byref(counters), ctypes.sizeof(counters))
    if not ret:
        raise ctypes.WinError()
    info = dict((name, getattr(counters, name)) for name, _ in counters._fields_)
    return info

def used_mem():
    """
    Return memory usage as (working set, private bytes).
    Windows only for now.
    """
    try:
        m = get_memory_info()
        return m["WorkingSetSize"], m["PrivateUsage"]
    except:
        return None


def tester():
    from .prt import prt
    prt("test: util")
    import doctest, sys
    this_mod = sys.modules[__name__]    
    doctest.testmod(this_mod, verbose=False)
