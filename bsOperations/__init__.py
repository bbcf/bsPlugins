try:
    from tw2 import core as twc
    from tw2 import forms as twf
    from tw2 import dynforms as twd
except ImportError: 
    from tw0 import twf, twc, twd

try:
    from bs.operations.base import *
except ImportError: 
    from base import *
