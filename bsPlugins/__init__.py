try:
    from tw2 import core as twc
    from tw2 import forms as twf
    from tw2 import dynforms as twd
except ImportError:
    from tw0 import twf, twc, twd


from base import OperationPlugin
try:
    from bs.operations.base import BaseForm, Multi, DynForm
except ImportError:
    try:
        class BaseForm(twf.TableForm):
            pass
    except AttributeError:
        class BaseForm():
            pass
    class DynForm():
        pass
    class Multi():
        pass

import os
filterplugins = lambda x: x.endswith('.py') and x not in ('__init__.py',)

PLUGINS_FILES = [str(f)[:-3] for p in __path__ for f in os.listdir(p) if filterplugins(f)]
