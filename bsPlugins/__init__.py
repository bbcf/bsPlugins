try:
    from tw2 import core as twc
    from tw2 import forms as twf
    from tw2 import dynforms as twd
    from tw2 import bs as twb
except ImportError:
    from tw0 import twf, twc, twd, twb


from base import BasePlugin
try:
    from bs.operations.base import BaseForm, Multi, DynForm
except ImportError:
    try:
        class BaseForm(twf.TableForm):
            pass
    except AttributeError:
        class BaseForm():
            pass
    try:
        class Multi(twd.GrowingGridLayout):
            pass
    except AttributeError:
        class Multi():
            pass

    class DynForm():
        pass

import os, glob
PLUGINS_FILES = [os.path.basename(f)[:-3] for p in __path__ for f in glob.glob(p+"/[a-zA-Z]*.py")]
