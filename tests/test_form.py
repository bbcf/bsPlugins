#!/usr/bin/env python

import sys
import bsPlugins

def formtest(plugin):
    imp = "from bsPlugins.%s import %sForm" % (plugin,plugin)
    test = "bsPlugins.base.test_form(%sForm, port=8080)" % plugin
    exec(imp)
    exec(test)

if __name__ == '__main__':
    plugin = sys.argv[1]
    print "Test plugin %s" % plugin
    sys.exit(formtest(plugin))

