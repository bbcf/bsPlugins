###########
# IMPORTS #
###########
import os
import random
from bsPlugins import OperationPlugin


############
# METADATA #
############

meta = {'version': "1.0.1",
        'author': "Yohan Jarosz",
        'contact': "webmaster-bbcf@epfl.ch"}

######################
# PLUGIN INFORMATION #
######################
parameters = {'in': [{'id': 'input', 'type': 'text', 'required': True}, ],
              'out': [{'id': 'output', 'type': 'file'}, ]}

plugin_information = {
    'title': 'WriteFile',

    'description': """As an exemple, this plugin writes the input you give in a text file.
    You can also put <a href='https://github.com/bbcf/bs.operations/blob/master/Examples.py'>links</a> to some documentation or
    other <br/> HTML tags.""",

    'path': ['Examples', 'Automatics forms', 'Write output to a file'],

    'meta': meta,

    'in': parameters['in'],
    'out': parameters['out']
}

########
# FORM # skipped
########

##########
# PLUGIN #
##########


class Simple(OperationPlugin):

    info = plugin_information

    def __call__(self, *args, **kw):
        text = kw.get('input', '')                    # get the parameter back

        path = self.temporary_path()                  # get a temporary path
        with open(path, 'w') as f:                    # open a file & write the input
            f.write(text)
        self.new_file(path, 'output')                 # add a file to the result

        return self.display_time()                    # display the time elapsed
