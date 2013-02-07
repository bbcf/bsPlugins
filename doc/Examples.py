from bs.operations.base import OperationPlugin  # import the base class to build your plugin
import os
import random


# some shared information about these exemples
meta = {'version': "1.0.1",
        'author': "Yohan Jarosz",
        'contact': "webmaster-bbcf@epfl.ch"}

__install_requires__ = []  # beta


class Simple(OperationPlugin):

    info = {
        'title': 'WriteFile',                                                                 # The title of your operation
        'description': """As an exemple, this plugin writes the input you give in a text file.
                          You can also put <a href='https://github.com/bbcf/bs.operations/blob/master/Examples.py'>
                          links</a> to some documentation.
        """,                                                                                  # Describe the operation's goal
        'path': ['Examples', 'Automatics forms', 'Write output to a file'],                   # Under which category the operation will be set
                                                                                              # First in the list mean higher category
        'in': [{'id': 'input', 'type': 'text', 'required': True}],    # All input parameters
        'out': [{'id': 'output', 'type': 'file'}],                    # All output parameters
        'meta': meta,                                # Meta information (authors, version, ...)
    }

    def __call__(self, *args, **kw):
        text = kw.get('input', '')                    # get the parameter back

        path = self.temporary_path()                  # get a temporary path
        with open(path, 'w') as f:                    # open a file & write the input
            f.write(text)
        self.new_file(path, 'output')                 # add a file to the result

        return 1


class ReadFile(OperationPlugin):

    info = {
        'title': 'ReadFile',
        'description': "Read a file to find its number of characters and write the result in an output file. "
                        "If randomize is checked, the result is between 0 and the file's number of characters ",
        'path': ['Examples', 'Automatics forms', 'Read a file'],
        'in': [{'id': 'fname', 'label': 'File Name', 'type': 'text'},
                {'id': 'input', 'label': 'File input', 'type': 'file', 'required': True},
                {'id': 'randomize', 'label': 'Randomize', 'type': 'boolean'}
        ],
        'out': [{'id': 'output', 'type': 'file'}],
        'meta': meta,
        }

    def __call__(self, *args, **kw):
        # get parameters
        fin = kw.get('input')
        fname = kw.get('fname', '')
        if fname == '':
            fname = os.path.split(fin)[1]
        randomize = kw.get('randomize', False)
        # read the input file
        nchar = 0
        with open(fin, 'r') as rin:
            for line in rin:
                nchar += len(line)
                # randomize if checked
        if randomize:
            nchar = random.randint(0, nchar - 1)
            # take a temporary path where to write the file
        fout = self.temporary_path(fname)
        # write the result
        with open(fout, 'w') as wout:
            wout.write(str(nchar))
            # add the file to the result
        self.new_file(fout, 'output')
        return self.display_time()

class Error(OperationPlugin):
    info = {
        'title': 'Error plugin',
        'description': 'This plugin ends up with an error.',
        'path': ['Examples', 'Automatics forms', 'Error'],
        'in': [{'id': 'a', 'label': 'parameter "a"', 'required': False, 'type': 'text'}],
        'out': [{'id': 'output', 'type': 'text'}],
        'meta': meta,
    }

    def __call__(self, *args, **kw):
        1 / 0
        return 1
