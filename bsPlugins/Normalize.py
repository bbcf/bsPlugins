from bsPlugins import *
from bbcflib.gfminer import common
from bbcflib.track import track
from bbcflib import genrep
from numpy import asarray,transpose
import os

__requires__ = ["numpy"]

class NormalizeForm(BaseForm):

    child = twd.HidingTableLayout()
    table = twb.BsFileField(
        label='Table: ',
        help_text='Select scores table',
        validator=twb.BsFileFieldValidator(required=True))
    method = twf.RadioButtonList(
        label='Method: ',
        options=['total','deseq','quantile'],
        help_text='Select the normalization method')
    submit = twf.SubmitButton(id="submit", value="Submit")

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id': 'table', 'type': 'txt', 'required': True, 'multiple': True},
        {'id': 'method', 'type': 'radio'}]

out_parameters = [{'id': 'normalized', 'type': 'file'}]

class NormalizePlugin(BasePlugin):
    """ Normalize the columns of a tab-delimited file using a specified method and returns a
normalized tab-delimited file."""
    info = {
        'title': 'Normalization',
        'description': __doc__,
        'path': ['Signal', 'Normalization'],
        'output': NormalizeForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):

        filename = kw.get('table')
        assert os.path.exists(str(filename)), "File not found: '%s'" % filename
        file = open(filename, 'r')
        header = file.readline()
        id = []
        matrix = []
        for line in file:
            newline = line.split()
            id.append(newline[0])
            matrix.append(map(int, newline[1:len(header)]))
        norm = common.normalize(asarray(matrix).transpose(), kw.get('method'))
        output = self.temporary_path(fname='output.tab')
        out = open(output, "w")
        out.write(header)
        for i in range(len(norm[0])):
            out.write(str(id[i])+"\t"+str(map(lambda x: "%.2g" % x, list(norm.transpose()[i]))).replace("'","").replace("[","").replace("]","").replace(", ","\t")+"\n")
        self.new_file(output, 'normalized')
        return self.display_time()
