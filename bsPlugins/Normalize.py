from bsPlugins import *
from bbcflib.btrack import track
from bbcflib import genrep
import os
from bbcflib.bFlatMajor import common
from numpy import asarray, transpose

ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100

__requires__ = ["numpy"]


class NormalizeForm(BaseForm):
    child = twd.HidingTableLayout()

    input_type = twd.HidingRadioButtonList(label='Input type',
        options=('Table', 'Signal'),
        mapping={'Table':  ['table'],
                 'Signal': ['SigMulti','feature_type'],},
        help_text='Select input type (Formatted table, or signal tracks)')
    table = twf.FileField(label='Table: ',
        help_text='Select scores table',
        validator=twf.FileValidator(required=True))
    feature_type = twd.HidingSingleSelectField(label='Feature type: ',
        options=ftypes, prompt_text=None,
        mapping={ftypes[-1][0]: ['features'],
                 1: ['upstream', 'downstream']},
        help_text='Choose a feature set or upload your own',
        validator=twc.Validator(required=True))
    class SigMulti(Multi):
        signals = twf.FileField(label='Signal: ',
            help_text='Select signal file (position and score, e.g. bedgraph)',
            validator=twf.FileValidator(required=True))

#    method = twf.RadioButton(label='Method',
#        options=('total','deseq','quantile'),
#        help_text='Select input type (Formatted table, or signal tracks)')

    features = twf.FileField(label='Custom feature set: ',
        help_text='Select a feature file (e.g. bed)',
        validator=twf.FileValidator())
    upstream = twf.TextField(label='Promoter upstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_up_def,
        help_text='Size of promoter upstream of TSS')
    downstream = twf.TextField(label='Promoter downstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_down_def,
        help_text='Size of promoter downstream of TSS')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        validator=twc.Validator(required=True),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Quantify")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id': 'input_type', 'type': 'radio'},
        {'id': 'signals', 'type': 'track', 'required': True, 'multiple': True},
        {'id': 'table', 'type': 'txt', 'required': True, 'multiple': True},
        {'id': 'feature_type', 'type': 'int'},
        {'id': 'upstream', 'type': 'int'},
        {'id': 'downstream', 'type': 'int'},
        {'id': 'assembly', 'type': 'assembly'},
        {'id': 'features', 'type': 'userfile'},
        {'id': 'method', 'type': 'radio'}
]
out_parameters = [{'id': 'normalized', 'type': 'file'}]


class NormalizePlugin(OperationPlugin):

    description = """
    """
    info = {
        'title': 'Normalization',
        'description': description,
        'path': ['Signal', 'DESeq'],
        'output': NormalizeForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):

        assembly = genrep.Assembly(kw.get('assembly'))
        chrmeta = assembly.chrmeta or "guess"

        if kw.get('input_type') == 'Table':
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
#        else:
#            from QuantifyTable import QuantifyTablePlugin
#            kw['score_op'] = 'sum'
#            table = QuantifyTablePlugin().quantify(**kw)
#            signals = kw.get('signals',[])

        norm = common.normalize(asarray(matrix).transpose(), 'total')
        output = self.temporary_path(fname='output.tab')
        out = open(output, "w")
        out.write(header)
        for i in range(len(norm[0])):
            out.write(str(id[i])+"\t"+str(map(lambda x: "%.2g" % x, list(norm.transpose()[i]))).replace("'","").replace("[","").replace("]","").replace(", ","\t")+"\n")
        self.new_file(output, 'normalized')
        return self.display_time()
