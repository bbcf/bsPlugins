from bsPlugins import *
from bbcflib.bFlatMajor.stream import neighborhood, score_by_feature
from bbcflib.btrack import track
from bbcflib import genrep
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import numpy
import os
import itertools


ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100

__requires__ = ["numpy"]


class NormalizeForm(BaseForm):
    child = twd.HidingTableLayout()

    input_type = twd.HidingRadioButtonList(label_text='Input type',
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

    method = twf.RadioButton(label_text='Method',
        options=('total','deseq','quantile'),
        help_text='Select input type (Formatted table, or signal tracks)')

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
out_parameters = [{'id': 'differential_expression', 'type': 'file'}]


class NormalizePlugin(OperationPlugin):

    description = """
    """
    info = {
        'title': 'Normalization',
        'description': description,
        'path': ['Signal', 'DESeq'],
        'output': DESeqForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    # from bsPlugins.Normalize import NormalizePlugin; DESeqPlugin()(**{'table':'tests/DESeq/table.tab'})

    def __call__(self, **kw):
        self.new_file(out, 'normalized')
        return 1
