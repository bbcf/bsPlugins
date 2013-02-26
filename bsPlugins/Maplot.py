from bsPlugins import *
from bbcflib.btrack import track
from bbcflib import genrep
import os, shutil
from maplot import MAplot

__requires__ = ["numpy"]

ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100


class MaplotForm(BaseForm):
    child = twd.HidingTableLayout()

    input_type = twd.HidingRadioButtonList(label='Input type',
        options=('Table', 'Signal'),
        mapping={'Table':  ['table'],
                 'Signal': ['SigMulti','feature_type','assembly'],},
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
    submit = twf.SubmitButton(id="submit", value="MA-plot")


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
]
out_parameters = [{'id': 'MA-plot', 'type': 'file'}]


class MaplotPlugin(OperationPlugin):

    description = """
    """
    info = {
        'title': 'MA-plot',
        'description': description,
        'path': ['Signal', 'MA-plot'],
        'output': MaplotForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):

        if kw.get('input_type') == 'Table':
            table = kw.get('table')
            assert os.path.exists(str(table)), "File not found: '%s'" % table
            with open(table) as t:
                line1 = t.readline()
                nscores = len(line1.split())-1
        else:
            from QuantifyTable import QuantifyTablePlugin
            kw['score_op'] = 'sum'
            table = QuantifyTablePlugin().quantify(**kw)
            signals = kw.get('signals',[])
            nscores = len(signals)

        #table = track(table, format='txt', fields=["name"]+['score'+str(i) for i in range(nscores)])
        output_filename = MAplot(table)
        output = self.temporary_path(fname='maplot.png')
        shutil.copy(output_filename,output)
        self.new_file(output, 'MA-plot')
        return self.display_time()
