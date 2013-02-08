from bsPlugins import *
from bbcflib.bFlatMajor.stream import neighborhood, score_by_feature
from bbcflib.btrack import track
from bbcflib import genrep

prom_up_def = 1000
prom_down_def = 100
ftypes = [(0, 'gene bodies'), (1, 'gene promoters'), (2, 'custom upload')]
funcs = ['mean', 'sum', 'median', 'min', 'max']

class QuantifyTableForm(BaseForm):

    class SigMulti(Multi):
        signals = twf.FileField(label='Signal: ',
                                help_text='Select signal file (e.g. bedgraph)',
                                validator=twf.FileValidator(required=True))

    score_op = twf.SingleSelectField(label_text='Score operation: ',
                                     options=funcs, prompt_text=None,
                                     help_text='Operation performed on scores within each feature')
    child = twd.HidingTableLayout()
    feature_type = twd.HidingSingleSelectField(label='Feature type: ',
        options=ftypes, prompt_text=None,
        mapping={ftypes[-1][0]: ['features'],
                 1: ['upstream', 'downstream']},
        help_text='Choose a feature set or upload your own',
        validator=twc.Validator(required=True))
    features = twf.FileField(label='Custom feature set: ',
        help_text='Select a feature file (e.g. bed)',
        validator=twf.FileValidator(required=True))
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    upstream = twf.TextField(label='Promoter upstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_up_def,
        help_text='Size of promoter upstream of TSS')
    downstream = twf.TextField(label='Promoter downstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_down_def,
        help_text='Size of promoter downstream of TSS')
    submit = twf.SubmitButton(id="submit", value="Quantify")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'signals', 'type': 'track', 'multiple': True, 'required': True},
                 {'id': 'feature_type', 'type': 'list'},
                 {'id': 'features', 'type': 'userfile'},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'upstream', 'type': 'int', 'required': True},
                 {'id': 'downstream', 'type': 'int', 'required': True}]
out_parameters = [{'id': 'features_quantification', 'type': 'track'}]


class QuantifyTablePlugin(OperationPlugin):
    info = {
        'title': 'Signal quantification',
        'description': 'Quantify a signal track on a set of regions',
        'path': ['Signal', 'Quantify features'],
        'output': QuantifyTableForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        feature_type = int(kw.get('feature_type', 0))
        func = str(kw.get('score_op', 'mean'))
        assembly_id = kw.get('assembly') or None
        chrmeta = "guess"
        if assembly_id:
            assembly = genrep.Assembly(assembly_id)
            chrmeta = assembly.chrmeta
            genes = assembly.gene_track
        elif not(feature_type == 2):
            raise ValueError("Please specify an assembly")
        signals = kw.get('signals', [])
        if not isinstance(signals, list):
            signals = [signals]
        signals = [track(sig, chrmeta=chrmeta) for sig in signals]
        if feature_type == 0:
            features = genes
        elif feature_type == 1:
            prom_pars = {'before_start': int(kw.get('upstream') or prom_up_def),
                          'after_start': int(kw.get('downstream') or prom_down_def),
                          'on_strand': True}
            features = lambda c: neighborhood(genes(c), **prom_pars)
        elif feature_type == 2:
            _t = track(kw.get('features'), chrmeta=chrmeta)
            chrmeta = _t.chrmeta
            features = _t.read
        else:
            return 2
        output = self.temporary_path(fname='features_quantification.sql')
        if len(signals) > 1:
            _f = ["score" + str(i) for i in range(len(signals))]
        else:
            _f = ["score"]
        tout = track(output, format='sql', fields=['start', 'end', 'name'] + _f,
                     chrmeta=chrmeta, info={'datatype': 'qualitative'})
        for chrom in chrmeta:
            sread = [sig.read(chrom) for sig in signals]
            tout.write(score_by_feature(sread, features(chrom), fn=func),
                       chrom=chrom, clip=True)
        tout.close()
        self.new_file(output, 'features_quantification')
        return 1
