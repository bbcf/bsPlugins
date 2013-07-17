from bsPlugins import *
from bbcflib.gfminer.common import sorted_stream
from bbcflib.gfminer.stream import neighborhood, score_by_feature
from bbcflib.gfminer.numeric import score_array, correlation
from bbcflib.gfminer.figure import pairs
from bbcflib.track import track
from bbcflib import genrep
from numpy import vstack, array

ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100
plot_types = [(0, 'correlations'), (1, 'density plots')]
cormax = 500

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'input_type', 'type': 'radio'},
                 {'id': 'signals', 'type': 'track', 'required': True, 'multiple':'SigMulti'},
                 {'id': 'feature_type', 'type': 'list'},
                 {'id': 'features', 'type': 'track'},
                 {'id': 'upstream', 'type': 'int', 'required': True},
                 {'id': 'downstream', 'type': 'int', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'mode', 'type': 'list', 'required': True}]
out_parameters = [{'id': 'plot_pairs', 'type': 'pdf'}]


class PairsPlotForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signal: '
        signals = twb.BsFileField(label=' ',
                                help_text='Select signal file (.g. bedgraph)',
                                validator=twb.BsFileFieldValidator(required=True))
    child = twd.HidingTableLayout()
    feature_type = twd.HidingSingleSelectField(label='Feature type: ',
                                               options=ftypes, prompt_text=None,
                                               mapping={ftypes[-1][0]: ['features'],
                                                        1: ['upstream', 'downstream']},
                                               help_text='Choose a feature set or upload your own',
                                               validator=twc.Validator(required=True))
    features = twb.BsFileField(label='Custom feature set: ',
                             help_text='Select a feature file (e.g. bed)',
                             validator=twb.BsFileFieldValidator())
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
                                     help_text='Reference genome')
    mode = twf.SingleSelectField(label='Plot type: ',
                                 options=plot_types,
                                 prompt_text=None)
    submit = twf.SubmitButton(id="submit", value="Plot")


class PairsPlotPlugin(BasePlugin):
    """Plots pairwise comparisons between signal tracks."""
    info = {
        'title': 'Pairwise plots',
        'description': __doc__,
        'path': ['Graphics', 'Plot pairs'],
        'output': PairsPlotForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        feature_type = int(kw.get('feature_type') or 0)
        assembly_id = kw.get('assembly') or None
        chrmeta = "guess"
        if assembly_id:
            assembly = genrep.Assembly(assembly_id)
            chrmeta = assembly.chrmeta
            genes = assembly.gene_track
            exons = assembly.exon_track
        elif not(feature_type == 3):
            raise ValueError("Please specify an assembly")
        signals = kw.get('SigMulti',{}).get('signals', [])
        if not isinstance(signals, list): signals = [signals]
        signals = [track(sig, chrmeta=chrmeta) for sig in signals]
        snames = [sig.name for sig in signals]
        if feature_type == 0: #bodies
            features = genes
        elif feature_type == 1: #promoters
            prom_pars = {'before_start': int(kw.get('upstream') or prom_up_def),
                         'after_start': int(kw.get('downstream') or prom_down_def),
                         'on_strand': True}
            features = lambda c: neighborhood(genes(c), **prom_pars)
        elif feature_type == 2: #exons
            features = exons
        elif feature_type == 3: #custom track
            _t = track(kw.get('features'), chrmeta=chrmeta)
            chrmeta = _t.chrmeta
            features = _t.read
        else:
            raise ValueError("Feature type not known: %i" % feature_type)
        pdf = self.temporary_path(fname='plot_pairs.pdf')
        narr = None
        if int(kw['mode']) == 0: #correl
            xarr = array(range(-cormax, cormax + 1))
            srtdchrom = sorted(chrmeta.keys())
            features = [x[:3] for chrom in srtdchrom
                        for x in sorted_stream(features(chrom))]
            _f = ['chr', 'start', 'end', 'score']
            narr = correlation([s.read(fields=_f) for s in signals],
                               features, (-cormax, cormax), True)
        elif int(kw['mode']) == 1: #density
            xarr = None
            for chrom in chrmeta:
                feat = features(chrom)
                means = score_by_feature([s.read(chrom) for s in signals], feat)
                mf = means.fields[len(feat.fields):]
                _n, _l = score_array(means, mf)
                if _n.size == 0: continue
                if narr is None: narr = _n
                else:            narr = vstack((narr, _n))
        else:
            raise ValueError("Mode not implemented: %s" % kw['mode'])
        if narr is None:
            raise ValueError("No data")
        pairs(narr, xarr, labels=snames, output=pdf)
        self.new_file(pdf, 'plot_pairs')
        return self.display_time()

