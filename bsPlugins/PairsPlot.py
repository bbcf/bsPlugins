from bsPlugins import *
from bbcflib.gfminer.common import sorted_stream, add_name_field
from bbcflib.gfminer.stream import neighborhood, score_by_feature
from bbcflib.gfminer.numeric import score_array, correlation
from bbcflib.gfminer.figure import pairs
from bbcflib.track import track
from bbcflib import genrep
from numpy import vstack, array

ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100
plot_types = [(0, 'density plots'), (1, 'correlations')]
_cormax = 500
_MAX_PLOTS_ = 100

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

fmap = {ftypes[-1][0]: ['features'], 1: ['upstream', 'downstream']},

in_parameters = [{'id': 'signals', 'type': 'track', 'required': True, 'multiple': True, 'label': 'Signal: ', 'help_text': 'Select signal file (e.g. bedgraph)'},
                 {'id': 'feature_type', 'type': 'list', 'required': True, 'label': 'Feature type: ', 'help_text': 'Choose a feature set or upload your own', 'options': ftypes, 'mapping': {ftypes[-1][0]: ['features'], 1: ['upstream', 'downstream']}, 'prompt_text': None},
                 {'id': 'features', 'type': 'track', 'label':'Custom feature set: ', 'help_text':'Select a feature file (e.g. bed)'},
                 {'id': 'upstream', 'type': 'int', 'required': True, 'label':'Promoter upstream distance: ', 'help_text':'Size of promoter upstream of TSS', 'value':prom_up_def},
                 {'id': 'downstream', 'type': 'int', 'required': True, 'label':'Promoter downstream distance: ', 'help_text':'Size of promoter downstream of TSS', 'value':prom_down_def},
                 {'id': 'assembly', 'type': 'assembly', 'label': 'Assembly: ', 'help_text':'Reference genome','options':genrep.GenRep().assemblies_available()},
                 {'id': 'highlights', 'type': 'track', 'multiple': 'HiMulti', 'label':'features to highlight: ', 'help_text':'Select a feature file (e.g. bed)'},
                 {'id': 'mode', 'type': 'list', 'required': True, 'label': 'Plot type: ', 'options': plot_types, 'mapping':{1: ['cormax','individual']},'prompt_text':None},
                 {'id': 'cormax', 'type': 'int', 'label':'Spatial range: ', 'help_text':'Maximum lag in bp to compute correlations', 'value': _cormax},
                 {'id': 'individual', 'type': 'boolean', 'label': 'Individual: ','help_text':'Plot each region individually (default: false)', 'value':False }]

out_parameters = [{'id':'table', 'type':'file'},
                  {'id': 'plot_pairs', 'type': 'pdf'}]

class PairsPlotForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signal: '
        signals = twb.BsFileField(label=' ',
                                  help_text='Select signal file (e.g. bedgraph)',
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
                             validator=twc.IntValidator(),
                             value=prom_up_def,
                             help_text='Size of promoter upstream of TSS')
    downstream = twf.TextField(label='Promoter downstream distance: ',
                               validator=twc.IntValidator(),
                               value=prom_down_def,
                               help_text='Size of promoter downstream of TSS')
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=genrep.GenRep().assemblies_available(),
                                     help_text='Reference genome')
    class HiMulti(twb.BsMultiple):
        label='features to highlight: '
        highlights = twb.BsFileField(label=' ',
                                     help_text='Select a feature file (e.g. bed)',
                                     validator=twb.BsFileFieldValidator())
    mode = twd.HidingSingleSelectField(label='Plot type: ',
                                       options=plot_types,
                                       mapping={1: ['cormax','individual']},
                                       prompt_text=None)
    cormax = twf.TextField(label='Spatial range: ',
                           validator=twc.IntValidator(),
                           value=_cormax,
                           help_text='Maximum lag in bp to compute correlations')
    individual = twf.CheckBox(label='Individual: ',
                              value=False,
                              help_text='Plot each region individually (default: false)')
    submit = twf.SubmitButton(id="submit", value="Plot")


class PairsPlotPlugin(BasePlugin):
    """Plots pairwise comparisons between signal tracks:

* For *density plots* each signal track is quantified at the selected features, and this data is represented as two-way scatter plots (above diagonal), histograms (on the diagonal), and quantile plots (below diagonal).
* *Correlation* plots show spatial auto- and cross-correlation of signals within the selected features. If the individual option is selected, each region is plotted separately (max 100 plots) and a table containing the maximum correlation values as well as their corresponding lag values is given.
"""
    info = {
        'title': 'Pairwise plots',
        'description': __doc__,
        'path': ['Graphics', 'Plot pairs'],
#        'output': PairsPlotForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        feature_type = int(kw.get('feature_type') or 0)
        individual = kw.get('individual',False)
        if isinstance(individual, basestring):
            individual = (individual.lower() in ['1', 'true', 't', 'on'])
        if individual and int(kw['mode']) != 1:
            raise ValueError("Only correlation plots can work with the 'individual' option.")

        assembly_id = kw.get('assembly') or None
        chrmeta = "guess"
        if assembly_id:
            assembly = genrep.Assembly(assembly_id)
            chrmeta = assembly.chrmeta
            genes = assembly.gene_track
            exons = assembly.exon_track
        elif not(feature_type == 3):
            raise ValueError("Please specify an assembly")
        #signals = kw.get('SigMulti',{}).get('signals', [])
        signals = kw.get('signals', [])
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
        #highlights = kw.get('HiMulti',{}).get('highlights', [])
        highlights = kw.get('highlights', [])
        if not isinstance(highlights, list): highlights = [highlights]
        if highlights is not None:
            highlights = [track(hi, chrmeta=chrmeta) for hi in highlights]
            hinames = [t.name for t in highlights]
        pdf = self.temporary_path(fname='plot_pairs.pdf')
        narr = None
        set_index = []
        set_labels = []
        _new = True
        if int(kw['mode']) == 1: #correl
            cormax = int(kw.get('cormax') or _cormax)
            xarr = array(range(-cormax, cormax + 1))
            _f = ['chr', 'start', 'end', 'score']
            features = [x[:3] for chrom in chrmeta
                        for x in sorted_stream(features(chrom))]
            table = self.temporary_path(fname='table.txt')
            with open(table,"w") as t:
                t.write("\t".join(["chr","start","end","max(correlation)","lag_max"])+"\n")
                if individual:
                    for nplot,feature in enumerate(features):
                        if (narr is not None and nplot < _MAX_PLOTS_):
                            pairs(narr, xarr, labels=snames, output=pdf, new=_new, last=False)
                            _new = False
                        narr = correlation([s.read(fields=_f) for s in signals], [feature], (-cormax, cormax), True)
                        list_corr = list(narr[0][0])
                        max_corr = max(list_corr)
                        lag_max = list_corr.index(max_corr)-cormax
                        t.write("\t".join([str(x) for x in
                                           feature[:3]+(max_corr,lag_max)])+"\n")
                else:
                    narr = correlation([s.read(fields=_f) for s in signals], features, (-cormax, cormax), True)
                    list_corr = list(narr[0][0])
                    max_corr = max(list_corr)
                    lag_max = list_corr.index(max_corr)-cormax
                    t.write("\t".join(["-","-","-"]+[str(max_corr),str(lag_max)])+"\n")
        elif int(kw['mode']) == 0: #density
            xarr = None
            for chrom in chrmeta:
                feat = features(chrom)
                if 'name' not in feat.fields:
                    feat = add_name_field(feat)
                means = score_by_feature([s.read(chrom) for s in signals], feat)
                mf = means.fields[len(feat.fields):]
                _n, _l = score_array(means, mf)
                if _n.size == 0: continue
                if narr is None: narr = _n
                else:            narr = vstack((narr, _n))
            set_index = [narr.shape[0]]
            for hitrack in highlights:
                for chrom in chrmeta:
                    hiread = hitrack.read(chrom)
                    if 'name' not in hiread.fields:
                        hiread = add_name_field(hiread)
                    means = score_by_feature([s.read(chrom) for s in signals], hiread)
                    mf = means.fields[len(hiread.fields):]
                    _n, _l = score_array(means, mf)
                    if _n.size == 0: continue
                    narr = vstack((narr, _n))
                    set_labels.extend(_l)
                set_index.append(narr.shape[0])
        else:
            raise ValueError("Mode not implemented: %s" % kw['mode'])
        if narr is None:
            raise ValueError("No data")
        pairs(narr, xarr, labels=snames, output=pdf, highlights=[set_index,set_labels], new=_new, last=True)
        if int(kw['mode']) == 1: self.new_file(table, 'table')
        self.new_file(pdf, 'plot_pairs')
        return self.display_time()

