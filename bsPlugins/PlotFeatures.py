from bsPlugins import *
from bbcflib.bFlatMajor.numeric import feature_matrix
from bbcflib.bFlatMajor.figure import heatmap, lineplot
from bbcflib.btrack import track
from numpy import vstack, concatenate

prom_up_def = 1000
prom_down_def = 100
plot_types = [(0, 'heatmap'), (1, 'average lineplot'), (2, 'mosaic plot')]
max_pages = 200

class PlotFeaturesForm(BaseForm):

    class SigMulti(Multi):
        signals = twf.FileField(label='Signal: ',
            help_text='Select signal file (e.g. bedgraph)',
            validator=twf.FileValidator(required=True))

    features = twf.FileField(label='Features: ',
        help_text='Select a feature file (e.g. bed)',
        validator=twf.FileValidator())

    mode = twf.SingleSelectField(label='Plot type: ',
        options=plot_types,
        prompt_text=None,
        validator=twc.Validator(required=True))
    upstream = twf.TextField(label='Upstream flank: ',
        validator=twc.IntValidator(),
        value=prom_up_def,
        help_text='Size of upstream flank in bp')
    downstream = twf.TextField(label='Downstream flank: ',
        validator=twc.IntValidator(),
        value=prom_down_def,
        help_text='Size of downstream flank in bp')
    submit = twf.SubmitButton(id="submit", value="Plot")

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'signals', 'type': 'track', 'multiple': True, 'required': True},
                 {'id': 'features', 'type': 'userfile'},
                 {'id': 'mode', 'type': 'list', 'required': True},
                 {'id': 'upstream', 'type': 'int'},
                 {'id': 'downstream', 'type': 'int'}]
out_parameters = [{'id': 'plot_features', 'type': 'pdf'}]


class PlotFeaturesPlugin(OperationPlugin):

    info = {
        'title': 'Plot signal',
        'description': 'Plot signal from a selection of regions',
        'path': ['Plot', 'Plot features'],
        'output': PlotFeaturesForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        chrmeta = "guess"
        features = track(kw.get('features'), chrmeta=chrmeta)
        signals = kw.get('signals', [])
        if not isinstance(signals, list): signals = [signals]
        signals = [track(sig) for sig in signals]
        labels = None
        data = None
        for chrom in features.chrmeta:
            _l, _d = feature_matrix([s.read(chrom) for s in signals],
                                    features.read(chrom), segment=True)
            if _d.size == 0:
                continue
            if data is None:
                labels = _l
                data = _d
            else:
                labels = concatenate((labels, _l))
                data = vstack((data, _d))
        pdf = self.temporary_path(fname='plot_features', ext='.pdf')
        if data is None:
            raise ValueError("No data")
        kw['mode'] = int(kw.get('mode', 0))
        if kw['mode'] == 0:
            new = True
            for n in range(data.shape[-1] - 1):
                heatmap(data[:, :, n], output=pdf, new=new, last=False,
                        rows=labels, orderRows=True, orderCols=False)
                new = False
            heatmap(data[:, :, -1], output=pdf, new=new, last=True,
                    rows=labels, orderRows=True, orderCols=False)
        elif kw['mode'] == 1:
            X = range(data.shape[1])
            Y = data.mean(axis=0)
            lineplot(X, [Y[:, n] for n in range(data.shape[-1])],
                     output=pdf, new=True, last=True)
        elif kw['mode'] == 2:
            X = range(data.shape[1])
            new = True
            mfrow = [4, 3]
            nplot = min(data.shape[0], max_pages * mfrow[0] * mfrow[1])
            for reg in range(nplot - 1):
                lineplot(X, [data[reg, :, n] for n in range(data.shape[-1])],
                         output=pdf, new=new, last=False, mfrow=mfrow)
                new = False
                mfrow = []
            lineplot(X, [data[nplot - 1, :, n] for n in range(data.shape[-1])],
                     output=pdf, new=new, last=True)
        else:
            raise ValueError("Mode not implemented: %s" % kw['mode'])
        self.new_file(pdf, 'plot_features')
        return 1
