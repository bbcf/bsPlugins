from bsPlugins import *
from bbcflib.gfminer.numeric import feature_matrix
from bbcflib.gfminer.figure import heatmap, lineplot
from bbcflib.track import track
from numpy import vstack, concatenate, array
import os, tarfile

_nbins = 50
_upstr = (.1,5)
_downstr = (.1,5)
prom_up_def = 1000
prom_down_def = 100
plot_types = [(0, 'heatmap'), (1, 'average lineplot'), (2, 'mosaic plot')]
max_pages = 200
output_list = ['pdf','archive']

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'signals', 'type': 'track', 'multiple': 'SigMulti', 'required': True},
                 {'id': 'features', 'type': 'track'},
                 {'id': 'mode', 'type': 'list', 'required': True},
                 {'id': 'upstream', 'type': 'int'},
                 {'id': 'downstream', 'type': 'int'},
                 {'id': 'nbins', 'type': 'int'},
                 {'id': 'output', 'type': 'list', 'required': True}]
out_parameters = [{'id': 'plot_features', 'type': 'pdf'},
                  {'id':'data_archive', 'type':'file'}]


class PlotFeaturesForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signal: '
        signals = twb.BsFileField(label=' ',
                                  help_text='Select signal file (e.g. bedgraph)',
                                  validator=twb.BsFileFieldValidator(required=True))
    features = twb.BsFileField(label='Features: ',
                               help_text='Select a feature file (e.g. bed)',
                               validator=twb.BsFileFieldValidator(required=True))

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
    nbins = twf.TextField(label='Number of bins: ',
                               validator=twc.IntValidator(),
                               value=_nbins,
                               help_text='Number of bins each feature is divided into')
    output = twf.SingleSelectField(label='Output: ', options=output_list,
                                   prompt_text=None, help_text='Pdf only or data+pdf')
    submit = twf.SubmitButton(id="submit", value="Plot")


class PlotFeaturesPlugin(BasePlugin):
    """Plot signals from a selection of regions."""
    info = {
        'title': 'Plot signals in genomic regions',
        'description': __doc__,
        'path': ['Graphics', 'Plot features'],
        'output': PlotFeaturesForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        chrmeta = "guess"
        features = track(kw.get('features'), chrmeta=chrmeta)
        signals = kw.get('SigMulti',{}).get('signals', [])
        if not isinstance(signals, list): signals = [signals]
        signals = [track(sig) for sig in signals]
        snames = [sig.name for sig in signals]
        labels = None
        data = None
        upstr = _upstr
        downstr = _downstr
        if kw.get("upstream") is not None:
            _up = int(kw["upstream"])
            if _up > 50: upstr = (_up,5)
            elif _up > 0: upstr = (_up,1)
            else: upstr = (0,0)
        if kw.get("downstream") is not None:
            _down = int(kw["downstream"])
            if _down > 50: downstr = (_down,5)
            elif _down > 0: downstr = (_down,1)
            else: downstr = (0,0)
        if kw.get("nbins") is not None: nbins = max(1,int(kw["nbins"]))
        else: nbins = _nbins
        for chrom in features.chrmeta:
            _l, _d = feature_matrix([s.read(chrom) for s in signals],
                                    features.read(chrom), segment=True,
                                    nbins=nbins, upstream=upstr, downstream=downstr)
            if _d.size == 0:
                continue
            if data is None:
                labels = _l
                data = _d
            else:
                labels = concatenate((labels, _l))
                data = vstack((data, _d))
        outf = str(kw.get('output'))
        if outf not in output_list:
            outf = output_list[0]
        pdf = self.temporary_path(fname='plot_features.pdf')
        if outf == 'archive':
            tarname = self.temporary_path(fname='plot_features.tar.gz')
            tarfh = tarfile.open(tarname, "w:gz")
        if data is None:
            raise ValueError("No data")
        mode = kw.get('mode', 0)
        if str(mode) in [str(x[0]) for x in plot_types]:
            mode = int(mode)
        X = array(range(-upstr[1]+1,nbins+downstr[1]+1))/(1.0*nbins)
        if mode in plot_types[0]: #heatmap
            new = True
            for n in range(data.shape[-1]-1):
                heatmap(data[:, :, n], output=pdf, new=new, last=False,
                        rows=labels, columns=X, main=snames[n],
                        orderRows=True, orderCols=False)
                new = False
            heatmap(data[:, :, -1], output=pdf, new=new, last=True,
                    rows=labels,  columns=X, main=snames[-1],
                    orderRows=True, orderCols=False)
            if outf == 'archive':
                for n,sn in enumerate(snames):
                    _datf = self.temporary_path(fname=sn+"_data.txt")
                    with open(_datf,"w") as dff:
                        dff.write("\t".join([""]+[str(x) for x in X])+"\n")
                        for k,l in enumerate(labels):
                            dff.write("\t".join([l]+[str(x) for x in data[k, :, n]])+"\n")
                    tarfh.add(_datf,arcname=os.path.basename(_datf))
        elif mode in plot_types[1]: #average lineplot
            Y = data.mean(axis=0)
            ymin = min([x.min() for x in Y]+[0])
            ymax = max([x.max() for x in Y])
            lineplot(X, [Y[:, n] for n in range(data.shape[-1])],
                     output=pdf, new=True, last=True, legend=snames, ylim=(ymin,ymax))
            if outf == 'archive':
                _datf = self.temporary_path(fname="lineplot_data.txt")
                with open(_datf,"w") as dff:
                    dff.write("\t".join([""]+[str(x) for x in X])+"\n")
                    for n,sn in enumerate(snames):
                        dff.write("\t".join([sn]+[str(x) for x in Y[:, n]])+"\n")
                tarfh.add(_datf,arcname=os.path.basename(_datf))
        elif mode in plot_types[2]: #mosaic
            new = True
            mfrow = [4, 3]
            nplot = min(data.shape[0], max_pages*mfrow[0]*mfrow[1])
            ymin = min([data.min(),0])
            ymax = data.max()
            for reg in range(nplot-1):
                lineplot(X, [data[reg, :, n] for n in range(data.shape[-1])],
                         output=pdf, new=new, last=False, mfrow=mfrow,
                         main=labels[reg], ylim=(ymin,ymax))
                new = False
                mfrow = []
            lineplot(X, [data[nplot-1, :, n] for n in range(data.shape[-1])],
                     output=pdf, new=new, last=True, main=labels[-1],
                     legend=snames, ylim=(ymin,ymax))
            if outf == 'archive':
                for n,sn in enumerate(snames):
                    _datf = self.temporary_path(fname=sn+"_data.txt")
                    with open(_datf,"w") as dff:
                        dff.write("\t".join([""]+[str(x) for x in X])+"\n")
                        for k,l in enumerate(labels):
                            dff.write("\t".join([l]+[str(x) for x in data[k, :, n]])+"\n")
                    tarfh.add(_datf,arcname=os.path.basename(_datf))
        else:
            raise ValueError("Mode not implemented: %s" % mode)
        if outf == 'archive':
            tarfh.add(pdf,arcname=os.path.basename(pdf))
            tarfh.close()
            self.new_file(tarname, 'data_archive')
        else:
            self.new_file(pdf, 'plot_features')
        return self.display_time()
