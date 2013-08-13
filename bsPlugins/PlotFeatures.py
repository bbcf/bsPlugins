from bsPlugins import *
from bbcflib.gfminer.numeric import feature_matrix
from bbcflib.gfminer.figure import heatmap, lineplot
from bbcflib.track import track, FeatureStream
from numpy import vstack, concatenate, array, where
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
                 {'id': 'noclust', 'type':'boolean', 'required':True},
                 {'id': 'ymin', 'type': 'float'},
                 {'id': 'ymax', 'type': 'float'},
                 {'id': 'output', 'type': 'list', 'required': True}]
out_parameters = [{'id': 'plot_features', 'type': 'pdf'},
                  {'id':'data_archive', 'type':'file'}]


class PlotFeaturesForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signal: '
        signals = twb.BsFileField(label=' ',
                                  help_text='Select signal file (e.g. bedgraph)',
                                  validator=twb.BsFileFieldValidator())
    features = twb.BsFileField(label='Features: ',
                               help_text='Select a feature file (e.g. bed)',
                               validator=twb.BsFileFieldValidator())

    child = twd.HidingTableLayout()
    mode = twd.HidingSingleSelectField(label='Plot type: ',
                                       options=plot_types, mapping={0:['noclust']},
                                       prompt_text=None)
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
    noclust = twf.CheckBox(label='Do not cluster features: ',
                           value=False,
                           help_text='Keep features in the same order as input')
    ymin = twf.TextField(label='Minimal signal value: ', 
                         validator=twb.FloatValidator(),
                         help_text='Minimum value displayed in graphs (optional)')
    ymax = twf.TextField(label='Maximum signal value: ', 
                         validator=twb.FloatValidator(),
                         help_text='Maximum value displayed in graphs (optional)')
    output = twf.SingleSelectField(label='Output: ', options=output_list,
                                   prompt_text=None, help_text='Pdf only or data+pdf')
    submit = twf.SubmitButton(id="submit", value="Plot")


class PlotFeaturesPlugin(BasePlugin):
    """Plot several genomic signals on a selection of regions (features). 

`Heatmap` make one heatmap per signal file, `mosaic plot` creates one plot for each feature showing each signal as a separate line on the plot, `average lineplot` calculates the average of those line plots over all features."""
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

        def make_X_labels(X,start,end,strand,down,up):
            flen = end-start
            i0 = where(X == 0)[0][0]+1
            i1 = where(X == 1)[0][0]+1
            i2 = len(X)-i1
            istep  = 0.5/(i1-i0)
            if down < 1: down *= flen
            if up < 1: up *= flen
            Xup = (array(range(-i0,0))+.5)*up/i0
            Xb = (X[i0:i1]+istep)*flen
            Xdown = flen+(array(range(i2))+.5)*down/i2
            if strand is None or strand > 0: return start+concatenate([Xup,Xb,Xdown])
            else:                            return end-concatenate([Xup,Xb,Xdown])

        def add_name(_s):
            """Adds a name field to a stream using 'chr:start-end'."""
            ci = _s.fields.index('chr')
            si = _s.fields.index('start')
            ei = _s.fields.index('end')
            _f = _s.fields+['name']
            return FeatureStream((r+("%s:%i-%i"%(r[ci],r[si],r[ei]),) for r in _s), fields=_f)

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
        if kw.get("noclust") is not None: noclust = str(kw["noclust"]).lower() in ['1','true','t']
        else: noclust = False
        ymin = kw.get('ymin')
        ymax = kw.get('ymax')
        for chrom in features.chrmeta:
            if 'name' in features.fields: _fread = features.read(chrom)
            else: _fread = add_name(features.read(chrom))
            _l, _d = feature_matrix([s.read(chrom) for s in signals], _fread,
                                    segment=True, nbins=nbins, 
                                    upstream=upstr, downstream=downstr)
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
            if 'name' in features.fields: _fread = features.read(fields=['chr','start','end','name'])
            else: _fread = add_name(features.read(fields=['chr','start','end']))
            order = [where(labels == feat[3])[0][0] for feat in _fread]
            for n in range(data.shape[-1]-1):
                heatmap(data[order, :, n], output=pdf, new=new, last=False,
                        rows=labels[order], columns=X, main=snames[n],
                        orderRows=not(noclust), orderCols=False, 
                        ymin=ymin, ymax=ymax)
                new = False
            heatmap(data[order, :, -1], output=pdf, new=new, last=True,
                    rows=labels[order],  columns=X, main=snames[-1],
                    orderRows=not(noclust), orderCols=False, 
                    ymin=ymin, ymax=ymax)
            if outf == 'archive':
                for n,sn in enumerate(snames):
                    _datf = self.temporary_path(fname=sn+"_data.txt")
                    with open(_datf,"w") as dff:
                        dff.write("\t".join([""]+[str(x) for x in X])+"\n")
                        for k in order:
                            dff.write("\t".join([labels[k]]+[str(x) for x in data[k, :, n]])+"\n")
                    tarfh.add(_datf,arcname=os.path.basename(_datf))
        elif mode in plot_types[1]: #average lineplot
            Y = data.mean(axis=0)
            if ymin is None: ymin = min([x.min() for x in Y]+[0])
            if ymax is None: ymax = max([x.max() for x in Y])
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
            mfrow = [4,3]
            nplot = min(data.shape[0], max_pages*mfrow[0]*mfrow[1])
            if ymin is None: ymin = min([data.min(),0])
            if ymax is None: ymax = data.max()
            _f = ['chr','start','end']
            _si = None
            if 'strand' in features.fields: 
                _f.append('strand')
                _si = 3
            if 'name' in features.fields: _fread = features.read(fields=_f+['name'])
            else: _fread = add_name(features.read(fields=_f))
            order = []
            for nf,feat in enumerate(_fread):
                reg = where(labels == feat[-1])[0][0]
                order.append(reg)
                X1 = make_X_labels(X, feat[1], feat[2], feat[_si] if _si else None, downstr[0], upstr[0])
                Y = [data[reg, :, n] for n in range(data.shape[-1])]
                if nf == 0:
                    lineplot(X1, Y,  output=pdf, new=True, last=False, mfrow=mfrow,
                             main=labels[reg], ylim=(ymin,ymax))
                elif nf < nplot-1:
                    lineplot(X1, Y, output=pdf, new=False, last=False, 
                             main=labels[reg], ylim=(ymin,ymax))
                else:
                    lineplot(X1, Y, output=pdf, new=False, last=True, legend=snames, 
                             main=labels[reg], ylim=(ymin,ymax))
                    break
            if outf == 'archive':
                for n,sn in enumerate(snames):
                    _datf = self.temporary_path(fname=sn+"_data.txt")
                    with open(_datf,"w") as dff:
                        dff.write("\t".join([""]+[str(x) for x in X])+"\n")
                        for k in order:
                            dff.write("\t".join([labels[k]]+[str(x) for x in data[k, :, n]])+"\n")
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
