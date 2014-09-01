from bsPlugins import *
from bbcflib.gfminer.figure import Vplot, Vplot2
from bbcflib.track import track
from numpy import asarray, concatenate, mean

nbin_x_def = 500; nbin_y_def = 500
bandwidth_x_def = 0.1; bandwidth_y_def = 0.1

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'bamfiles', 'type': 'bam', 'required': True, 'multiple': 'BamMulti'},
                 {'id': 'features', 'type': 'track'},
                 {'id': 'left_right', 'type': 'boolean'},
                 {'id': 'linear', 'type': 'boolean'},
                 {'id': 'nbins_x', 'type': 'int'},
                 {'id': 'nbins_y', 'type': 'int'},
                 {'id': 'bandwidth_x', 'type': 'float'},
                 {'id': 'bandwidth_y', 'type': 'float'},
                 {'id': 'ymin', 'type': 'int'},
                 {'id': 'ymax', 'type': 'int'}]
out_parameters = [{'id': 'Vplot', 'type': 'pdf'}]


class VplotForm(BaseForm):
    class BamMulti(twb.BsMultiple):
        label='Paired-end BAM files: '
        bamfiles = twb.BsFileField(label=' ',
                                   validator=twb.BsFileFieldValidator(required=True))
    features = twb.BsFileField(label='Features: ',
                                help_text='Select a feature file (e.g. bed) in which all regions have the same length',
                                validator=twb.BsFileFieldValidator(required=True))
    left_right = twf.CheckBox(label='Left-right: ',
                                value=False,
                                help_text='Plot the mean fragment length associated to left and right fragment ends (default: false)')
    linear = twf.CheckBox(label='Linear scale: ',
                                value=False,
                                help_text='Plot the mean fragment length in linear scale (default: log scale)')
    nbin_x = twf.TextField(label='Number of bins along x axis: ',
                                validator=twc.IntValidator(required=False),
                                value=nbin_x_def,
                                help_text='Number of equally spaced grid points for the density estimation (default: 500)')
    nbin_y = twf.TextField(label='Number of bins along y axis: ',
                                validator=twc.IntValidator(required=False),
                                value=nbin_y_def,
                                help_text='Number of equally spaced grid points for the density estimation (default: 500)')
    bandwidth_x = twf.TextField(label='Smoothing bandwidth along x axis: ',
                                validator=twb.FloatValidator(min=0,max=1000),
                                value=bandwidth_x_def,
                                help_text='The smoothing bandwidth must be between 0 and 1000 (default: 0.1)')
    bandwidth_y = twf.TextField(label='Smoothing bandwidth along y axis: ',
                                validator=twb.FloatValidator(min=0,max=1000),
                                value=bandwidth_y_def,
                                help_text='The smoothing bandwidth must be between 0 and 1000 (default: 0.1)')
    ymin = twf.TextField(label='Minimum y value: ',
                                validator=twc.IntValidator(required=False),
                                help_text='The default values: ymin=0 in lin scale and ymin=50 in log scale')
    ymax = twf.TextField(label='Maximum y value: ',
                                validator=twc.IntValidator(required=False),
                                help_text='The default value: ymax=maximum fragment length in the selected regions')
    submit = twf.SubmitButton(id="submit", value="Plot")


class VplotPlugin(BasePlugin):
    """Draw the dotplot of the mean fragment lengths (associated to the fragment middle points) corresponding to a paired-end BAM file and a set of regions (features)."""
    info = {
        'title': 'Make Vplot for a paired-end BAM file in a selection of regions',
        'description': __doc__,
        'path': ['Graphics', 'Paired-end Vplots'],
        'output': VplotForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        bamfiles = kw.get('BamMulti',{}).get('bamfiles',[])
        if not isinstance(bamfiles, (tuple,list)): bamfiles = [bamfiles]
        bamfiles = [track(bam) for bam in bamfiles]
        features = track(kw.get('features'), chrmeta=bamfiles[0].chrmeta)
        if 'strand' in features.fields:
            strandi = features.fields.index('strand')
        else:
            strandi = -1
        left_right = kw.get('left_right',False)
        if isinstance(left_right, basestring):
            left_right = (left_right.lower() in ['1', 'true', 't','on'])
        linscale = kw.get('linear',False)
        if isinstance(linscale, basestring):
            linscale = (linscale.lower() in ['1', 'true', 't','on'])
        if linscale:
            ymin_def = 0
            log = ''
        else:
            ymin_def = 50
            log = 'y'
        nbin = (int(kw.get('nbin_x') or nbin_x_def),
                int(kw.get('nbin_y') or nbin_y_def))
        bandwidth = (float(kw.get('bandwidth_x') or bandwidth_x_def),
                     float(kw.get('bandwidth_y') or bandwidth_y_def))
        xlab = "Position in window [bp]"
        ylab = "Fragment size [bp]"
        pdf = self.temporary_path(fname='Vplot.pdf')
        new = True
        last = False
        extra_window = 1000
        for bam_nb, bam in enumerate(bamfiles):
            if bam_nb == len(bamfiles)-1: last = True
            if left_right:
                XL = None; YL = None; XR = None; YR = None
                for region_nb, region in enumerate(features.read()):
                    if strandi > -1:
                        strand = region[strandi]
                    else:
                        strand = 1
                    _XL = []; _YL = []; _XR = []; _YR = []
                    region_extended = (region[0],max(0,region[1]-extra_window))+region[2:]
                    for _s in bam.PE_fragment_size(region_extended,end="left"):
                        for pos in range(_s[1],_s[2]):
                            if pos < region[1]: continue
                            if pos >= region[2]: break
                            _YL.append(int(_s[3]))
                            if strand < 0:
                                _XL.append(region[2]-pos-1)
                            else:
                                _XL.append(pos-region[1])
                    if XL is None:
                        XL = asarray(_XL); YL = asarray(_YL)
                    else:
                        XL = concatenate((XL,asarray(_XL))); YL = concatenate((YL,asarray(_YL)))
                    for _s in bam.PE_fragment_size(region_extended,end="right"):
                        for pos in range(_s[1],_s[2]):
                            if pos < region[1]: continue
                            if pos >= region[2]: break
                            _YR.append(int(_s[3]))
                            if strand < 0:
                                _XR.append(region[2]-pos-1)
                            else:
                                _XR.append(pos-region[1])
                    if XR is None:
                        XR = asarray(_XR); YR = asarray(_YR)
                    else:
                        XR = concatenate((XR,asarray(_XR))); YR = concatenate((YR,asarray(_YR)))
                ylims = (int(kw.get('ymin') or ymin_def), int(kw.get('ymax') or max(concatenate((YL,YR)))))
                Vplot2( XL, YL, XR, YR, output=pdf, new=new, last=last, main=bam.name,
                       xlab=xlab, ylab=ylab, ylim=ylims, log=log, nbin=nbin,
                       bandwidth=bandwidth )
            else:
                X = None; Y = None
                for region_nb, region in enumerate(features.read()):
                    if strandi > -1:
                        strand = region[strandi]
                    else:
                        strand = 1
                    _X = []; _Y = []
                    region_extended = (region[0],max(0,region[1]-extra_window))+region[2:]
                    for _s in bam.PE_fragment_size(region_extended,midpoint=True):
                        for pos in range(_s[1],_s[2]):
                            if pos < region[1]: continue
                            if pos >= region[2]: break
                            _Y.append(int(_s[3]))
                            if strand < 0:
                                _X.append(region[2]-pos-1)
                            else:
                                _X.append(pos-region[1])
                    if X is None:
                        X = asarray(_X); Y = asarray(_Y)
                    else:
                        X = concatenate((X,asarray(_X))); Y = concatenate((Y,asarray(_Y)))
                ylims = (int(kw.get('ymin') or ymin_def), int(kw.get('ymax') or max(Y)))
                Vplot( X, Y, output=pdf, new=new, last=last, main=bam.name,
                       xlab=xlab, ylab=ylab, ylim=ylims, log=log, nbin=nbin,
                       bandwidth=bandwidth )
            new = False
        self.new_file(pdf, 'Vplot')
        return self.display_time()
