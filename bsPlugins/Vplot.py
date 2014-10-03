from bsPlugins import *
from bbcflib.gfminer.figure import smoothScatter
from bbcflib.track import track
from numpy import asarray, concatenate, mean

nbin_x_def = 500
nbin_y_def = 500
bandwidth_x_def = 0.1
bandwidth_y_def = 0.1

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
out_parameters = [{'id': 'Vplot', 'type': 'png'},
                  {'id': 'Vplots_archive', 'type': 'file'}]


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
            left_right = (left_right.lower() in ['1', 'true', 't', 'on'])
        linscale = kw.get('linear',False)
        if isinstance(linscale, basestring):
            linscale = (linscale.lower() in ['1', 'true', 't', 'on'])
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
        new = True
        extra_window = 1000
        strand = 1
        pnglist = []
        for bam_nb, bam in enumerate(bamfiles):
            XL = None; XR = None; Y = None
            for region_nb, region in enumerate(features.read()):
                if strandi > -1: strand = region[strandi]
                chrom,start,end = region[:3]
                _XL = []; _XR = []; _Y = []
                for read in bam.fetch(chrom,max(0,start-extra_window),end+extra_window):
                    if read.is_proper_pair and read.isize>0 and not read.is_reverse:
                        _rs = read.isize
                        if strand < 0: rpos = end-read.pos-_rs
                        else:          rpos = read.pos-start
                        if rpos < -_rs: continue
                        if rpos >= end-start+_rs: break
                        _Y.append(_rs)
                        if left_right:
                            _XL.append(rpos)
                            _XR.append(rpos+_rs)
                        else:
                            _XR.append(rpos+_rs/2)
                if Y is None:
                    if left_right: XL = asarray(_XL)
                    XR = asarray(_XR)
                    Y = asarray(_Y)
                else:
                    if left_right: XL = concatenate((XL,asarray(_XL)))
                    XR = concatenate((XR,asarray(_XR)))
                    Y = concatenate((Y,asarray(_Y)))
            ylims = (int(kw.get('ymin') or ymin_def), int(kw.get('ymax') or max(Y)))
            xlims = (0,end-start)
            colrs = ["white","blue","red"]
            mlabel = bam.name
            png = self.temporary_path(fname='Vplot_%s.png'%mlabel)
            if left_right:
                smoothScatter( XL, Y, output=png, new=True, last=False, main=mlabel+" left fragment end",
                               xlab=xlab, ylab=ylab, xlim=xlims, ylim=ylims, log=log, color=["white","red"],
                               mfrow=[2,1], nbin=nbin, bandwidth=bandwidth )
                colrs = ["white","blue"]
                mlabel += " right fragment end"
                new = False
            smoothScatter( XR, Y, output=png, new=new, last=True, main=mlabel,
                           xlab=xlab, ylab=ylab, xlim=xlims, ylim=ylims, log=log, color=colrs,
                           nbin=nbin, bandwidth=bandwidth )
            new = True
            pnglist.append(png)
        if len(pnglist) > 1:
            tar_png_name = self.temporary_path('Vplots.tgz')
            tar_png = tarfile.open(tar_png_name, "w:gz")
            [tar_png.add(f,arcname=os.path.basename(f)) for f in pnglist]
            tar_png.close()
            self.new_file(tar_png_name, 'Vplots_archive')
        elif len(pnglist) == 1:
            self.new_file(pnglist[0], 'Vplot')
        return self.display_time()
