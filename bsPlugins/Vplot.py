from bsPlugins import *
from bbcflib.gfminer.figure import Vplot
from bbcflib.track import track

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'bamfiles', 'type': 'bam', 'required': True, 'multiple': 'BamMulti'},
                 {'id': 'features', 'type': 'track'},
                 {'id': 'linear', 'type': 'boolean'}]
out_parameters = [{'id': 'Vplot', 'type': 'pdf'}]


class VplotForm(BaseForm):
    class TrackMulti(twb.BsMultiple):
        label='Paired-end BAM files: '
        bamfiles = twb.BsFileField(label=' ',
                                  validator=twb.BsFileFieldValidator(required=True))
    features = twb.BsFileField(label='Features: ',
                               help_text='Select a feature file (e.g. bed) in which all regions have the same length',
                               validator=twb.BsFileFieldValidator(required=True))
    linear = twf.CheckBox(label='Linear scale: ',
                             value=False,
                             help_text='Plot the mean fragment length in linear scale (default: log scale)')
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
        linscale = kw.get('linear',False)
        if isinstance(linscale, basestring):
            linscale = (linscale.lower() in ['1', 'true', 't','on'])
        if linscale: 
            ymin = 0
            log = ''
        else: 
            ymin = 50
            log = 'y'
        xlab = "Position in window [bp]"
        ylab = "Fragment size [bp]"
        pdf = self.temporary_path(fname='Vplot.pdf')
        new = True
        last = False
        for bam_nb, bam in enumerate(bamfiles):
            if bam_nb == len(bamfiles)-1: last = True
            list_regions = features.read()
            Y = []
            X = []
            for region_nb, region in enumerate(list_regions):
                if strandi > -1:
                    strand = region[strandi]
                else:
                    strand = 1
                for _s in bam.PE_fragment_size(region,midpoint="midpoint"):
                    for pos in range(_s[1],_s[2]):
                        if pos < region[1]: continue
                        if pos >= region[2]: break
                        Y.append(_s[3])
                        if strand < 0:
                            X.append(region[2]-pos-1)
                        else:
                            X.append(pos-region[1])
            ymax = max(Y)
            Vplot( X, Y, output=pdf, new=new, last=last, main=bam.name, 
                   xlab=xlab, ylab=ylab, ylim=(ymin,ymax), log=log )
            new = False
        self.new_file(pdf, 'Vplot')
        return self.display_time()
