from bsPlugins import *
from bbcflib.gfminer.common import unroll, add_name_field
from bbcflib.gfminer.numeric import feature_matrix
#from bbcflib.gfminer.figure import Vplot, heatmap, lineplot
from bbcflib.track import track, FeatureStream
from numpy import vstack, concatenate, array, where
import os, tarfile, pysam

dot_size_def = 4

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'bamfiles', 'type': 'bam', 'required': True, 'multiple': 'BamMulti'},
                 {'id': 'features', 'type': 'track'},
                 {'id': 'dot_size', 'type': 'int'},
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
    dot_size = twf.TextField(label='Dot size: ',
                             validator=twc.IntValidator(required=False),
                             value=dot_size_def,
                             help_text='Size of the dots used to draw the Vplot (default: 4)')
    linear = twf.CheckBox(label='Linear scale: ',
                             value=False,
                             help_text='Plot the mean fragment length in linear scale (default: log scale)')
    submit = twf.SubmitButton(id="submit", value="Plot")


class VplotPlugin(BasePlugin):
    """Draw the dotplot of the mean fragment lengths (associated to the fragment middle points) corresponding to a paired-end BAM file and a set of regions (features)."""
    info = {
        'title': 'Make Vplot for a paired-end BAM file in a selection of regions',
        'description': __doc__,
        'path': ['Graphics', 'Plot features'],
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
        dot_size = int(kw.get('dot_size'))
        if dot_size <= 0:
            dot_size = dot_size_def
        scale = kw.get('linear')
        if scale:
            scale = "linear"
        else:
            scale = "log"
        pdf = self.temporary_path(fname='Vplot.pdf')
        #robjects.r('pdf("%s",paper="a4",height=11,width=8)' %pdf)
        for bam in bamfiles:
            bam_name = bam.name.split(".")[0]
            list_regions = features.read()
            data = []
            for region in list_regions:
                density = unroll(bam.PE_fragment_size(region,midpoint="midpoint"),regions=(region[1],region[2]),fields=["score"])
                scores = []
                print "OK1"
                for score in density:
                    print "OK2"
                    scores.append(score[2])
                data.append(scores)
                print data
            # data[i][j] = mean frag length in region i at position j
            Vplot(data,output=pdf,scale=scale,title=bam_name)
        #robjects.r('dev.off()')
        self.new_file(pdf, 'Vplot')
        return self.display_time()
