from bsPlugins import *
from bbcflib.gfminer.common import unroll
from bbcflib.gfminer.figure import Vplot
from bbcflib.track import track
from math import floor, ceil

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
        nb_plots = len(bamfiles)
        features = track(kw.get('features'), chrmeta=bamfiles[0].chrmeta)
        scale = kw.get('linear')
        pdf = self.temporary_path(fname='Vplot.pdf')
        bam_nb = 0
        for bam in bamfiles:
            bam_nb += 1
            bam_name = bam.name.split(".")[0]
            list_regions = features.read()
            scores = []; nb_regions = 0; X = []
            for region in list_regions:
                density = unroll(bam.PE_fragment_size(region,midpoint="midpoint"),regions=(region[1],region[2]))
                for score in density:
                    scores.append(score[0])
                nb_regions += 1
                if nb_regions == 1:
                    interval_length = len(scores)
                    xmin = -int(floor(interval_length/2.0))
                    xmax = int(ceil(interval_length/2.0))
                try:
                    strand = region[5]
                    if strand == 1: X = X+range(xmin,xmax)
                    else: X = X+list(reversed(range(xmin,xmax)))
                except:
                    X = X+range(xmin,xmax)
            Y = ["NA" if x==0 else x for x in scores]
            xmin = min(X); xmax = max(X); ymax = max(scores)
            if bam_nb == 1:
                if scale:
                    Vplot(X,Y,output=pdf,new=True,last=False,main=bam_name,xlab="Distance to the region centers",ylab="Mean non-zero fragment length",xlim=(xmin,xmax),ylim=(0,ymax))
                else:
                    Vplot(X,Y,output=pdf,new=True,last=False,main=bam_name,xlab="Distance to the region centers",ylab="Mean non-zero fragment length",xlim=(xmin,xmax),ylim=(50,ymax),log="y")
            elif bam_nb < nb_plots:
                if scale:
                    Vplot(X,Y,output=pdf,new=False,last=False,main=bam_name,xlab="Distance to the region centers",ylab="Mean non-zero fragment length",xlim=(xmin,xmax),ylim=(0,ymax))
                else:
                    Vplot(X,Y,output=pdf,new=False,last=False,main=bam_name,xlab="Distance to the region centers",ylab="Mean non-zero fragment length",xlim=(xmin,xmax),ylim=(50,ymax),log="y")
            else:
                if scale:
                    Vplot(X,Y,output=pdf,new=False,last=True,main=bam_name,xlab="Distance to the region centers",ylab="Mean non-zero fragment length",xlim=(xmin,xmax),ylim=(0,ymax))
                else:
                    Vplot(X,Y,output=pdf,new=False,last=True,main=bam_name,xlab="Distance to the region centers",ylab="Mean non-zero fragment length",xlim=(xmin,xmax),ylim=(50,ymax),log="y")
        self.new_file(pdf, 'Vplot')
        return self.display_time()
