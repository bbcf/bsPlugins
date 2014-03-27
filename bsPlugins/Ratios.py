from bsPlugins import *
from bbcflib.gfminer.stream import merge_scores
from bbcflib.gfminer.stream import window_smoothing
from bbcflib.gfminer.figure import boxplot
from bbcflib.gfminer.common import unroll
from bbcflib.track import track
from bbcflib import genrep
from math import log
from numpy import array
import random

size_def = 1
pseudo_def = 0.5

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'numerator', 'type': 'track', 'required': True},
                 {'id': 'denominator', 'type': 'track', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'window_size', 'type': 'int'},
                 {'id': 'pseudo', 'type': 'float'},
                 {'id': 'log', 'type':'boolean', 'required':True},
                 {'id': 'distribution', 'type':'boolean', 'required':True}]
out_parameters = [{'id': 'ratios', 'type': 'track'}, {'id': 'boxplot', 'type': 'pdf'}]


class RatiosForm(BaseForm):
    numerator = twb.BsFileField(
        label='File 1: ',
        help_text='Select the track with the numerators',
        validator=twb.BsFileFieldValidator(required=True))
    denominator = twb.BsFileField(
        label='File 2: ',
        help_text='Select the track with the denominators',
        validator=twb.BsFileFieldValidator(required=True))
    assembly = twf.SingleSelectField(
        label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    format = twf.SingleSelectField(
        label='Output format ',
        prompt_text=None,
        options=["bedGraph","sql","wig","bigWig","sga"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    window_size = twf.TextField(label='Window size: ',
        validator=twc.IntValidator(),
        value=size_def,
        help_text='Size of the sliding window in bp (default: 1)')
    pseudo = twf.TextField(label='Pseudo-count: ',
        validator=twc.IntValidator(),
        value=pseudo_def,
        help_text='Value to be added to both signals (default: 0.5)')
    log = twf.CheckBox(label='Log ratios: ',
        value=False,
        help_text='Returns the log2 of the ratios')
    distribution = twf.CheckBox(label='Boxplot: ',
        value=False,
        help_text='Returns the boxplot of log2 of the ratios from a random sample of 10^4 values taken genomewise. If only a few regions are highly differentiated, then one expects the median to be zero')
    submit = twf.SubmitButton(id="submit", value="Submit")


class RatiosPlugin(BasePlugin):
    """Divides the mean scores of the first track by the mean scores of the second over a sliding window, and returns a single track with the ratios as new scores associated to the center of the window. Uses pseudo-counts (1/2), applies a log transform if `log` is True, and makes a boxplot of the log2 of the ratios if `Boxplot` is True."""
    info = {
        'title': 'Score ratios',
        'description': __doc__,
        'path': ['Signal', 'Ratios'],
        'output': RatiosForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def _divide(self,x):
        if len(x) > 1:
            if self.log:
                return log(self.pseudo+x[0],2)-log(self.pseudo+x[1],2)
            else:
                return (self.pseudo+x[0])/(self.pseudo+x[1])
        if self.log:
            return log(self.pseudo+x[0],2)+self.baseline
        return (self.pseudo+x[0])/self.pseudo

    def __call__(self,**kw):
        t1 = track(kw.get('numerator'),chrmeta=kw.get('assembly') or None)
        t2 = track(kw.get('denominator'),chrmeta=kw.get('assembly') or None)
        format = kw.get('format',t1.format)
        wsize = int(kw.get('window_size', size_def) or size_def)
        self.log = kw.get('log',False)
        if isinstance(self.log, basestring):
            self.log = (self.log.lower() in ['1', 'true', 't','on'])
        self.pseudo = float(kw.get('pseudo', pseudo_def) or pseudo_def)
        self.baseline = -log(self.pseudo,2)

        output = self.temporary_path(fname='ratios_%s-%s.%s'%(t1.name,t2.name,format))
        with track(output, chrmeta=t1.chrmeta, fields=t1.fields, info={'datatype': "qualitative"}) as tout:
            if wsize > 1:
                for chrom in t1.chrmeta.keys():
                    s1 = window_smoothing(t1.read(chrom),window_size=wsize,step_size=1,featurewise=False)
                    s2 = window_smoothing(t2.read(chrom),window_size=wsize,step_size=1,featurewise=False)
                    s = merge_scores([s1,s2],method=self._divide)
                    tout.write(s, chrom=chrom)
            else:
                for chrom in t1.chrmeta.keys():
                    s = merge_scores([t1.read(chrom),t2.read(chrom)],method=self._divide)
                    tout.write(s, chrom=chrom)
        self.new_file(output, 'ratios')

        self.distribution = kw.get('distribution',False)
        if isinstance(self.distribution, basestring):
            self.distribution = (self.distribution.lower() in ['1', 'true', 't','on'])
        if self.distribution:
            pdf = self.temporary_path(fname='boxplot.pdf')
            genome_length = 0
            for chrom in t1.chrmeta.keys():
                genome_length += t1.chrmeta[chrom]['length']
            positions = random.sample(range(genome_length), 10000)
            self.log = "True"
            p = -1; ratios = []; labels = []
            if wsize > 1:
                for chrom in t1.chrmeta.keys():
                    s1 = window_smoothing(t1.read(chrom),window_size=wsize,step_size=1,featurewise=False)
                    s2 = window_smoothing(t2.read(chrom),window_size=wsize,step_size=1,featurewise=False)
                    s = merge_scores([s1,s2],method=self._divide)
                    for l in unroll(s,regions=(0,t1.chrmeta[chrom]['length']),fields=['score']):
                        p += 1
                        if p == positions[0]:
                            ratios.append(l[0])
                            labels.append("G")
                            positions.remove(p)
            else:
                for chrom in t1.chrmeta.keys():
                    s = merge_scores([t1.read(chrom),t2.read(chrom)],method=self._divide)
                    for l in unroll(s,regions=(0,t1.chrmeta[chrom]['length']),fields=['score']):
                        p += 1
                        if p == positions[0]:
                            ratios.append(l[0])
                            labels.append("G")
                            positions.remove(p)
            boxplot(ratios,labels,output=pdf)
            self.new_file(pdf, 'boxplot')
            return self.display_time()

