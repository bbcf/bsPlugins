from bsPlugins import *
from bbcflib.gfminer.stream import merge_scores, window_smoothing
from bbcflib.gfminer.figure import density_boxplot
from bbcflib.gfminer.common import unroll
from bbcflib.track import track, FeatureStream
from bbcflib import genrep
from math import log
from numpy.random import poisson
from numpy import median

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
        label='Numerator: ',
        help_text='Select the track with the numerators',
        validator=twb.BsFileFieldValidator(required=True))
    denominator = twb.BsFileField(
        label='Denominator: ',
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
        validator=twb.FloatValidator(min=0,max=1000),
        value=pseudo_def,
        help_text='Value to be added to both signals (default: 0.5)')
    log = twf.CheckBox(label='Log ratios: ',
        value=False,
        help_text='Computes the log2 of the ratios')
    distribution = twf.CheckBox(label='Plot distribution: ',
        value=False,
        help_text='Creates a graphical representation of the distributions of the ratios based on a sample of genomic regions')
    submit = twf.SubmitButton(id="submit", value="Submit")


class RatiosPlugin(BasePlugin):
    """Divides the mean scores of the first track by the mean scores of the second over a sliding window, and returns a single track with the ratios as new scores associated to the center of the window. Uses pseudo-counts (1/2), applies a log transform if `Log ratios` is True, and makes a plot of the log2 of the ratios if `Plot distribution` is True."""

    sample_length = 200
    sample_num = 50000

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

    def _sample_stream(self, stream, limit):
        start = 0
        end = -1
        scores = []
        ist = stream.fields.index('start')
        ien = stream.fields.index('end')
        isc = stream.fields.index('score')
        for x in stream:
            yield x
            if start > limit: continue
            if x[ist] >= end:
                self.ratios.extend(scores)
                scores = [0]*self.sample_length
                start += self.shifts[0]
                end = start+self.sample_length
                if end <= limit:
                    self.shifts.pop(0)
                else:
                    start = limit+1
            elif x[ien] > start:
                _s = max(0,x[ist]-start)
                _e = min(self.sample_length,x[ien]-start)
                scores[_s:_e] = [x[isc]]*(_e-_s)


    def __call__(self,**kw):
        assembly = kw.get('assembly') or 'guess'
        t1 = track(kw['numerator'],chrmeta=assembly)
        t2 = track(kw['denominator'],chrmeta=assembly)
        format = kw.get('format') or t1.format
        wsize = int(kw.get('window_size') or size_def)
        self.log = kw.get('log',False)
        if isinstance(self.log, basestring):
            self.log = (self.log.lower() in ['1', 'true', 't','on'])
        try:
            self.pseudo = float(kw.get('pseudo'))
        except:
            self.pseudo = pseudo_def
        self.baseline = -log(self.pseudo,2)
        distribution = kw.get('distribution',False)
        if isinstance(distribution, basestring):
            distribution = (distribution.lower() in ['1', 'true', 't','on'])
        if distribution:
            genome_length = sum((v['length'] for v in t1.chrmeta.values()))
            self.shifts = list(poisson(float(genome_length)/float(self.sample_num),
                                       self.sample_num))
            self.ratios = []

        output = self.temporary_path(fname='ratios_%s-%s.%s'%(t1.name,t2.name,format))
        with track(output, chrmeta=t1.chrmeta, fields=t1.fields,
                   info={'datatype': 'quantitative',
                         'log': self.log,
                         'pseudocounts': self.pseudo,
                         'window_size': wsize}) as tout:
            for chrom,vchr in t1.chrmeta.iteritems():
                if wsize > 1:
                    s1 = window_smoothing(t1.read(chrom),window_size=wsize,step_size=1,featurewise=False)
                    s2 = window_smoothing(t2.read(chrom),window_size=wsize,step_size=1,featurewise=False)
                else:
                    s1 = t1.read(chrom)
                    s2 = t2.read(chrom)
                s3 = merge_scores([s1,s2],method=self._divide)
                if distribution:
                    s3 = FeatureStream(self._sample_stream(s3,vchr['length']),fields=s3.fields)
                tout.write(s3, chrom=chrom, clip=True)
        self.new_file(output, 'ratios')

        if distribution:
            pdf = self.temporary_path(fname='%s-%s_ratios_distribution.pdf'%(t1.name,t2.name))
            density_boxplot(self.ratios,output=pdf,name=t1.name+"/"+t2.name+" (median="+str(round(median(self.ratios),2))+")")
            self.new_file(pdf, 'boxplot')
        return self.display_time()
