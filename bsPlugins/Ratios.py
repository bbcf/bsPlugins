from bsPlugins import *
from bbcflib.gfminer.stream import merge_scores
from bbcflib.track import track
from math import log

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'numerator', 'type': 'track', 'required': True},
                 {'id': 'denominator', 'type': 'track', 'required': True},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'log', 'type':'boolean', 'required':True}]
out_parameters = [{'id': 'ratios', 'type': 'track'}]


class RatiosForm(BaseForm):
    numerator = twb.BsFileField(
        label='File 1: ',
        help_text='Select the track with the numerators',
        validator=twb.BsFileFieldValidator(required=True))
    denominator = twb.BsFileField(
        label='File 2: ',
        help_text='Select the track with the denominators',
        validator=twb.BsFileFieldValidator(required=True))
    format = twf.SingleSelectField(
        label='Output format ',
        prompt_text=None,
        options=["bedGraph","sql","wig","bigWig","sga"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    log = twf.CheckBox(label='Log ratios: ',
                       value=False,
                       help_text='Returns the log2 of the ratios')
    submit = twf.SubmitButton(id="submit", value="Submit")


class RatiosPlugin(BasePlugin):
    """Divides the scores of the first track by the scores of the second, and returns a single track with the ratios as new scores. Uses pseudo-counts (1/2) and applies a log transform is `log` is True."""
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
        numerator = kw['numerator']
        denominator = kw['denominator']
        t1 = track(numerator)
        t2 = track(denominator)
        format = kw.get('format',t1.format)
        self.log = kw.get('log',False)
        self.pseudo = kw.get('pseudo')
        if self.pseudo is not None:
            self.pseudo = float(self.pseudo)
            self.baseline = -log(self.pseudo,2)
        else:
            self.pseudo = 0.5
            self.baseline = 1.0
        
        output = self.temporary_path(fname='ratios_%s-%s.%s'%(t1.name,t2.name,format))
        with track(output, chrmeta=t2.chrmeta, fields=t1.fields) as tout:
            for chrom in t2.chrmeta.keys():
                s = merge_scores([t1.read(chrom),t2.read(chrom)],method=self._divide)
                tout.write(s, chrom=chrom)
        self.new_file(output, 'ratios')
        return self.display_time()

