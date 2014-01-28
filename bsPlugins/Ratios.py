from bsPlugins import *
from bbcflib.gfminer.stream import merge_scores
from bbcflib.track import track

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
        options=["bedGraph","sql","wig","sga"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'numerator', 'type': 'track', 'required': True},
                 {'id': 'denominator', 'type': 'track', 'required': True},
                 {'id': 'format', 'type': 'text'}]
out_parameters = [{'id': 'ratios', 'type': 'track'}]


def _divide(x):
    if abs(x[0]) < 1e-5:
        return 0.0
    if len(x) == 2:
        if abs(x[1]) > 1e-5:
            return float(x[0])/x[1]
    return 1e10 # or 0?

class RatiosPlugin(BasePlugin):
    """Divides the scores of the first track by the scores of the second,
and returns a single track with the ratios as new scores."""
    info = {
        'title': 'Score ratios',
        'description': __doc__,
        'path': ['Signal', 'Ratios'],
        'output': RatiosForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self,**kw):
        numerator = kw['numerator']
        denominator = kw['denominator']
        t1 = track(numerator)
        t2 = track(denominator)
        name1 = t1.name
        name2 = t2.name
        format = kw.get('format',t2.format)
        output = self.temporary_path(fname='ratios_%s-%s.%s'%(name1,name2,format))
        with track(output, chrmeta=t2.chrmeta, fields=t2.fields) as tout:
            for chrom in t2.chrmeta.keys():
                s = merge_scores([t1.read(chrom),t2.read(chrom)],method=_divide)
                tout.write(s, chrom=chrom)
        self.new_file(output, 'ratios')
        return self.display_time()

