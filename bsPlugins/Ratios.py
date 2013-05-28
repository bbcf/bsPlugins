from bsPlugins import *
from bbcflib.bFlatMajor.stream import merge_scores
from bbcflib.btrack import track
from bbcflib import genrep
import os
import math


class RatiosForm(BaseForm):
    child = twd.HidingTableLayout()
    track1 = twf.FileField(
        label='File 1: ',
        help_text='Select the track with the numerators',)
    track2 = twf.FileField(
        label='File 2: ',
        help_text='Select the track with the denominators',)
    format = twf.SingleSelectField(label='Output format: ',
        prompt_text=None,
        options=["bedGraph","sql","wig","sga"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    assembly = twf.SingleSelectField(label='Assembly: ',
        prompt_text=None,
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'track1', 'type': 'track', 'required': True},
                 {'id': 'track2', 'type': 'track', 'required': True},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'assembly', 'type': 'assembly'},]
out_parameters = [{'id': 'ratios', 'type': 'track'}]


class RatiosPlugin(BasePlugin):
    """Divides the scores of the first track by the scores of the second,
and returns a single track with the ratios as new scores."""
    info = {
        'title': 'Score ratios',
        'description': __doc__,
        'path': ['Analysis', 'Ratios'],
        'output': RatiosForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def _divide(x1,x2):
        if x2 != 0:
            return float(x1)/x2
        else:
            return math.log(x1) # arbitrary, make it an option? "0"?

    def divide(self,**kw):
        assembly = genrep.Assembly(kw['assembly'])
        format = kw['format']
        track1 = kw['track1']
        track2 = kw['track2']
        name1,ext1 = os.path.splitext(os.path.basename(track1))
        name2,ext2 = os.path.splitext(os.path.basename(track2))
        t1 = track(track1, chrmeta=assembly.chrmeta)
        t2 = track(track2, chrmeta=assembly.chrmeta)
        output = self.temporary_path(fname='ratios_%s-%s'%(name1,name2)+format)
        tout = track(output, format=kw['format'], chrmeta=assembly.chrmeta)
        for chrom in chrmeta:
            s1 = t1.read(chrom)
            s2 = t2.read(chrom)
            s = merge_scores([s1,s2],method=_divide)
            tout.write(s)
        t1.close(); t2.close(); tout.close()
        self.new_file(output, 'ratios')
        return self.display_time()

