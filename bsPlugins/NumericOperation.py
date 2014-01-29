from bsPlugins import *
from bbcflib import genrep
from bbcflib.track import track
from bbcflib.gfminer.common import score_threshold, apply
from math import log10, sqrt, log
import os, tarfile

class NumericOperationForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signals: '
        track = twb.BsFileField(label=' ',
                                help_text='Select files (e.g. bedgraph)',
                                validator=twb.BsFileFieldValidator(required=True))
    function =  twf.SingleSelectField(label='Operation: ',
                                      prompt_text=None,
                                      options=["log2","log10","sqrt"],
                                      help_text='Select a function')
    format = twf.SingleSelectField(label='Output format: ',
                                   options=["sql","bedgraph","bigwig","wig"],
                                   validator=twc.Validator(required=False),
                                   help_text='Output file(s) format, by default: same format as input file(s) format(s)')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'track', 'type': 'track', 'required': True, 'multiple':'SigMulti'},
                {'id': 'function', 'type': 'function'},
                {'id': 'format', 'type': 'format'}]
out_parameters = [{'id': 'converted_track_tar', 'type': 'file'},
                  {'id': 'converted_track', 'type': 'track'}]

def log2(x):
    return log(x,2)

class NumericOperationPlugin(BasePlugin):
    """Apply a numeric transformation to the track scores - such as logarithm or square root."""
    info = {
        'title': 'Numeric Operation',
        'description': __doc__,
        'path': ['Signal', 'Numeric Operation'],
        'output': NumericOperationForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        func = kw.get('function',"log2")
        l_track = kw.get('SigMulti', {}).get('track',[])
        if not isinstance(l_track, list): l_track = [l_track]
        outall = []
        for tname in l_track :
            tinput = track(tname)
            if 'score' not in tinput.fields: continue
            format = kw.get('format',tinput.format)
            out_name = tinput.name+'_'+func+'.'+format
            outtemp = self.temporary_path(out_name)
            out_track = track(outtemp,chrmeta=tinput.chrmeta)
            filtered = score_threshold(tinput, strict=(func[:3] == "log"))
            out_track.write(apply(filtered,'score',eval(func)), mode='write')
            out_track.close()
            outall.append(outtemp)
            tinput.close()
        if len(outall) == 1:
            self.new_file(outall[0], 'converted_track')
        elif len(outall) > 1:
            tar_name = self.temporary_path(fname="numeric_operation_out.tgz")
            tar = tarfile.open(tar_name, "w:gz")
            [tar.add(f,arcname=os.path.basename(f)) for f in outall]
            tar.close()
            self.new_file(tar_name, 'converted_track_tar')
        return self.display_time()

# nosetests --logging-filter=-tw2 test_NumericOperation.py
