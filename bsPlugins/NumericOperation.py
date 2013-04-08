from bsPlugins import *
from bbcflib import genrep
from bbcflib.btrack import track
from bbcflib.bFlatMajor import common
import math
import tw2.forms as twf


class NumericOperationForm(BaseForm):
    class SigMulti(Multi):
        label='Signals: '
        track = twf.FileField(label=' ',
        help_text='Select files (e.g. bedgraph)',
        validator=twf.FileValidator(required=True))
    function =  twf.SingleSelectField(label='Operation: ',
        options=["log2","log10","sqrt"],
        validator=twc.Validator(required=True),
        help_text='Select a method, by default: log2')
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bedgraph","bigwig","wig"],
        validator=twc.Validator(required=False),
        help_text='Output file(s) format, by default: same format as input file(s) format(s)')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [{'id': 'track', 'type': 'track', 'required': True},
                {'id': 'assembly', 'type': 'assembly', 'required': True},
                {'id': 'function', 'type': 'function', 'required': True},
                {'id': 'format', 'type': 'format'}]
out_parameters = [{'id': 'output', 'type': 'file'}]


class NumericOperationPlugin(BasePlugin):
    description = """Apply a numeric transformation to the track scores - such as logarithm or square root."""
    info = {
        'title': 'Numeric Operation',
        'description': description,
        'path': ['Signal', 'Numeric Operation'],
        'output': NumericOperationForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        def method(x):    # the function that is applied to the scores
            if kw['function']=="log2":
                return math.log(x,2) ;
            elif kw['function']=="log10":
                return math.log10(x) ;
            elif kw['function']=="sqrt":
                return math.sqrt(x) ;
            else:
                return math.log(x) ;

        assembly = genrep.Assembly(kw.get('assembly'))
        for i in range(len(kw.get('track'))) :
            tname = kw['track'][i]
            tinput = track(tname, chrmeta=kw.get('assembly'))
            (filepath, filename) = os.path.split(tname)
            (shortname, extension) = os.path.splitext(filename)
            modif = kw['function']
            if not isinstance(tinput, (tuple,list)):
                tinput = [tinput]
            for t in tinput:
                if  "score" in t.fields:
                    if kw['format']=="":
                        out_name = shortname+'_'+modif+str(extension)
                    else:
                        out_name = shortname+'_'+modif +'.'+kw['format']
                    output_name = self.temporary_path(out_name)
                    out_track = track(output_name,chrmeta=assembly.chrmeta)
                    out_track.write(common.apply(t.read(),'score',method), mode='write')
                    out_track.close()
                    t.close()
        return self.display_time()
