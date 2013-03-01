from bsPlugins import *
from bbcflib.bFlatMajor.stream import combine
from bbcflib.btrack import track
from bbcflib import genrep


class CombineForm(BaseForm):
    child = twd.HidingTableLayout()
    class SigMulti(Multi):
        signals = twf.FileField(label='Signals: ',
                                help_text='Select files to combine',
                                validator=twf.FileValidator(required=True))
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bed","sga"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Quantify")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'signals', 'type': 'track', 'multiple': True, 'required': True},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'assembly', 'type': 'assembly'},
                ]
out_parameters = [{'id': 'combined', 'type': 'track'}]


def _get_chrmeta(**kw):
    chrmeta = "guess"
    assembly_id = kw.get('assembly')
    if assembly_id:
        assembly = genrep.Assembly(assembly_id)
        chrmeta = assembly.chrmeta
    return chrmeta

def _combine(func,output,**kw):
    chrmeta = _get_chrmeta(**kw)
    format = kw.get('format','sql')
    output += format
    signals = kw.get('signals', [])
    if not isinstance(signals, list): signals = [signals]
    signals = [track(sig, chrmeta=chrmeta) for sig in signals]
    tout = track(output, chrmeta=chrmeta, info={'datatype':'qualitative'})
    for chrom in chrmeta:
        trackList = [sig.read(chrom) for sig in signals]
        res = combine(trackList, fn=func)
        tout.write(res, chrom=chrom, clip=True)
    tout.close()
    return output


class IntersectPlugin(OperationPlugin):
    info = {
        'title': 'Intersection of a set of tracks',
        'description': 'Returns a new track with only regions covered in every input track.',
        'path': ['Features', 'Intersect'],
        'output': CombineForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        func = all
        output = self.temporary_path(fname='combined.')
        output = _combine(func,output,**kw)
        self.new_file(output, 'combined')
        return self.display_time()

class UnionPlugin(OperationPlugin):
    info = {
        'title': 'Intersection of a set of tracks',
        'description': 'Returns a new track with only regions covered in every input track.',
        'path': ['Features', 'Union'],
        'output': CombineForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        func = any
        output = self.temporary_path(fname='combined.')
        output = _combine(func,output,**kw)
        self.new_file(output, 'combined')
        return self.display_time()
