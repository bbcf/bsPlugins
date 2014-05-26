from bsPlugins import *
from bbcflib.gfminer.stream import combine
from bbcflib.track import track, FeatureStream
from bbcflib import genrep


class CombineForm(BaseForm):
    child = twd.HidingTableLayout()
    class TrackMulti(twb.BsMultiple):
        label = 'Tracks: '
        tracks = twb.BsFileField(label=' ',
                                  help_text='Select files to combine',
                                  validator=twb.BsFileFieldValidator(required=True))
    format = twf.SingleSelectField(label='Output format: ',
                                   options=["sql","bed","sga"],
                                   prompt_text=None,
                                   validator=twc.Validator(required=True),
                                   help_text='Format of the output file')
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=genrep.GenRep().assemblies_available(),
                                     help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Quantify")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'tracks', 'type': 'track', 'multiple': 'TrackMulti', 'required': True},
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
    format = kw.get('format') or 'sql'
    output += format
    tracks = kw['TrackMulti']['tracks']
    if not isinstance(tracks, list):
        tracks = [tracks]
    tracks = [track(sig, chrmeta=chrmeta) for sig in tracks]
    chrmeta = tracks[0].chrmeta
    tout = track(output, chrmeta=chrmeta, info={'datatype': 'qualitative'})
    for chrom in chrmeta:
        trackList = [sig.read(chrom) for sig in tracks]
        res = combine(trackList, fn=func)
        tout.fields = res.fields
        tout.write(res, chrom=chrom, clip=True)
    tout.close()
    return output


class IntersectPlugin(BasePlugin):
    info = {
        'title': 'Intersection of a set of tracks',
        'description': 'Returns a new track with only regions covered in every input track.',
        'path': ['Intervals', 'Intersect'],
        'output': CombineForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def _func(self,X):
        return all(X)
    def __call__(self, **kw):
        output = self.temporary_path(fname='combined.')
        output = _combine(self._func,output,**kw)
        self.new_file(output, 'combined')
        return self.display_time()

class UnionPlugin(BasePlugin):
    info = {
        'title': 'Union of a set of tracks',
        'description': 'Returns a new track with regions covered in at least one of the input tracks.',
        'path': ['Intervals', 'Union'],
        'output': CombineForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def _func(self,X):
        return any(X)
    def __call__(self, **kw):
        output = self.temporary_path(fname='combined.')
        output = _combine(self._func,output,**kw)
        self.new_file(output, 'combined')
        return self.display_time()

class SubtractPlugin(BasePlugin):
    info = {
        'title': 'Subtract',
        'description': 'Returns a new track with regions present in the first input track, but not in the others.',
        'path': ['Intervals', 'Subtract'],
        'output': CombineForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def _func(self,X):
        return X[0] and not any(X[1:])
    def __call__(self, **kw):
        output = self.temporary_path(fname='combined.')
        output = _combine(self._func,output,**kw)
        self.new_file(output, 'combined')
        return self.display_time()

class ComplementPlugin(BasePlugin):
    info = {
        'title': 'Complement',
        'description': 'Returns a new track with all regions not covered by a set of input tracks.',
        'path': ['Intervals', 'Complement'],
        'output': CombineForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def _func(self,X):
        """Same as for Subtract, since Complement is subtraction from whole chromosome of a set of features.
           Same as the others with func = 'not any(X)' would forget the extremities of the chromosome."""
        return X[0] and not any(X[1:])
    def __call__(self, **kw):
        # Create a track with the whole chromosome
        chrmeta = _get_chrmeta(**kw)
        sig0 = track(kw['TrackMulti']['tracks'][0])
        fields = sig0.fields
        format = sig0.format
        is_chr = 'chr' in fields
        _f0 = ('chr','start','end') if is_chr else ('start','end')
        _f1 = [f for f in fields if f not in _f0]
        whole_chr = []
        if is_chr:
            for chr in chrmeta:
                whole_chr.append( (chr,0,chrmeta[chr]['length'])+('0',)*len(_f1) )
        else:
            fields = [f for f in fields if f not in ['start','end']]
            fields = ['start','end']+fields
            for chr in chrmeta:
                whole_chr.append( (0,chrmeta[chr]['length'])+('0',)*len(_f1) )
        whole_chr = FeatureStream(whole_chr,fields=fields)
        temp = self.temporary_path()+'.'+format
        with track(temp,fields=fields) as wc:
            wc.write(whole_chr)

        kw['TrackMulti']['tracks'] = [temp] + kw['TrackMulti']['tracks']
        output = self.temporary_path(fname='combined.')
        output = _combine(self._func,output,**kw)
        self.new_file(output, 'combined')
        return self.display_time()

