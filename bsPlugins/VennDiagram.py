from bsPlugins import *
from bbcflib.btrack import track
from bbcflib.bFlatMajor.common import *
from bbcflib.bFlatMajor.stream import concatenate
from bbcflib.bFlatMajor.figure import venn
from bbcflib import genrep
from itertools import combinations
import os

# nosetests --logging-filter=-tw2 test_VennDiagram.py

class VennDiagramForm(BaseForm):
    child = twd.HidingTableLayout()
    class SigMulti(Multi):
        label = "Files: "
        files = twb.BsFileField(label=' ',
            help_text='Select your track files',
            validator=twb.BsFileFieldValidator(required=True))
    format = twf.SingleSelectField(label='Format: ',
        options=['pdf','png','jpeg'],
        help_text='Output file format')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        validator=twc.Validator(required=True),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id':'files', 'type':'track', 'required':True, 'multiple':True},
        {'id':'format', 'type':'list'},
        {'id':'assembly', 'type':'assembly'},
]
out_parameters = [{'id':'venn_diagram', 'type':'file'}]


class VennDiagramPlugin(BasePlugin):
    description = """

    """
    info = {
        'title': 'Venn Diagram',
        'description': description,
        'path': ['Graphics', 'Venn Diagram'],
        'output': VennDiagramForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        assembly = genrep.Assembly(kw['assembly'])
        filenames = kw.get('files',[])
        if not isinstance(filenames,(list,tuple)): filenames = [filenames]
        for f in filenames: assert os.path.exists(f), "Fie not found: %s ." % filename
        tracks = [track(f) for f in filenames]
        track_names = [chr(i+65) for i in range(len(tracks))] # file name?, or 'A','B','C',...
        combn = [combinations(track_names,k) for k in range(1,len(tracks)+1)]
        combn = ['|'.join(sorted(y)) for x in combn for y in x]
        sets = dict(zip(combn,[0]*len(combn)))
        def _f(i): # hack
            return lambda x:track_names[i]
        for chrom in assembly.chrmeta:
            streams = [t.read(chrom) for t in tracks]
            streams = [duplicate(s,'chr','track_name') for s in streams]
            streams = [apply(s,'track_name',_f(i)) for i,s in enumerate(streams)]
            s = concatenate(streams, aggregate={'track_name':lambda x:'|'.join(x)})
            s = cobble(s)
            name_idx = s.fields.index('track_name')
            for x in s:
                # Add 1 to each sub-category piece x belongs to
                sub = sorted(x[name_idx].split('|'))
                cb = [combinations(sub,k) for k in range(1,len(sub)+1)]
                cb = ['|'.join(sorted(y)) for x in cb for y in x]
                for c in cb: sets[c] += 1
        venn_options = {} # tune it here
        output = self.temporary_path(fname='venn_diagram.'+kw['format'])
        venn(sets,options=venn_options,output=output,format=kw['format'])
        self.new_file(output, 'venn_diagram')
        return self.display_time()
