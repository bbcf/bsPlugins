from bsPlugins import *
from bbcflib.btrack import track
from bbcflib.bFlatMajor.common import *
from bbcflib.bFlatMajor.stream import concatenate
from bbcflib import genrep
import os


class VennDiagramForm(BaseForm):
    child = twd.HidingTableLayout()
    class SigMulti(Multi):
        label = "Files: "
        files = twb.BsFileField(label=' ',
            help_text='Select your track files',
            validator=twb.BsFileFieldValidator(required=True))
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
        tracks = [track(f) for f in filenames]
        track_names = [chr(i) for i in range(65,65+len(tracks))] # 'A','B','C',...
        sets = {}
        def _f(i): # hack
            return lambda x:track_names[i]
        for chrom in assembly.chrmeta:
            streams = [t.read(chrom) for t in tracks]
            s0 = [duplicate(s,'chr','track_name') for s in streams]
            s1 = [apply(s,'track_name',_f(i)) for i,s in enumerate(s0)]
            s2 = concatenate(s1, aggregate={'track_name':lambda x:'|'.join(x)})
            s3 = cobble(s2)
            name_idx = s3.fields.index('track_name')
            for x in s3:
                names = '|'.join(sorted(x[name_idx].split('|')))
                sets.setdefault(names,[]).append(x)
            if chrom=='chr1':print sets


        output = self.temporary_path(fname='venn_diagram.png')
        self.new_file(output, 'venn_diagram')
        return self.display_time()
