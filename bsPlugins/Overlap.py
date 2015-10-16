from bsPlugins import *
from bbcflib.gfminer.stream import overlap
from bbcflib.track import track
from bbcflib import genrep

class OverlapForm(BaseForm):
    child = twd.HidingTableLayout()
    features = twb.BsFileField(label='Features file: ',
        help_text='Upload your own file',
        validator=twb.BsFileFieldValidator(required=True))
    filter = twb.BsFileField(label='Filter file: ',
        help_text='Upload your own file',
        validator=twb.BsFileFieldValidator(required=True))
    format = twf.SingleSelectField(label='Output format: ',
        options=["txt","bed","sql","bedGraph","bigWig"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'filter', 'type': 'userfile', 'required': True, 'label': 'Filter file: ', 'help_text': 'Upload your own file'},
                 {'id': 'features', 'type': 'track', 'required': True, 'label': 'Features file: ', 'help_text': 'Upload your own file'},
                 {'id': 'output', 'type': 'listing', 'label': 'Output format: ', 'help_text': 'Format of the output file','options': ["txt","bed","sql","bedGraph","bigWig"]},
                 {'id': 'assembly', 'type': 'assembly', 'label': 'Assembly: ', 'help_text': 'Reference genome', 'options': genrep.GenRep().assemblies_available()}]
out_parameters = [{'id': 'filtered', 'type': 'track'}]


class OverlapPlugin(BasePlugin):
    """Returns only the regions of the first input file that overlap
(or contain) some feature from the second ('filter')."""
    info = {
        'title': 'Overlap',
        'description': __doc__,
        'path': ['Intervals', 'Overlap'],
#        'output': OverlapForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        # Set assembly
        assembly_id = kw.get('assembly')
        chrmeta = "guess"
        if assembly_id:
            assembly = genrep.Assembly(assembly_id)
            chrmeta = assembly.chrmeta
        # Set features track
        features = track(kw['features'], chrmeta=chrmeta or None )
        chrmeta = features.chrmeta
        # Set filter track
        filter = track(kw.get('filter'), chrmeta=chrmeta or None)
        # Main
        format = kw.get('format',features.format)
        output = self.temporary_path(fname=features.name+'_filtered.'+format)
        tout = track(output, format, fields=filter.fields,
                     chrmeta=chrmeta, info={'datatype':'qualitative'})
        for chrom in chrmeta:
            tout.write(overlap(features.read(chrom),filter.read(chrom)), chrom=chrom,clip=True)
        tout.close()
        self.new_file(output, 'filtered')
        return self.display_time()
