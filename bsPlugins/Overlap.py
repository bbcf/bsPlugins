from bsPlugins import *
from bbcflib.bFlatMajor.stream import overlap
from bbcflib.btrack import track
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

in_parameters = [{'id': 'filter', 'type': 'userfile', 'required': True},
                 {'id': 'features', 'type': 'userfile', 'required': True},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'assembly', 'type': 'assembly'}]
out_parameters = [{'id': 'filtered', 'type': 'track'}]


class OverlapPlugin(BasePlugin):
    info = {
        'title': 'Overlap',
        'description': "Returns only the regions of the first input file that overlap \
                        (or contain) some feature from the second ('filter').",
        'path': ['Intervals', 'Overlap'],
        'output': OverlapForm,
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
        assert os.path.exists(str(kw.get('features'))), "Features file not found: '%s'"%kw.get("features")
        features = track(kw.get('features'), chrmeta=chrmeta or None )
        chrmeta = _t.chrmeta
        # Set filter track
        assert os.path.exists(str(kw.get('filter'))), "Filter file not found: '%s'"%kw.get("filter")
        filter = track(kw.get('filter'), chrmeta=chrmeta or None)
        # Main
        format = kw.get('format',filter.format)
        output = self.temporary_path(fname='filtered.'+format)
        tout = track(output, format, fields=filter.fields,
                     chrmeta=chrmeta, info={'datatype':'qualitative'})
        for chrom in chrmeta:
            tout.write(overlap(features.read(chrom),filter.read(chrom)), chrom=chrom,clip=True)
        tout.close()
        self.new_file(output, 'filtered')
        return self.display_time()
