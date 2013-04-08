from bsPlugins import *
from bbcflib.bFlatMajor.stream import overlap
from bbcflib.btrack import track
from bbcflib import genrep

prom_up_def = 1000
prom_down_def = 100
ftypes = [(0,'gene bodies'), (1,'gene promoters'), (2,'exons'), (3,'custom upload')]

class OverlapForm(BaseForm):
    child = twd.HidingTableLayout()
    feature_type = twd.HidingSingleSelectField(label='Input type: ',
       options=ftypes, prompt_text=None,
       mapping={ftypes[-1][0]: ['features'],
                            1: ['upstream_1','downstream_1']},
       help_text='Choose a feature set or upload your own',
       validator=twc.Validator(required=True))
    features = twf.FileField(label='Custom file: ',
        help_text='Upload your own file',
        validator=twf.FileValidator(required=True))
    filter_type = twd.HidingSingleSelectField(label='Filter type: ',
       options=ftypes, prompt_text=None,
       mapping={ftypes[-1][0]: ['features'],
                            1: ['upstream_2','downstream_2']},
       help_text='Choose a filter set or upload your own',
       validator=twc.Validator(required=True))
    filter = twf.FileField(label='Custom Filter file: ',
                           help_text='Upload your own file',
                           validator=twf.FileValidator(required=True))
    format = twf.SingleSelectField(label='Output format: ',
        options=["txt","bed","sql","bedGraph","bigWig"],
        validator=twc.Validator(required=True),
        help_text='Format of the output file')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    upstream_1 = twf.TextField(label='Promoter upstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_up_def,
        help_text='Size of promoter upstream of TSS')
    downstream_1 = twf.TextField(label='Promoter downstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_down_def,
        help_text='Size of promoter downstream of TSS')
    upstream_2 = twf.TextField(label='Promoter upstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_up_def,
        help_text='Size of promoter upstream of TSS')
    downstream_2 = twf.TextField(label='Promoter downstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_down_def,
        help_text='Size of promoter downstream of TSS')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'filter', 'type': 'userfile', 'required': True},
                 {'id': 'filter_type', 'type': 'list'},
                 {'id': 'feature_type', 'type': 'list'},
                 {'id': 'features', 'type': 'userfile', 'required': True},
                 {'id': 'format', 'type': 'text'},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'upstream_1', 'type': 'int'},
                 {'id': 'downstream_1', 'type': 'int'},
                 {'id': 'upstream_2', 'type': 'int'},
                 {'id': 'downstream_2', 'type': 'int'}]
out_parameters = [{'id': 'filtered', 'type': 'track'}]


class OverlapPlugin(BasePlugin):
    info = {
        'title': 'Overlap',
        'description': "Returns only the regions of the first input file that overlap \
                        (or contain) some feature from the second ('filter').",
        'path': ['Features', 'Overlap'],
        'output': OverlapForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        feature_type = int(kw.get('feature_type', 0))
        assembly_id = kw.get('assembly')
        filter = kw.get("filter")
        assert os.path.exists(str(filter)), "File file not found: '%s'." % filter
        filter = os.path.abspath(filter)
        # Set assembly
        chrmeta = "guess"
        if assembly_id:
            assembly = genrep.Assembly(assembly_id)
            chrmeta = assembly.chrmeta
            genes = assembly.gene_track
            exons = assembly.exon_track
        elif not(feature_type == 3):
            raise ValueError("Please specify an assembly")
        # Set features track
        if feature_type == 0: features = genes
        elif feature_type == 1:
            prom_pars = {'before_start': int(kw.get('upstream_1') or prom_up_def),
                         'after_start': int(kw.get('downstream_1') or prom_down_def),
                         'on_strand': True}
            features = lambda c: neighborhood(genes(c), **prom_pars)
        elif feature_type == 2: features = exons
        elif feature_type == 3:
            assert os.path.exists(str(kw.get('features'))), "Features file not found: '%s'"%kw.get("features")
            _t = track(kw.get('features'), chrmeta=chrmeta)
            chrmeta = _t.chrmeta
            features = _t.read
        else: raise ValueError("Take feature_type in %s."%ftypes)
        # Set filter track
        if filter_type == 0: filter = genes
        elif filter_type == 1:
            prom_pars = {'before_start': int(kw.get('upstream_2') or prom_up_def),
                         'after_start': int(kw.get('downstream_2') or prom_down_def),
                         'on_strand': True}
            filter = lambda c: neighborhood(genes(c), **prom_pars)
        elif filter_type == 2: filter = exons
        elif filter_type == 3:
            assert os.path.exists(str(kw.get('filter'))), "Filter file not found: '%s'"%kw.get("filter")
            _t = track(kw.get('filter'), chrmeta=chrmeta)
            chrmeta = _t.chrmeta
            features = _t.read
        else: raise ValueError("Take filter_type in %s."%ftypes)
        # Main
        tfilter = track(filter,chrmeta=chrmeta)
        format = kw.get('format',tfilter.format)
        output = self.temporary_path(fname='filtered.'+format)
        tout = track(output, format, fields=tfilter.fields,
                     chrmeta=chrmeta, info={'datatype':'qualitative'})
        for chrom in chrmeta:
            tout.write(overlap(features(chrom),tfilter.read(chrom)), chrom=chrom,clip=True)
        tout.close()
        self.new_file(output, 'filtered')
        return self.display_time()
