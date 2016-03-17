from bsPlugins import *
from bbcflib.gfminer.stream import getNearestFeature
from bbcflib.track import track
from bbcflib import genrep

prom_def = 2000
inter_def = 20000
utr_def = 10


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'track', 'type': 'track', 'required': True, 'label': 'Features', 'help_text': 'Select features file (e.g. bed)'},
                 {'id': 'assembly', 'type': 'assembly', 'label': 'Assembly', 'help_text': 'Reference genome', 'options': genrep.GenRep().assemblies_available(), 'prompt_text': None},
                 {'id': 'promoter', 'type': 'int', 'required': True, 'label': 'Promoter size: ', 'help_text': 'Upstream distance from TSS in bp to be included in the promoter', 'value': prom_def},
                 {'id': 'intergenic', 'type': 'int', 'required': True, 'label': 'Intergenic distance: ', 'help_text': 'Maximum distance to be associated with a gene', 'value': inter_def},
                 {'id': 'UTR', 'type': 'int', 'required': True, 'label': "3' UTR ratio: ", 'help_text': "3' UTR to promoter ratio in %", 'value': utr_def}]
out_parameters = [{'id': 'table', 'type': 'file'}]


class AnnotateForm(BaseForm):
    track = twb.BsFileField(label='Features: ',
                            help_text='Select features file (e.g. bed)',
                            validator=twb.BsFileFieldValidator(required=True))
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=genrep.GenRep().assemblies_available(),
                                     prompt_text=None,
                                     help_text='Reference genome')
    promoter = twf.TextField(label='Promoter size: ',
                             validator=twc.IntValidator(required=True),
                             value=prom_def,
                             help_text='Upstream distance from TSS in bp to be included in the promoter')
    intergenic = twf.TextField(label='Intergenic distance: ',
                               validator=twc.IntValidator(required=True),
                               value=inter_def,
                               help_text='Maximum distance to be associated with a gene')
    UTR = twf.TextField(label="3' UTR ratio: ",
                        validator=twc.IntValidator(required=True),
                        value=utr_def,
                        help_text="3' UTR to promoter ratio in %")
    submit = twf.SubmitButton(id="submit", value="Annotate")


class AnnotatePlugin(BasePlugin):
    """Searches closest gene to each feature and returns associated distance and inclusion informations"""
    info = {
        'title': 'Associate features with genome annotations',
        'description': __doc__,
        'path': ['Analysis', 'Annotate'],
#        'output': AnnotateForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        assembly_id = kw.get('assembly') or None
        assembly = genrep.Assembly(assembly_id)
        tinput = track(kw.get('track'), chrmeta=assembly.chrmeta)
        try:
            thPromot = int(kw.get("promoter"))
        except (ValueError, TypeError):
            thPromot = prom_def
        try:
            thInter = int(kw.get("intergenic"))
        except (ValueError, TypeError):
            thInter = inter_def
        try:
            thUTR = int(kw.get("UTR"))
        except (ValueError, TypeError):
            thUTR = utr_def
        output = self.temporary_path(fname=tinput.name+'_annotated.txt')
        _fields = tinput.fields+['gene', 'location_type', 'distance']
        tout = track(output, format='txt', fields=_fields)
        tout.make_header("#"+"\t".join(tout.fields))
        for chrom in assembly.chrnames:
            tout.write(getNearestFeature(
                    tinput.read(selection=chrom),
                    assembly.gene_track(chrom),
                    thPromot, thInter, thUTR), mode='append')
        tout.close()
        self.new_file(output, 'table')
        return self.display_time()

