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

in_parameters = [{'id': 'track', 'type': 'track', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'promoter', 'type': 'int', 'required': True},
                 {'id': 'intergenic', 'type': 'int', 'required': True},
                 {'id': 'UTR', 'type': 'int', 'required': True}]
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
        'output': AnnotateForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        assembly_id = kw.get('assembly') or None
        assembly = genrep.Assembly(assembly_id)
        tinput = track(kw.get('track'), chrmeta=assembly.chrmeta)
        if kw.get("promoter") is None: thPromot = prom_def
        else:                          thPromot = int(kw["promoter"])
        if kw.get("intergenic") is None: thInter = inter_def
        else:                            thInter = int(kw["intergenic"])
        if kw.get("UTR") is None: thUTR = utr_def
        else:                     thUTR = int(kw["UTR"])
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

