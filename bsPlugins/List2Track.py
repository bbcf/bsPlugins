from bsPlugins import *
from bbcflib.btrack import track,FeatureStream
from bbcflib import genrep
from operator import itemgetter
import os


class List2TrackForm(BaseForm):
    child = twd.HidingTableLayout()
    ids_list = twf.FileField(
        label='IDs list: ',
        help_text='Select the file with the list of IDs',
        validator=twf.FileValidator(required=True))
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bed"],
        prompt_text=None,
        help_text='Format of the output file',
        required=True, )
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome',
        validator=twc.Validator(required=True), )
    feature_type = twf.SingleSelectField(label='Feature type: ',
        options=['genes','exons','transcripts'],
        help_text='Choose the kind of genomic features yo want to annotate',
        validator=twc.Validator(required=True), )
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [
        {'id': 'ids_list', 'type': 'txt', 'required': True, 'multiple': True},
        {'id': 'format', 'type': ''}]
out_parameters = [{'id': 'fulltrack', 'type': 'file'}]


class List2TrackPlugin(BasePlugin):
    description = """Transforms a list of Ensembl IDs into a fully annotated track file."""
    info = {
        'title': 'List2Track',
        'description': description,
        'path': ['Files', 'List2Track'],
        'output': List2TrackForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    # extract (chr,start,end,gene_id|gene_name,score,strand) - bed format
    def genes_annot(self,id,x):
        return (x[5],x[1],x[2],id+'|'+x[0],0.0,x[4])
    def exons_annot(self,id,x):
        return (x[6],x[3],x[4],x[1]+'|'+x[2],0.0,x[5])
    def trans_annot(self,id,x):
        return (x[6],x[2],x[3],x[0]+'|'+x[1],0.0,x[4])
    def __call__(self, **kw):
        assembly = genrep.Assembly(kw.get('assembly'))
        format = kw['format']
        ids_list = kw.get('ids_list')
        assert os.path.exists(str(ids_list)), "File not found: '%s'" % ids_list
        output = self.temporary_path(fname='output.'+format)
        out = track(output,chrmeta=assembly)
        if kw['feature_type'] == 'genes':
            map = assembly.get_gene_mapping()
            get_info = self.genes_annot
        elif kw['feature_type'] == 'exons':
            map = assembly.get_exon_mapping()
            get_info = self.exons_annot
        elif kw['feature_type'] == 'transcripts':
            map = assembly.get_transcript_mapping()
            get_info = self.trans_annot
        def _annotate(ids_list):
            with open(ids_list) as ids_file:
                for id in ids_file:
                    id = id.strip()
                    if map.get(id):
                        yield get_info(id,map.get(id))
                    else:
                        yield ('NA','0','0',id,0.0,'0')
        fields = ['chr','start','end','name','score','strand']
        fulltrack = FeatureStream(_annotate(ids_list),fields=fields)
        out.write(fulltrack)
        self.new_file(output, 'fulltrack')
        return self.display_time()
