from bsPlugins import *
from bbcflib.track import track,FeatureStream
from bbcflib import genrep
import os


class List2TrackForm(BaseForm):
    child = twd.HidingTableLayout()
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome',
        validator=twc.Validator(required=True), )
    feature_type = twf.SingleSelectField(label='Feature type: ',
        options=['genes','exons','transcripts'],
        prompt_text=None,
        help_text='Choose the kind of genomic features yo want to annotate',
        validator=twc.Validator(required=True), )
    ids_list = twf.FileField(
        label='IDs list: ',
        help_text='Select the file with the list of IDs',)
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bed"],
        prompt_text=None,
        help_text='Format of the output file', )
    submit = twf.SubmitButton(id="submit", value="Submit")


out_opts=["sql","bed"],
meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [
        {'id': 'ids_list', 'type': 'txt', 'required': True, 'label': 'IDs list: ', 'help_text': 'Select the file with the list of IDs'},
        {'id': 'output', 'type': 'listing', 'label': 'Output format: ', 'help_text': 'Format of the output file', 'options': ["sql","bed"], 'prompt_text': None}]
out_parameters = [{'id': 'fulltrack', 'type': 'track'}]


class List2TrackPlugin(BasePlugin):
    """Create a fully annotated track file from a features type or a subset of Ensembl IDs.
    
Either upload a raw text file with one Ensembl ID on each line, or choose a feature type to fetch them all."""
    info = {
        'title': 'Genome track from IDs',
        'description': __doc__,
        'path': ['Files', 'List2Track'],
#        'output': List2TrackForm,
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
        ids_list = kw.get('ids_list')
        fields = ['chr','start','end','name','score','strand']
        if ids_list:
            assert os.path.exists(str(ids_list)), "File not found: '%s'" % ids_list
            fulltrack = FeatureStream(_annotate(ids_list),fields=fields)
            fname = os.path.splitext(os.path.basename(ids_list))[0]
        else:
            fulltrack = FeatureStream((get_info(g,map[g]) for g in map),fields=fields)
            fname = kw['feature_type']
        output = self.temporary_path(fname=fname+'.'+format)
        out = track(output,chrmeta=assembly)
        out.write(fulltrack)
        self.new_file(output, 'fulltrack')
        return self.display_time()

# nosetests --logging-filter=-tw2 test_List2Track.py
