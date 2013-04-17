from bbcflib.btrack import track
from bbcflib import genrep
import os

class Table2TracksForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
        help_text='Select table',
        validator=twb.BsFileFieldValidator(required=True))
    id_columns = twf.TextField(label='columns id: ',
        validator=twc.IntValidator(required=True),
        value='',
        help_text='comma separated list of columns id for which signal tracks will be generated (e.g. 3,5)')
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bedgraph","bigwig","wig"],
        validator=twc.Validator(required=False),
        help_text='Output file(s) format (default: bedGraph)')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'txt', 'required': True},
                {'id': 'id_columns', 'type': 'txt', 'required': True},
                {'id': 'assembly', 'type': 'assembly', 'required': True},
                {'id': 'format', 'type': 'format'}
]

out_parameters = [{'id': 'output', 'type': 'file'}]

#kw= {'tableFile':"Dup_vs_Ctrl_resSelectedFrags_fromSmoothed_KCTD13_chr16.txt",
#	'colnames':['baseMeanA_Ctrl','baseMeanB_Dup','log2FoldChange_Dup_vs_Ctrl','pval'],
#	'assembly':"hg19",
#	'format':"bedGraph"
#}

class Table2TracksPlugin(BasePlugin):
    description = """generate signal tracks from a tab-delimited table"""
    info = {
        'title': 'Table2Tracks',
        'description': description,
        'path': ['Files', 'Table2Tracks'],
        'output': Table2TracksForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        assembly = genrep.Assembly(kw.get('assembly'))
	chrmeta = assembly.chrmeta or "guess"

        with open(kw['tableFile'],"rb") as f:
        h=f.readline().strip().replace('#','').split('\t')

        colnames=[]
        for i in kw['id_columns'].split(','):
            indice = int(i)-1
        if indice <= len(h) and indice > 2: #columns 0,1,2 are for chr,start,end
            colnames.append(h[indice])

        t=track(kw['tableFile'],chrmeta=chrmeta, fields=h)
        (filepath, filename) = os.path.split(kw['tableFile'])
        (shortname, extension) = os.path.splitext(filename)

        for _f in colnames:
        if kw['format']=="":
        out_name = shortname+'_'+_f+'.bedGraph'
        else:
        out_name = shortname+'_'+_f+'.'+kw['format']
        output_name = self.temporary_path(out_name)
        print output_name
        print _f
        out_track = track(output_name,chrmeta=chrmeta)
        s = t.read(fields=['chr','start','end',_f])
        s.fields[3] = "score"
        out_track.write(s, mode='write')
	out_track.close()

	return self.display_time()

