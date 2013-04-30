from bsPlugins import *
from bbcflib.btrack import track
from bbcflib import genrep
from bbcflib.bFlatMajor.figure import venn
import rpy2.robjects as robjects
import os

default_path = "/mnt/common/epfl/share"

class VennDiagramWithFilterForm(BaseForm):
    table = twb.BsFileField(label='table: ',
        help_text='Select table',
        validator=twb.BsFileFieldValidator(required=True))
    id_columns = twf.TextField(label='columns id: ',
        validator=twc.Validator(required=True),
        value='',
        help_text='comma separated list of columns id for which Venn diagram will be generated (e.g. 3,5)')
    filters = twf.TextField(label='filters: ',
        validator=twc.Validator(required=True),
        value='',
        help_text='comma separated list of simple filters which will be applied to each corresponding column id before doing the Venn diagram (e.g. >2,<0.05)')
    format = twf.SingleSelectField(label='Output format: ',
        options=["png","pdf"],
        validator=twc.Validator(required=False),
        help_text='Output figure format')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'txt', 'required': True},
                {'id': 'id_columns', 'type': 'txt', 'required': True},
                {'id': 'filters', 'type': 'txt', 'required': True},
                {'id': 'format', 'type': 'format'}
]

out_parameters = [{'id':'venn_diagram', 'type':'file'}]


class VennDiagramWithFilterPlugin(BasePlugin):
    description = """Creates a Venn diagram of a table preliminarly filtered on given columns.
    """
    info = {
        'title': 'Venn Diagram With Filter',
        'description': description,
        'path': ['Graphics', 'Venn Diagram With Filter'],
        'output': VennDiagramWithFilterForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        infile = kw.get('table')
        assert os.path.exists(infile),"File not found: %s ." % infile
        s_cols = kw.get('id_columns')
        s_filters = kw.get('filters')
        names = [chr(i+65) for i in range(len(s_cols.split(',')))] # 'A','B','C',...

        with open(kw['table'],"rb") as f:
            h=f.readline().split('\t')

        colnames=[]
        for i in s_cols.split(','):
            indice = int(i)-1
            if indice <= len(h):
                colnames.append(h[indice])

        script_path = kw.get("script_path",default_path)
        out=robjects.r("""
            source("%s/filterVenn.R")
            filterVenn("%s","%s","%s")
        """ %(script_path,infile,s_cols,s_filters))

        D={}
        for x in out[0].split(','): D[x.split(':')[0]]=int(x.split(':')[1])
        path=script_path
        output = self.temporary_path(fname=path+'/venn_diagram.'+kw['format'])
        print(D)
        venn(D,output=output,legend=colnames,format=kw['format'])
        self.new_file(output, 'venn_diagram')

        return self.display_time()


