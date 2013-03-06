from bsPlugins import *
from bbcflib.btrack import track
from bbcflib import genrep
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import os, tarfile

mart_map = [("GRCh37.p5",'hg19'), ("NCBIM37",'mm9'), ("EF3","sacCer2"), 
            ("BDGP5.25",'dm3'),("Zv9",'zv9')]

default_path = "/mnt/common/epfl/share"

class TopGoForm(BaseForm):
    gene_list = twf.FileField(label='Genes: ',
                              help_text='Provide a list of ensmbl IDs',
                              validator=twf.FileValidator(required=True))
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=mart_map,
                                     help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="TopGo analysis")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'gene_list', 'type': 'userfile', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'}]
out_parameters = [{'id': 'TopGO_table_tar', 'type': 'file'},
                  {'id': 'TopGO_plots_tar', 'type': 'file'},
                  {'id': 'TopGO_table', 'type': 'txt'},
                  {'id': 'TopGO_plots', 'type': 'pdf'}]


class TopGoPlugin(OperationPlugin):

    info = {
        'title': 'TopGo',
        'description': 'Makes a GO analysis on a list of Ensembl IDs',
        'path': ['Features', 'TopGo'],
        'output': TopGoForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }


    def __call__(self, **kw):
        assembly_id = kw.get('assembly') or None
        filename = kw.get('gene_list')
        assert os.path.exists(str(filename)), "File not found: '%s'" %filename
        script_path = kw.get("script_path",default_path)
        pdf = self.temporary_path(fname='TopGO_plots.pdf')
        table = self.temporary_path(fname='TopGO_tables.txt')
       
        robjects.r("""
source("%s/TopGo.R")
out = multi_topGo("%s","%s","%s","%s")
"""%(script_path,filename,assembly_id,pdf,table))

        pdf_list = list(robjects.r('out')[0])
        table_list = list(robjects.r('out')[1])
        if len(pdf_list) > 1:
            tar_pdf_name = self.temporary_path('TopGO_plots.tgz')
            tar_pdf = tarfile.open(tar_pdf_name, "w:gz")
            [tar_pdf.add(f) for f in pdf_list]
            tar_pdf.close()
            
            tar_table_name = self.temporary_path(fname='TopGO_tables.tgz')
            tar_table = tarfile.open(tar_table_name, "w:gz")
            [tar_table.add(f) for f in table_list]
            tar_table.close()
            self.new_file(tar_pdf_name, 'TopGO_plots_tar')
            self.new_file(tar_table_name, 'TopGO_table_tar')
        else:
            self.new_file(pdf,'TopGO_plots')
            self.new_file(table,'TopGO_table')
                
        return self.display_time()
