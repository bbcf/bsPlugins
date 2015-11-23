from bsPlugins import *
import rpy2.robjects as robjects
import os, tarfile

mart_map = [("GRCh37.p6",'hg19'), ("GRCh38",'hg38'), 
            ("NCBIM37",'mm9'), ("GRCm38.p2",'mm10'), 
            ("EF3","sacCer2"), ("BDGP5.25",'dm3'),("Zv9",'zv9')]

default_path = "/mnt/common/epfl/share"

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'gene_list', 'type': 'userfile', 'required': True, 'label': 'Genes: ', 'help_text': 'Provide a list of Ensembl IDs'},
                 {'id': 'assembly', 'type': 'assembly', 'label': 'Assembly: ', 'help_text': 'Reference genome', 'options': mart_map, 'prompt_text': None},
                 {'id': 'num_terms', 'type': 'int', 'label': 'Number of significant terms: ', 'help_text': 'Maximum number of significant terms to return', 'value': 10},
                 {'id': 'pval', 'type': 'float', 'label': 'P-value threshold: ', ' help_text': 'p-value threshold for significance', 'value': 0.05},
                 {'id': 'txid', 'type':'boolean', 'label': 'Transcript identifiers: ', 'help_text': 'Check if this a list of transcript ids (default is gene ids)', 'value': False}]
out_parameters = [{'id': 'TopGO_table_tar', 'type': 'file'},
                  {'id': 'TopGO_plots_tar', 'type': 'file'},
                  {'id': 'TopGO_table', 'type': 'txt'},
                  {'id': 'TopGO_plots', 'type': 'pdf'}]


class TopGoForm(BaseForm):
    gene_list = twb.BsFileField(label='Genes: ',
                              help_text='Provide a list of Ensembl IDs',
                              validator=twb.BsFileFieldValidator(required=True))
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=mart_map,
                                     prompt_text=None,
                                     help_text='Reference genome')
    num_terms = twf.TextField(label='Number of significant terms: ',
                              validator=twc.IntValidator(required=False),
                              value=10,
                              help_text='Maximum number of significant terms to return')
    pval = twf.TextField(label='P-value threshold: ',
                         validator=twb.FloatValidator(min=0,max=1),
                         value=.05,
                         help_text='p-value threshold for significance')
    txid = twf.CheckBox(label='Transcript identifiers: ',
                        value=False,
                        help_text='Check if this a list of transcript ids (default is gene ids)')
    submit = twf.SubmitButton(id="submit", value="TopGo analysis")


class TopGoPlugin(BasePlugin):
    """Makes a GO analysis on a list of Ensembl Gene IDs.

Given a file with one Ensembl Gene ID 
on each line (these are ids like: 'ENSG00000111640', 'ENSMUSG00000057666', 'ENSDARG00000043457', 'YJL153C', ...), it returns a summary table (.txt) and GO networks in a multi-page pdf. 
Multiple lists can be processed in a single job by separating them in the file with a comment line (starting with a "#").

The first summarizes the most significant terms concerning
Biological Processes (BP), Cellular Components (CC) and Molecular Function (MF).
Users can limit the maximum number of terms per category to be displayed in the output,
at a given p-value threshold.
    """
    info = {
        'title': 'Gene Ontology analysis (TopGO)',
        'description': __doc__,
        'path': ['Analysis', 'TopGo'],
#        'output': TopGoForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        assembly_id = kw.get('assembly') or None
        for k,v in mart_map:
            if assembly_id == v:
                assembly_id = k
                break
        if assembly_id is None:
            raise ValueError("Please specify an assembly")
        filename = kw.get('gene_list')
        assert os.path.exists(str(filename)), "File not found: '%s'" %filename
        txid = kw.get('txid',False)
        if isinstance(txid, basestring): txid = (txid.lower() in ['1', 'true', 't','on'])
        use_txids = "transcript" if txid else "gene"
        script_path = kw.get("script_path",default_path)
        fname = os.path.splitext(os.path.basename(filename))[0]
        pdf = self.temporary_path(fname='TopGO_plots.pdf')
        table = self.temporary_path(fname='TopGO_tables.txt')
        num_terms = int(kw.get('num_terms') or 10)
        pval = float(kw.get('pval') or .05)
        robjects.r("""
source("%s/TopGo.R")
out = multi_topGo("%s","%s","%s","%s",%i,%f,"%s")
"""%(script_path,filename,assembly_id,pdf,table,num_terms,pval,use_txids))

        pdf_list = [f[0] for f in robjects.r('out')[0]]
        table_list = [f[0] for f in robjects.r('out')[1]]
        if len(pdf_list) > 1:
            tar_pdf_name = self.temporary_path('TopGO_plots_'+fname+'.tgz')
            tar_pdf = tarfile.open(tar_pdf_name, "w:gz")
            [tar_pdf.add(f,arcname=os.path.basename(f)) for f in pdf_list]
            tar_pdf.close()

            tar_table_name = self.temporary_path(fname='TopGO_tables_'+fname+'.tgz')
            tar_table = tarfile.open(tar_table_name, "w:gz")
            [tar_table.add(f,arcname=os.path.basename(f)) for f in table_list]
            tar_table.close()
            self.new_file(tar_pdf_name, 'TopGO_plots_tar')
            self.new_file(tar_table_name, 'TopGO_table_tar')
        elif len(pdf_list) > 0:
            self.new_file(pdf_list[0],'TopGO_plots')
            self.new_file(table_list[0],'TopGO_table')
        else:
            raise ValueError("Gene list is empty")

        return self.display_time()
