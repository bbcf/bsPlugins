from bsPlugins import *
from bbcflib.btrack import track
from bbcflib import genrep
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri

mart_map = [("GRCh37.p5",'hg19'), ("NCBIM37",'mm9'), ("EF3","sacCer2"), 
            ("BDGP5.25",'dm3'),("Zv9",'zv9')]

ensembl_url = "sep2011.archive.ensembl.org"
biomart = "ENSEMBL_MART_ENSEMBL"
attribute_go = "go_id"
attribute_gene = "external_gene_id"
filter_go = "with_go"


class TopGoForm(BaseForm):
    gene_list = twf.FileField(label='Genes: ',
                              help_text='Provide a list of ensmbl IDs',
                              validator=twf.FileValidator(required=True))
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=mart_map().keys(),
                                     help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="TopGo analysis")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'gene_list', 'type': 'userfile', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'}]
out_parameters = [{'id': 'TopGO_table', 'type': 'txt'},
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
       assert os.path.exists(str(filename)), "File not found: '%s'" % filename
       pdf = self.temporary_path(fname='TopGO_plots.pdf')
       table = self.temporary_path(fname='TopGO_tables.txt')
       robjects.r("""
id_set = scan("%s",what=character())
pdf("%s",paper="a4",height=8,width=11)
output = "%s"
library(biomaRt)
library(topGO)
ensembl = useMart("%s",host="%s")
lsds = listDatasets(ensembl)
dataset = lsds$dataset[which(lsds$version == "%s")]
ensembl = useDataset(as.character(dataset),mart=ensembl)
attr1 = "%s"
attr2 = "%s"
filt = "%s"
genome = getBM(attr=c("ensembl_gene_id",attr1,attr2),
  filter=c(filt,"biotype"),
  values=list(TRUE,"protein_coding"),mart=ensembl)
gene2GO = split(genome[,attr1],genome[,"ensembl_gene_id"])
I = which(!duplicated(genome[,"ensembl_gene_id"]))
genome = data.frame(gene_name=genome[I,attr2],row.names=genome[I,"ensembl_gene_id"])
allGenes = row.names(genome)
genome[which(!nchar(genome[,1])),1] = allGenes[which(!nchar(genome[,1]))]
tab = list()
geneList = factor(as.integer(allGenes %in% id_set))
names(geneList) = allGenes
append = FALSE
for (ontol in c("BP","CC","MF")) {
    data = new("topGOdata",description=ontol,ontology=ontol,
               allGenes=geneList,annot=annFUN.gene2GO,gene2GO=gene2GO)
    result = list(classicFisher=runTest(data, statistic="fisher", algo="classic"),
                  elimFisher=runTest(data, statistic="fisher", algo="elim"))
    tab = GenTable(data, classicFisher=result$classicFisher, 
                   elimFisher=result$elimFisher, 
                   orderBy="elimFisher", ranksOf="elimFisher", topNodes=10)
    showSigOfNodes(data,score(result$elimFisher), first=10, useInfo="all")
    tab = cbind(tab,genes=sapply(tab$GO.ID, function(x){
             paste(sort(genome[unlist(genesInTerm(data,x)),1]),collapse=', ')
               }))
    write.table(tab,file=output,quote=F,sep="\t",row.names=F,append=append)
    append = TRUE
}
dev.off()
         """ % (filename,pdf,table,biomart,ensembl_url,assembly_id,
                attribute_go,attribute_gene,filter_go))
       self.new_file(pdf, 'TopGO_plots')
       self.new_file(table, 'TopGO_table')
       return self.display_time()
