from bsPlugins import *
from bbcflib.btrack import track
from bbcflib import genrep
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import numpy
import os
import itertools


ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100

__requires__ = ["ryp2", "numpy"]


class DESeqForm(BaseForm):
    child = twd.HidingTableLayout()

    input_type = twd.HidingRadioButtonList(label='Input type',
        options=('Table', 'Signals'),
        mapping={'Table':  ['table'],
                 'Signals': ['SigMulti','feature_type','assembly'],},
        help_text='Select input type (Formatted table, or signal tracks)')
    table = twf.FileField(label='Table: ',
        help_text='Select scores table',
        validator=twf.FileValidator(required=True))
    feature_type = twd.HidingSingleSelectField(label='Feature type: ',
        options=ftypes, prompt_text=None,
        mapping={ftypes[-1][0]: ['features'],
                 1: ['upstream', 'downstream']},
        help_text='Choose a feature set or upload your own',
        validator=twc.Validator(required=True))
    class SigMulti(Multi):
        signals = twf.FileField(label='Signals: ',
            help_text='Select signal file (position and score, e.g. bedgraph)',
            validator=twf.FileValidator(required=True))
    features = twf.FileField(label='Custom feature set: ',
        help_text='Select a feature file (e.g. bed)',
        validator=twf.FileValidator())
    upstream = twf.TextField(label='Promoter upstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_up_def,
        help_text='Size of promoter upstream of TSS')
    downstream = twf.TextField(label='Promoter downstream distance: ',
        validator=twc.IntValidator(required=True),
        value=prom_down_def,
        help_text='Size of promoter downstream of TSS')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        validator=twc.Validator(required=True),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Quantify")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id': 'input_type', 'type': 'radio'},
        {'id': 'signals', 'type': 'track', 'required': True, 'multiple': True},
        {'id': 'table', 'type': 'txt', 'required': True, 'multiple': True},
        {'id': 'feature_type', 'type': 'int'},
        {'id': 'upstream', 'type': 'int'},
        {'id': 'downstream', 'type': 'int'},
        {'id': 'assembly', 'type': 'assembly'},
        {'id': 'features', 'type': 'userfile'},
]
out_parameters = [{'id': 'differential_expression', 'type': 'file'}]


class DESeqPlugin(OperationPlugin):

    description = """Gets the score associated to each feature in each sample and runs DESeq
for differential analysis within them. <br /><br />

The input can be of two different types: <br />
* A set of 'signal' files, i.e. bedGraph-type text files,
  and a list of genomic features - either from a pre-defined list such as Ensembl genes,
  or a custom bed-like file. For every feature, a score is given for each of the signal samples,
  and DESeq is run on the resulting table. The name of each sample is the one given in the track
  definition line ("track name=... description=... etc."), if specified, otherwise the name of
  the file (without extension). <br />
* A tab-delimited table with feature names in the first column, then one column of respective
  scores per sample. The first line is a header of the type "id  sample1  sample2 ...". <br />

If sample names are in the format 'group_name.run_id', all samples with
the same group_name will be considered as replicates of the same group/condition.
Else they are considered as belonging to different groups.
    """
    info = {
        'title': 'Differential expression analysis',
        'description': description,
        'path': ['Signal', 'DESeq'],
        'output': DESeqForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):

        assembly = genrep.Assembly(kw.get('assembly'))
        chrmeta = assembly.chrmeta or "guess"

        if kw.get('input_type') == 'Table':
            filename = kw.get('table')
            assert os.path.exists(str(filename)), "File not found: '%s'" % filename
            colnames = numpy.asarray(open(filename).readline().split()[1:])
            robjects.r.assign('col_names', numpy2ri.numpy2ri(colnames))
            robjects.r("""
            Mdata <- read.table('%s',sep='\t',header=T,row.names=1)
            conds <- unlist(strsplit(col_names,".",fixed=T))
            conds <- colnames(Mdata)
            """ % filename)
        else:
            from QuantifyTable import QuantifyTablePlugin
            kw['score_op'] = 'sum'
            table = QuantifyTablePlugin().quantify(**kw)
            signals = kw.get('signals',[])
            stracks = []
            norm_factors = []
            for sig in signals:
                assert os.path.exists(str(sig)), "Signal file not found: '%s'." % sig
                _t = track(sig, chrmeta=chrmeta)
                if 'normalization' in _t.info:
                    print 'normalized'
                    _nf = float(_t.info['normalization'])
                elif 'nreads' in _t.info:
                    print 'nreads'
                    _nf = float(_t.info['nreads']) * 1e-7 / float(_t.info.get('read_extension', 1))
                else:
                    _nf = 1
                stracks.append(_t)
                norm_factors.append(_nf)
            t = track(table,chrmeta=chrmeta)
            _f = [f for f in t.fields if f.startswith('score')]
            de_list = list(t.read(fields=['name']+_f))
            t.close(); os.remove(table)
            # Turn all scores into integers
            de_matrix = numpy.asarray([[int(float(s) * norm_factors[k] + .5) for k,s in enumerate(x[1:])]
                                       for x in de_list], dtype=numpy.float)
            rownames = numpy.asarray([x[0] for x in de_list])
            colnames = numpy.asarray([s.info.get('name',os.path.splitext(os.path.basename(s.path))[0])
                                      for s in stracks])
            robjects.r.assign('Mdata', numpy2ri.numpy2ri(de_matrix))
            robjects.r.assign('row_names', numpy2ri.numpy2ri(rownames))
            robjects.r.assign('col_names', numpy2ri.numpy2ri(colnames))
            robjects.r("""
            Mdata <- as.data.frame(Mdata,row.names=row_names)
            conds <- unlist(strsplit(col_names,".",fixed=T))
            colnames(Mdata) <- conds
            """)

        robjects.r("""
        ### Still need to check that replicates are not identical - lfproc would fail
        groups <- unique(conds)
        couples <- combn(groups,2)
        if (any(table(conds)>1)){ method = 'normal' # if replicates
        } else { method = 'blind' }
        """)

        robjects.r("""
        library(DESeq)
        cds <- newCountDataSet(Mdata, conds)
        cds <- estimateSizeFactors(cds)
        cds <- estimateVarianceFunctions(cds,method='blind')
        """)

        groups = list(set(colnames))
        couples = itertools.combinations(groups, 2)
        output = self.temporary_path(fname='DE')
        for c in couples:
            out = output + '_' + c[0] + '-' + c[1] + '.txt'
            r_cmd = """
            res <- nbinomTest(cds, '%s', '%s')
            res <- res[order(res[,8]),]
            write.table(res, '%s', row.names=F)
            """ % (c[0], c[1], out)
            robjects.r(r_cmd)
            self.new_file(out, 'differential_expression')
        return self.display_time()
