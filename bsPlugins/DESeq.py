from bsPlugins import *
from bbcflib.track import track
from bbcflib import genrep
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
import numpy
import os,shutil
import itertools

ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100

__requires__ = ["ryp2", "numpy"]


class DESeqForm(BaseForm):
    child = twd.HidingTableLayout()
    input_type = twd.HidingRadioButtonList(label='Input type: ',
                                           options=['Table', 'Signals'],
                                           mapping={'Table': ['table'],
                                                    'Signals': ['Group1','Group2','feature_type','assembly'],},
                                           value='Table',
                                           help_text='Select input type (Formatted table, or signal tracks)')
    table = twb.BsFileField(label='Table: ',
        help_text='Select scores table',
        validator=twb.BsFileFieldValidator(required=True))
    class Group1(twb.BsMultiple):
        label = 'Signals group 1: '
        signals1 = twb.BsFileField(label=' ',
            help_text='Select signal files (position and score, e.g. bedgraph)',
            validator=twb.BsFileFieldValidator(required=True))
    class Group2(twb.BsMultiple):
        label = 'Signals group 2: '
        signals2 = twb.BsFileField(label=' ',
            help_text='Select signal files (position and score, e.g. bedgraph)',
            validator=twb.BsFileFieldValidator(required=True))
    feature_type = twd.HidingSingleSelectField(label='Feature type: ',
        options=ftypes, prompt_text=None,
        mapping={ftypes[-1][0]: ['features'],
                 1: ['upstream', 'downstream']},
        help_text='Choose a feature set or upload your own',
        validator=twc.Validator(required=True))
    features = twb.BsFileField(label='Custom feature set: ',
        help_text='Select a feature file (e.g. bed)',
        validator=twb.BsFileFieldValidator(required=True))
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
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id': 'input_type', 'type': 'radio'},
        {'id': 'signals1', 'type': 'track', 'required': True, 'multiple': 'Group1'},
        {'id': 'signals2', 'type': 'track', 'required': True, 'multiple': 'Group2'},
        {'id': 'table', 'type': 'txt', 'required': True},
        {'id': 'feature_type', 'type': 'int'},
        {'id': 'upstream', 'type': 'int'},
        {'id': 'downstream', 'type': 'int'},
        {'id': 'assembly', 'type': 'assembly'},
        {'id': 'features', 'type': 'track'},
]
out_parameters = [{'id': 'differential_expression', 'type': 'file'}]


class DESeqPlugin(BasePlugin):
    """Gets the score associated to each genomic feature in each sample and runs DESeq
for differential analysis within them. It returns a tab-delimited file with the following fields:

Name, MeanA, MeanB, fold change, adjusted p-value.

The input can be of two different types:

* Two sets of 'signal' files - i.e. bedGraph-type text files - one for each of the two groups to compare -,
  and a list of genomic features - either from a pre-defined list such as Ensembl genes,
  or a custom bed-like file. For every feature, a score is given for each of the signal samples,
  and DESeq is run on the resulting table. The name of each sample is the one given in the track
  definition line ("track name=... description=... etc."), if specified, otherwise the name of
  the file (without extension).
* A tab-delimited table with feature names in the first column, then one column of respective
  scores per sample. The first line is a header of the type "id  sample1  sample2 ...".
  If sample names are in the format 'group_name.run_id', all samples with the same group_name
  will be considered as replicates of the same group/condition. Else they are considered as belonging
  to different groups. If there are more than 2 groups, all different pairs of comparisons
  will be performed and output in separate files.

    """
    info = {
        'title': 'Differential expression analysis',
        'description': __doc__,
        'path': ['Analysis', 'DE analysis'],
        'output': DESeqForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def clean_deseq_output(self,filename,contrast):
        """Delete all lines of *filename* with NA's everywhere, add 0.5 to zero scores
        before recalculating the fold change, remove row numbers, and keep only the following
        fields: Name, MeanA, MeanB, fold change, adjusted p-value. Return the new file name."""
        filename_clean = self.temporary_path()
        with open(filename,"rb") as f:
            with open(filename_clean,"wb") as g:
                f.readline() # header
                A = 'Mean.'+contrast[0].strip()
                B = 'Mean.'+contrast[1].strip()
                g.write('-'.join(contrast)+'\n')
                g.write('\t'.join(['Name',A,B,'foldChange','padj'])+'\n')
                for line in f:
                    line = line.strip().split("\t")
                    if not (line[2]=="0" and line[3]=="0") and not (line[2]=='NA' or line[3]=='NA'):
                        meanA = float(line[2]) or 0.5
                        meanB = float(line[3]) or 0.5
                        fold = meanB/meanA
                        line = '\t'.join([line[0],str(meanA),str(meanB),str(fold),line[7]])+'\n'
                        g.write(line)
        return filename_clean

    def __call__(self, **kw):
        if kw.get('input_type') == 'Table':
            filename = kw.get('table')
            assert os.path.exists(str(filename)), "File not found: '%s'" % filename
            colnames = numpy.asarray(open(filename).readline().split()[1:])
            robjects.r.assign('col_names', numpy2ri.numpy2ri(colnames))
            robjects.r("""
Mdata = read.table('%s',sep='\t',header=T,row.names=1)
conds = unlist(strsplit(col_names,".",fixed=T))
""" % filename)
        else:
            from QuantifyTable import QuantifyTablePlugin
            assembly = genrep.Assembly(kw.get('assembly'))
            chrmeta = assembly.chrmeta or "guess"
            kw['score_op'] = 'sum'
            signals1 = kw['Group1']['signals1']
            signals2 = kw['Group2']['signals2']
            if not isinstance(signals1,(list,tuple)): signals1 = [signals1]
            if not isinstance(signals2,(list,tuple)): signals2 = [signals2]
            signals = signals1+signals2
            kw['SigMulti'] = {'signals': signals} # to pass it to QuantifyTable plugin
            table = QuantifyTablePlugin().quantify(**kw)
            stracks = []
            norm_factors = []
            for sig in signals:
                assert os.path.exists(str(sig)), "Signal file not found: '%s'." % sig
                _t = track(sig, chrmeta=chrmeta)
                if 'normalization' in _t.info:
                    _nf = float(_t.info['normalization'])
                elif 'nreads' in _t.info:
                    _nf = float(_t.info['nreads']) * 1e-7 / float(_t.info.get('read_extension', 1))
                else:
                    _nf = 1
                stracks.append(_t)
                norm_factors.append(_nf)
            t = track(table,chrmeta=chrmeta,fields=['chr','start','end','name']+ \
                                                   ['score%d'%x for x in range(len(signals))])
            _f = [f for f in t.fields if f.startswith('score')]
            de_list = list(t.read(fields=['name']+_f))
            t.close(); os.remove(table)
            # Turn all scores into integers
            de_matrix = numpy.asarray([[int(float(s) * norm_factors[k] + .5) for k,s in enumerate(x[1:])]
                                       for x in de_list], dtype=numpy.float)
            rownames = numpy.asarray([x[0] for x in de_list])
            colnames = numpy.asarray([s.name for s in stracks])
             # if all prefixes are identical within a group, keep this prefix as group identifier.
            if len(list(set( [x.split('.')[0] for x in colnames[:len(signals1)]] ))) == 1 \
                    and len(list(set( [x.split('.')[0] for x in colnames[len(signals1):]] ))) == 1:
                group1 = colnames[0].split('.')[0]
                group2 = colnames[-1].split('.')[0]
            else:
                group1 = "Group1"
                group2 = "Group2"
            conds = [group1]*len(signals1) + [group2]*len(signals2)
            robjects.r.assign('Mdata', numpy2ri.numpy2ri(de_matrix))
            robjects.r.assign('row_names', numpy2ri.numpy2ri(rownames))
            robjects.r.assign('col_names', numpy2ri.numpy2ri(colnames))
            robjects.r.assign('conds', numpy2ri.numpy2ri(conds))
            robjects.r("""
Mdata = as.data.frame(Mdata,row.names=row_names)
colnames(Mdata) = unlist(col_names)
conds = unlist(conds)
""")

        robjects.r("""
library(DESeq)
if (all(table(conds)>=3)){        # if >3 replicates in all conditions
    method = 'per-condition'        # for each group estimate the variance from its replicates
    sharingMode = 'gene-est-only'   # use the per-gene variance estimates only
} else if (any(table(conds)>1)){ # if few replicates
    method = 'pooled'               # use all groups with replicates to estimate the variance
    sharingMode = 'maximum'         # use the max of the GLM fit and the estimated variance
} else {                         # if no replicates
    method = 'blind'                # pools all groups together to estimate the variance
    sharingMode='fit-only'          # use only the GLM fit across the pooled variance
}
cds = newCountDataSet(Mdata, conds)
cds = estimateSizeFactors(cds)
test = try({
    cds = estimateDispersions(cds, method=method, fitType='parametric', sharingMode=sharingMode)
})
if(class(test) == "try-error") {
    cds = estimateDispersions(cds, method=method, fitType='local', sharingMode=sharingMode)
}
""")

        groups = list(set(conds))
        couples = itertools.combinations(groups, 2)
        output = self.temporary_path(fname='DE')
        for c in couples:
            out = "%s_%s-%s.txt" %((output,)+tuple(c))
            robjects.r("""
            res = nbinomTest(cds, '%s', '%s')
            res = res[order(res[,8]),]
            write.table(res, '%s', row.names=F, quote=F, sep='\t')
            """ % (c[0], c[1], out))
            if kw.get('complete') is None:
                clean = self.clean_deseq_output(out,c)
                shutil.move(clean,out)
            self.new_file(out, 'differential_expression')
        return self.display_time()

