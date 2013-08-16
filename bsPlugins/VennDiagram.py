from bsPlugins import *
from bbcflib.track import track, FeatureStream
from bbcflib.gfminer.common import cobble
from bbcflib.gfminer.stream import concatenate
from bbcflib.gfminer.figure import venn
from bbcflib import genrep
from itertools import combinations
import os, re, sys

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'input_type', 'type': 'radio'},
                 {'id': 'files', 'type':'track', 'multiple':'TrMulti'},
                 {'id': 'type', 'type':' list'},
                 {'id': 'format', 'type': 'list'},
                 {'id': 'table', 'type': 'track'},
                 {'id': 'id_columns', 'type': 'text'},
                 {'id': 'filters', 'type': 'text'}]

out_parameters = [{'id':'venn_diagram', 'type':'file'},
                  {'id':'venn_summary', 'type':'file'}]



class VennDiagramForm(BaseForm):
    child = twd.HidingTableLayout()

    input_type = twd.HidingRadioButtonList(label='Input from: ',
                                           options=['Table', 'Tracks'],
                                           mapping={'Table': ['table','id_columns','filters'],
                                                    'Tracks': ['files','type']},
                                           value='Table',
                                           help_text='Select input type (Formatted table, or genomic tracks)')


    class TrMulti(Multi):
        label = "Files: "
        files = twb.BsFileField(label=' ',
                                help_text='Select your track files',
                                validator=twb.BsFileFieldValidator(required=True))
    type = twf.SingleSelectField(label='Type: ',
                                 prompt_text=None,
                                 options=['intervals','score'],
                                 help_text='Output figure format')
    format = twf.SingleSelectField(label='Format: ',
                                   prompt_text=None,
                                   options=['png','pdf'],
                                   help_text='Output figure format')

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
                            help_text='comma separated list of simple filters which will be applied to each corresponding column id before doing the Venn diagram (e.g. >2,<0.05,>=2 OR <=-2,>=-2 AND <2,==2,!=2) - one filter per column should be given - leave an empty string if no filter should be applied to a given column (e.g., >2,,<0.05)')

    submit = twf.SubmitButton(id="submit", value="Run")


class VennDiagramPlugin(BasePlugin):
    """
Creates a Venn diagram of the proportions of
total coverage/total score attributed to each of the given tracks.

If the parameter 'type' has the value 'intervals', the diagram will show the percent
of the genome covered by each possible combination of the input tracks. For instance,
If tracks A and B are given, it will show the portion covered by A only, B only, or
A and B (where they intersect).

If it has the value 'score', the diagram will show the percent of the total score
due to each combination of the input tracks, as above.

The output includes the figure of the Venn diagram and a text summary of the different statistics.
If more than 4 samples are given, no graph is produced, but the text summary still contains
all the information.
"""

    info = {
        'title': 'Venn Diagram',
        'description': __doc__,
        'path': ['Graphics', 'Venn Diagram'],
        'output': VennDiagramForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):

        def _parse_logic(string):
            s = re.sub(r'[^\w\d!=><\. ]', '', string)
            s = re.sub(r' OR ', ')or(%f ', s)
            s = re.sub(r' AND ', ')and(%f ', s)
            return "(%f "+s+")"

        def _run_test(row, indx, cond):
            num = float(row[col_ind[indx]])
            num = max(-sys.maxint,min(sys.maxint,num))
            num = (num,)*c.count("%f")
            return eval(cond % (num))

        def _add_label(s,x):
            _f = s.fields+['track_name']
            return FeatureStream((y+(x,) for y in s), fields=_f)

        venn_options = {} # tune it here
        tracks = []
        intype = kw.get("input_type") or "Table"
        if intype == "Table":
            s_cols = kw.get('id_columns','')
            s_filters = kw.get('filters','')
            infile = track(kw.get('table',''),format='txt',header=True)
            col_ind = [int(i)-1 for i in s_cols.split(",")]
            legend = [infile.fields[i] if i<len(infile.fields) else str(i) for i in col_ind]
            conds = [_parse_logic(x) for x in s_filters.split(",")]
            tlabels = [chr(k+65) for k in range(len(col_ind))]
            conds += ["1"]*(len(col_ind)-len(conds))
            combn = [tuple(sorted(x)) for k in range(len(tlabels)) 
                     for x in combinations(tlabels,k+1)]
            c1 = dict(("|".join(c),0) for c in combn)
            c2 = dict(("|".join(c),0) for c in combn)
            indx = dict((c,[tlabels.index(x) for x in c]) for c in combn)
            for row in infile:
                tests = [_run_test(row,i,c) for i,c in enumerate(conds)]
                for c in combn:
                    c1["|".join([tlabels[n] for n,t in enumerate(tests) if t])] += 1
                    c2["|".join(c)] += all([tests[i] for i in indx[c]])
            tracks = [infile]
            combn = ['|'.join(y) for x in combn for y in x]
        elif intype == "Tracks":
            filenames = kw['TrMulti']['files']
            if not isinstance(filenames,(list,tuple)): filenames = [filenames]
            for f in filenames: assert os.path.exists(f), "File not found: %s ." % f
            tracks = [track(f,chrmeta='guess') for f in filenames]
            tlabels = [chr(k+65) for k in range(len(tracks))]
            combn = [combinations(tlabels,k+1) for k in range(len(tlabels))]
            combn = ['|'.join(sorted(y)) for x in combn for y in x]
            c1 = dict(zip(combn,[0]*len(combn)))
            c2 = dict(zip(combn,[0]*len(combn)))
            total_cov = 0.0
            _scored = (kw.get('type') == 'score')
            chromset = set([c for t in tracks for c in t.chrmeta])
            for chrom in chromset:
                streams = [_add_label(t.read(chrom),tlabels[n]) for n,t in enumerate(tracks)]
                s = cobble(concatenate(streams),scored=_scored)
                name_idx = s.fields.index('track_name')
                start_idx = s.fields.index('start')
                end_idx = s.fields.index('end')
                if _scored: score_idx = s.fields.index('score')
                for x in s:
                    length = x[end_idx]-x[start_idx]
                    total_cov += length
                    sub = sorted(list(set(x[name_idx].split('|')))) # avoid 'A|A'
                    cb = [combinations(sub,k) for k in range(1,len(sub)+1)]
                    cb = ['|'.join(sorted(y)) for c in cb for y in c]
                    if _scored:
                        c1['|'.join(sub)] += x[score_idx]
                        for c in cb: c2[c] += x[score_idx]
                    else:
                        c1['|'.join(sub)] += length
                        for c in cb: c2[c] += length
            if total_cov < 1:
                output = self.temporary_path(fname='venn_summary.txt')
                with open(output,'wb') as summary:
                    summary.write("Empty content (no coverage) on %s." %(",".join(chromset)))
                self.new_file(output, 'venn_summary')
                return
            legend = [t.name for t in tracks]
            if _scored:
                for c in combn:
                    c2[c] = round(c2[c])
            else:
                for c in combn:
                    c2[c] = round((100*c2[c])/total_cov)
                    c1[c] = (100*c1[c])/total_cov
        else:
            raise ValueError("Input type '%s' not supported." %intype)


        format = kw.get('format') or 'pdf'
        output = self.temporary_path(fname='venn_diagram.'+format)
        if len(tracks) <= 4:
            venn(c2,legend=legend,options=venn_options,output=output,format=format)
        self.new_file(output, 'venn_diagram')

        # Text summary
        output = self.temporary_path(fname='venn_summary.txt')
        with open(output,'w') as summary:
            summary.write("%s\t%s\t%s\n" % ("Group","Coverage", "Cumulative coverage"))
            record = "%s\t%.2f\t%d\n"
            for c in sorted(combn, key=lambda x:(len(x),x)):
                summary.write(record%(c,c1[c],c2[c]))
        self.new_file(output, 'venn_summary')
        return self.display_time()

