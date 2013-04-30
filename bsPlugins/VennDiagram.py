from bsPlugins import *
from bbcflib.btrack import track
from bbcflib.bFlatMajor.common import *
from bbcflib.bFlatMajor.stream import concatenate
from bbcflib.bFlatMajor.figure import venn
from bbcflib import genrep
from itertools import combinations
import os


class VennDiagramForm(BaseForm):
    child = twd.HidingTableLayout()
    class SigMulti(Multi):
        label = "Files: "
        files = twb.BsFileField(label=' ',
            help_text='Select your track files',
            validator=twb.BsFileFieldValidator(required=True))
    names = twf.TextArea(label='Sample names: ',
        placeholder="(One name per line)",
        help_text="Sample names (facultative)" )
    type = twf.SingleSelectField(label='Type: ',
        prompt_text=None,
        options=['coverage %','tag count'],
        help_text='Output figure format')
    format = twf.SingleSelectField(label='Format: ',
        prompt_text=None,
        options=['png','pdf'],
        help_text='Output figure format')
    assembly = twf.SingleSelectField(label='Assembly: ',
        prompt_text=None,
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id':'files', 'type':'track', 'required':True, 'multiple':'SigMulti'},
        {'id':'names', 'type':'text'},
        {'id':'type', 'type':'list'},
        {'id':'format', 'type':'list'},
        {'id':'assembly', 'type':'assembly'},
]
out_parameters = [{'id':'venn_diagram', 'type':'file'},
                  {'id':'venn_summary', 'type':'file'}]


class VennDiagramPlugin(BasePlugin):
    description = """Creates a Venn diagram of the proportions of
    total coverage/total tag count due to each of the given tracks.
    """
    info = {
        'title': 'Venn Diagram',
        'description': description,
        'path': ['Graphics', 'Venn Diagram'],
        'output': VennDiagramForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        self.debug(**kw)
        assembly = genrep.Assembly(kw['assembly'])
        filenames = kw['SigMulti']['files']
        if not isinstance(filenames,(list,tuple)): filenames = [filenames]
        for f in filenames: assert os.path.exists(f), "File not found: %s ." % f
        tracks = [track(f) for f in filenames]
        group_names = [n.strip() for n in kw['names'].split('\n')]
        track_names = group_names if len(group_names) == len(tracks) \
                      else [chr(i+65) for i in range(len(tracks))] # 'A','B','C',...
        combn = [combinations(track_names,k) for k in range(1,len(tracks)+1)]
        combn = ['|'.join(sorted(y)) for x in combn for y in x]
        cnt = dict(zip(combn,[0]*len(combn)))
        cumcnt = dict(zip(combn,[0]*len(combn)))
        cov = dict(zip(combn,[0]*len(combn)))
        cumcov = dict(zip(combn,[0]*len(combn)))
        def _f(i): # hack
            return lambda x:track_names[i]
        total_cov = 0.0
        for chrom in assembly.chrmeta:
            streams = [t.read(chrom) for t in tracks]
            streams = [duplicate(s,'chr','track_name') for s in streams]
            streams = [apply(s,'track_name',_f(i)) for i,s in enumerate(streams)]
            s = concatenate(streams, aggregate={'track_name':lambda x:'|'.join(x)})
            s = cobble(s,scored=True)
            name_idx = s.fields.index('track_name')
            start_idx = s.fields.index('start')
            end_idx = s.fields.index('end')
            if 'score' in s.fields:
                score_idx = s.fields.index('score')
                _score = lambda x: x[score_idx]
            else:
                if kw['type']=='tag count': raise ValueError("'score' field not found.")
                _score = lambda x: 0
            for x in s:
                score = _score(x)
                length = x[end_idx]-x[start_idx]
                total_cov += length
                sub = sorted(x[name_idx].split('|'))
                cb = [combinations(sub,k) for k in range(1,len(sub)+1)]
                cb = ['|'.join(sorted(y)) for c in cb for y in c]
                for c in cb:
                    cumcnt[c] += score   # 'cumulative', for the plot
                    cumcov[c] += length
                cnt['|'.join(sub)] += score  # 'separate', for the stats
                cov['|'.join(sub)] += length
        venn_options = {} # tune it here
        output = self.temporary_path(fname='venn_diagram.'+kw['format'])
        legend = None if len(group_names)==len(tracks) \
                 else [os.path.basename(f) for i,f in enumerate(filenames)]
        if kw['type']=='tag count':
            for c in cnt:
                cumcnt[c] = round(cumcnt[c])
                cnt[c] = cnt[c]
            if len(tracks) <= 4:
                venn(cumcnt,legend=legend,options=venn_options,output=output,format=kw['format'])
        elif kw['type']=='coverage %':
            for c in cnt:
                cumcnt[c] = round(cumcov[c]/total_cov * 100)
                cov[c] = cov[c]/total_cov * 100
            if len(tracks) <= 4:
                venn(cumcov,legend=legend,options=venn_options,output=output,format=kw['format'])
        self.new_file(output, 'venn_diagram')
        # Text summary
        output = self.temporary_path(fname='venn_summary.txt')
        with open(output,'wb') as summary:
            summary.write("%s\t%s\t%s\n" % ("Group","Coverage", "Cumulative coverage"))
            if kw['type']=='tag count':
                for c in sorted(cumcnt.keys(), key=lambda x:(len(x),x)):
                    summary.write("%s\t%.2f\t%d\n" % (c,cnt[c],cumcnt[c]))
            elif kw['type']=='coverage %':
                for c in sorted(cumcnt.keys(), key=lambda x:(len(x),x)):
                    summary.write("%s\t%.2f %%\t%d %%\n" % (c,cov[c],cumcov[c]))
        self.new_file(output, 'venn_summary')
        return self.display_time()

