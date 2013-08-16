from bsPlugins import *
from bbcflib.track import track, FeatureStream
from bbcflib.gfminer.common import cobble
from bbcflib.gfminer.stream import concatenate
from bbcflib.gfminer.figure import venn
from bbcflib import genrep
from itertools import combinations
import os

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id':'files', 'type':'track', 'required':True, 'multiple':'SigMulti'},
        {'id':'type', 'type':'list'},
        {'id':'format', 'type':'list'},
]
out_parameters = [{'id':'venn_diagram', 'type':'file'},
                  {'id':'venn_summary', 'type':'file'}]



class VennDiagramForm(BaseForm):
    child = twd.HidingTableLayout()
    class SigMulti(Multi):
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
    submit = twf.SubmitButton(id="submit", value="Submit")


class VennDiagramPlugin(BasePlugin):
    """Creates a Venn diagram of the proportions of
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
        self.debug(**kw)
        filenames = kw['SigMulti']['files']
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
        def _add_label(s,x):
            _f = s.fields+['track_name']
            return FeatureStream((y+(x,) for y in s), fields=_f)
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
                summary.write("Empty content (no coverage). Check the chromosome names.")
            self.new_file(output, 'venn_summary')
            return
        venn_options = {} # tune it here
        format = kw.get('format','pdf')
        output = self.temporary_path(fname='venn_diagram.'+format)
        legend = [t.name for t in tracks]
        if _scored:
            for c in combn:
                c2[c] = round(c2[c])
        else:
            for c in combn:
                c2[c] = round((100*c2[c])/total_cov)
                c1[c] = (100*c1[c])/total_cov
        if len(tracks) <= 4:
            venn(c2,legend=legend,options=venn_options,output=output,format=format)
        self.new_file(output, 'venn_diagram')
        # Text summary
        output = self.temporary_path(fname='venn_summary.txt')
        with open(output,'wb') as summary:
            summary.write("%s\t%s\t%s\n" % ("Group","Coverage", "Cumulative coverage"))
            record = "%s\t%.2f\t%d\n"
            for c in sorted(combn, key=lambda x:(len(x),x)):
                summary.write(record%(c,c1[c],c2[c]))
        self.new_file(output, 'venn_summary')
        return self.display_time()

