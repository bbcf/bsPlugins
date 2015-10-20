from bsPlugins import *
from bbcflib.track import track
from bbcflib import genrep
import os, shutil
from bbcflib.maplot import MAplot

__requires__ = ["numpy"]

ftypes = [(0, 'genes bodies'), (1, 'gene promoters'), (2, 'exons'), (3, 'custom upload')]
prom_up_def = 1000
prom_down_def = 100

input_opts=['Table', 'Signals']
input_map={'Table': ['table'],
        'Signals': ['Group1','Group2','feature_type','assembly']}
f_map={ftypes[-1][0]: ['features'],
                 1: ['upstream', 'downstream']}

class MaplotForm(BaseForm):
    child = twd.HidingTableLayout()
    input_type = twd.HidingRadioButtonList(label='Input type: ',
                                           options=('Table', 'Signals'),
                                           mapping={'Table':  ['table'],
                                                    'Signals': ['Group1','Group2','feature_type','assembly'],},
                                           help_text='Select input type (Formatted table, or signal tracks)')
    table = twb.BsFileField(label='Table: ',
        help_text='Select scores table',
        validator=twb.BsFileFieldValidator(required=True))
    class Group1(twb.BsMultiple):
        label = "Signals group 1: "
        signals1 = twb.BsFileField(label=' ',
            help_text='Select signal files (position and score, e.g. bedgraph)',
            validator=twb.BsFileFieldValidator(required=True))
    class Group2(twb.BsMultiple):
        label = "Signals group 2: "
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
        validator=twb.BsFileFieldValidator())
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
        {'id': 'input_type', 'type': 'radio', 'label': 'Input type: ','help_text': 'Select input type (Formatted table, or signal tracks)', 'options': input_opts, 'value': 'Table', 'mapping': input_map },
        {'id': 'signals1', 'type': 'track', 'required': True, 'multiple': 'Group1', 'label': 'Signals group 1: ', 'help_text': 'Select signal files (position and score, e.g. bedgraph)'},
        {'id': 'signals2', 'type': 'track', 'required': True, 'multiple': 'Group2', 'label': 'Signals group 2: ', 'help_text': 'Select signal files (position and score, e.g. bedgraph)'},
        {'id': 'table', 'type': 'txt', 'required': True, 'lable': 'Table: ', 'help_text': 'Select scores table'},
        {'id': 'feature_type', 'type': 'list', 'required': True, 'label': 'Feature type: ', 'help_text': 'Choose a feature set or upload your own', 'options': ftypes, 'prompt_text': None, 'mapping': f_map},
        {'id': 'upstream', 'type': 'int', 'required': True, 'label': 'Promoter upstream distance: ', 'help_text': 'Size of promoter upstream of TSS', 'value': prom_up_def},
        {'id': 'downstream', 'type': 'int', 'required': True, 'label': 'Promoter downstream distance: ', 'help_text': 'Size of promoter downstream of TSS', 'value': prom_down_def},
        {'id': 'assembly', 'type': 'assembly', 'required': True, 'label': 'Assembly: ', 'help_text': 'Reference genome','options': genrep.GenRep().assemblies_available()},
        {'id': 'features', 'type': 'track', 'required': True, 'label': 'Custom feature set: ', 'help_text': 'Select a feature file (e.g. bed)'}]

out_parameters = [{'id': 'MA-plot', 'type': 'file'}]


class MaplotPlugin(BasePlugin):
    """Creates an MA-plot to compare levels of expression of genomic features
across two samples.

The input can be of two different types:

* Two 'signal' files, i.e. bedGraph-type text files,
  and a list of genomic features - either from a pre-defined list such as Ensembl genes,
  or a custom bed-like file. The name of each sample is the one given in the track
  definition line ("track name=... description=... etc."), if specified, otherwise the name of
  the file (without extension). </li>
* A tab-delimited table with feature names in the first column, then one column of respective
  scores per sample. The first line is a header of the type "id  sample1  sample2 ...".
  """
    info = {
        'title': 'MA-plot',
        'description': __doc__,
        'path': ['Graphics', 'MA-plot'],
#        'output': MaplotForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):

        if kw.get('input_type') == 'Table':
            table = kw.get('table')
            assert os.path.exists(str(table)), "File not found: '%s'" % table
            with open(table) as t:
                colnames = t.readline()
                _f = colnames.strip().split()
                nscores = len(_f)-1
            groups = len(list(set([x.split('.')[0] for x in _f])))
            if nscores == 2: # 3 columns, cols 2 and 3 contain the scores
                sample1 = [2]
                sample2 = [3]
            elif len(groups) == 2: # more columns, look if there are two groups of prefixes
                sample1 = [_f.index(x) for x in _f if x.split('.')==groups[0]]
                sample2 = [_f.index(x) for x in _f if x.split('.')==groups[1]]
            else: # not implemented yet, ask the user to choose the columns he wants? Checkboxes...
                raise ValueError("For the moment, either have only 2 columns of scores, \
                                 or use names of the form <group_name>.<run_id>")
        else:
            # Use QuantifyTablePlugin to build a table from score tracks
            from QuantifyTable import QuantifyTablePlugin
            # Set QuantifyTablePlugin options
            kw['score_op'] = 'sum'
            kw['format'] = 'txt'
            #signals1 = kw['Group1']['signals1']
            signals1 = kw['signals1']
            #signals2 = kw['Group2']['signals2']
            signals2 = kw['signals2']
            if not isinstance(signals1,(list,tuple)): signals1 = [signals1]
            if not isinstance(signals2,(list,tuple)): signals2 = [signals2]
            kw['signals'] = signals1 + signals2
            signals = kw['signals']
            nscores = len(signals)
            qtable = QuantifyTablePlugin().quantify(**kw)
            # Remove useless fields and add header based on file names
            qtable = track(qtable, format='txt', fields=['chr','start','end','name']+ \
                                                        ['score'+str(i) for i in range(nscores)])
            table = self.temporary_path('scores_table.txt')
            _f = ['score'+str(i) for i in range(nscores)]
            strack = track(table, fields=['name']+_f)
            signal_tracks = [track(s) for s in signals]
            signames = [s.name for s in signal_tracks]
            strack.write([('Name',signames[0],signames[1])])
            strack.write(qtable.read(fields=strack.fields))
            sample1 = range(len(signals1))
            sample2 = range(nscores-len(signals1))

        output_filename = MAplot(table, cols={1:sample1, 2:sample2})
        output = self.temporary_path(fname='maplot.png')
        shutil.move(output_filename,output)
        self.new_file(output, 'MA-plot')
        return self.display_time()

# nosetests --logging-filter=-tw2 test_Maplot.py
