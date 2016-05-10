from bsPlugins import *
import re
from bein import execution
from bbcflib.gfminer.bedtools import *

all_tools = ["annotateBed", "bamToBed", "bamToFastq", "bed12ToBed6",
             "bedpeToBam", "bedToBam", "bedToIgv", "closestBed",
             "clusterBed", "complementBed", "coverageBed", "expandCols",
             "fastaFromBed", "flankBed", "genomeCoverageBed", "getOverlap",
             "groupBy", "intersectBed", "linksBed", "mapBed", "maskFastaFromBed",
             "mergeBed", "multiBamCov", "multiIntersectBed", "nucBed",
             "pairToBed", "pairToPair", "randomBed", "shuffleBed",
             "slopBed", "sortBed", "subtractBed", "tagBam",
             "unionBedGraphs", "windowBed", "windowMaker"]

tools_map = {
    0: ["bedfile", "files"],
    1: ["bamfile"],
    2: ["bamfile"],
    3: ["bedfile"],
    4: ["bedfile", "genomefile"],
    5: ["bedfile", "genomefile"],
    6: ["bedfile"],
    7: ["afile", "bfile"],
    8: ["bedfile"],
    9: ["bedfile", "genomefile"],
    10: ["bfile", "afile", "bamfile"],
    11: ["bedfile", "column"],
    12: ["bedfile", "fastafile"],
    13: ["bedfile", "genomefile"],
    14: ["genomefile", "bedfile", "bamfile"],
    15: ["bedfile"],
    16: ["bedfile", "groupcol", "opcol", "operation"],
    17: ["bfile", "afile", "bamfile"],
    18: ["bedfile"],
    19: ["afile", "bfile"],
    20: ["bedfile", "fastafile"],
    21: ["bedfile"],
    22: ["bedfile", "bamfiles"],
    23: ["bedfiles"],
    24: ["bedfile", "fastafile"],
    25: ["bfile", "afile", "bamfile"],
    26: ["afile", "bfile"],
    27: ["genomefile"],
    28: ["bedfile", "genomefile"],
    29: ["bedfile", "genomefile"],
    30: ["bedfile"],
    31: ["afile", "bfile"],
    32: ["bedfiles", "labels", "bamfile"],
    33: ["files"],
    34: ["bfile", "afile", "bamfile"],
    35: ["bedfile", "genomefile"]}

all_params = dict((y, '') for x in tools_map.values() for y in x).keys()

file_params = {"simple": [x for x in all_params if x[-4:] == "file"],
               "multiple": [x for x in all_params if x[-5:] == "files"]}

all_file_params = [{'id': x[0:], 'type': 'track', 'multiple': x} for x in all_params if x[-5:] == 'files']+\
                  [{'id': x, 'type': 'track'} for x in all_params if x[-4:] == 'file']

gr_operations = ["sum", "count", "count_distinct", "min", "max",
                 "mean", "median", "mode", "antimode", "stdev", "sstdev",
                 "collapse", "distinct", "concat", "freqdesc", "freqasc"]


other_params = [
    {'id': 'tool', 'type': 'list', 'label': 'Tool: ', 'help_text': 'Select BedTool', 'value': 0, 'options':list(enumerate(all_tools)), 'mapping': tools_map},
    {'id': 'column', 'type': 'int', 'required': True, 'label': 'column: ', 'value': 1},
    {'id': 'groupcol', 'type': 'int', 'required': True, 'label': 'groupcol: ', 'value': 1},
    {'id': 'opcol', 'type': 'int', 'required': True, 'label': 'opcol: ', 'value': 1},
    {'id': 'labels', 'type': 'text', 'label': 'labels: '},
    {'id': 'operation', 'type': 'list', 'label': 'operations: ', 'options': gr_operations},
    {'id': 'useropts', 'type': 'text', 'label': 'Options: ', 'help_text': 'Additional command-line options'}]



class BedToolsForm(BaseForm):
    child = twd.HidingTableLayout()
    tool = twd.HidingSingleSelectField(label='Tool: ',
                                       prompt_text=None,
                                       options=list(enumerate(all_tools)),
                                       mapping=tools_map,
                                       value=0,
                                       help_text='Select BedTool')
    for it in file_params['simple']:
        vars()[it] = twb.BsFileField(label=it+': ', validator=twb.BsFileFieldValidator())

    class Mfiles(twb.BsMultiple):
        label = 'files: '
        files = twb.BsFileField(label='', validator=twb.BsFileFieldValidator())

    class Mbamfiles(twb.BsMultiple):
        label = 'bamfiles: '
        bamfiles = twb.BsFileField(label=' ', validator=twb.BsFileFieldValidator())

    class Mbedfiles(twb.BsMultiple):
        label = 'bedfiles: '
        bedfiles = twb.BsFileField(label=' ', validator=twb.BsFileFieldValidator())

    column = twf.TextField(label='column: ',
                           validator=twc.IntValidator(min=1, max=100, required=True),
                           value=1)
    groupcol = twf.TextField(label='groupcol: ',
                             validator=twc.IntValidator(min=1, max=100, required=True),
                             value=1)
    opcol = twf.TextField(label='opcol: ',
                          validator=twc.IntValidator(min=1, max=100, required=True),
                          value=1)
    labels = twf.TextField(label='labels: ')
    operation = twf.SingleSelectField(label='operation: ',
                                      options=gr_operations, prompt_text=None)
    useropts = twf.TextField(label='Options: ',
                             help_text='Additional command-line options')
    submit = twf.SubmitButton(id="submit", value="Run")

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = all_file_params + other_params

out_parameters = [{'id': tool+'_result', 'type': 'track'} for tool in all_tools]


class BedToolsPlugin(BasePlugin):
    """Bedtools collection."""
    info = {
        'title': 'BedTools',
        'description': __doc__,
        'path': ['Intervals', 'Bedtools'],
#        'output': BedToolsForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        _tool = kw.pop('tool')
        try:
            selected_tool = all_tools[int(_tool)]
        except ValueError:
            selected_tool = str(_tool)
        _toolid = all_tools.index(selected_tool)
        kw = dict((k,v) for k,v in kw.iteritems() 
                  if k in ['outfile']+tools_map[_toolid])
        kw['outfile'] = self.temporary_path(fname=selected_tool+'.txt')
        if kw.get('useropts'):
            reo = re.search(r'([\w\s\-,.=]+)', kw.pop('useropts'))
            if reo:
                key = None
                for x in reo.groups()[0].split():
                    if x.startswith('-'):
                        key = x
                        kw[key] = ''
                    elif key:
                        kw[key] = str(x)
        if kw.get('labels'):
            kw['labels'] = kw['labels'].split(",")
        for x in all_params:
            if x[-5:] == "files" and kw.get(x):
                kw[x[1:]] = kw.pop(x)[x[1:]]
        for k in kw.keys():
            if kw[k] in (None,'',u'',[],{}):
                kw.pop(k)
        with execution(None) as ex:
            output = eval(selected_tool)(ex, **kw)
        self.new_file(output, selected_tool+'_result')
        return self.display_time()

