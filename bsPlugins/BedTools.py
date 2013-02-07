from . import *
import re
from bein import execution

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
    0: ["bedfile", "AllFilesMulti.files"],
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
    22: ["bedfile", "AllFilesMulti.bamfiles"],
    23: ["AllFilesMulti.bedfiles"],
    24: ["bedfile", "fastafile"],
    25: ["bfile", "afile", "bamfile"],
    26: ["afile", "bfile"],
    27: ["genomefile"],
    28: ["bedfile", "genomefile"],
    29: ["bedfile", "genomefile"],
    30: ["bedfile"],
    31: ["afile", "bfile"],
    32: ["AllFilesMulti.bedfiles", "labels", "bamfile"],
    33: ["AllFilesMulti.files"],
    34: ["bfile", "afile", "bamfile"],
    35: ["bedfile", "genomefile"]}

all_params = dict((y, '') for x in tools_map.values() for y in x).keys()

file_params = {"simple": [x for x in all_params if x[-4:] == "file"],
               "multiple": [x for x in all_params if x[-5:] == "files"]}

all_file_params = [{'id': x[14:], 'type': 'file', 'multiple': True}
                   for x in all_params if x[-5:] == 'files']+\
                   [{'id': x, 'type': 'file', 'multiple': False}
                    for x in all_params if x[-4:] == 'file']

other_params = [
    {'id': 'tool', 'type': 'list'},
    {'id': 'column', 'type': 'int', 'required': True},
    {'id': 'groupcol', 'type': 'int', 'required': True},
    {'id': 'opcol', 'type': 'int', 'required': True},
    {'id': 'labels', 'type': 'text'},
    {'id': 'operation', 'type': 'list'},
    {'id': 'useropts', 'type': 'text'},
]


gr_operations = ["sum", "count", "count_distinct", "min", "max",
                 "mean", "median", "mode", "antimode", "stdev", "sstdev",
                 "collapse", "distinct", "concat", "freqdesc", "freqasc"]


class BedToolsForm(BaseForm):
    child = twd.HidingTableLayout()
    tool = twd.HidingSingleSelectField(label_text='Tool: ',
        prompt_text=None,
        options=list(enumerate(all_tools)),
        mapping=tools_map,
        value=0,
        help_text='Select BedTool')
    for it in file_params['simple']:
        vars()[it] = twf.FileField(label=it + ': ', validator=twf.FileValidator())
    class AllMultiFiles(Multi):
        for it in file_params['multiple']:
            vars()[it[14:]] = twf.FileField(label=it[14:]+': ', validator=twf.FileValidator())
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
    operation = twf.SingleSelectField(label_text='operation: ',
        options=gr_operations, prompt_text=None)
    useropts = twf.TextField(label_text='Options: ',
        help_text='Additional command-line options')
    submit = twf.SubmitButton(id="submit", value="Run")

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = all_file_params + other_params

out_parameters = [{'id': 'bedtools_result', 'type': 'file'}]


class BedToolsPlugin(OperationPlugin):

    info = {
        'title': 'BedTools',
        'description': 'Bedtools collection',
        'path': ['Features', 'Bedtools'],
        'output': BedToolsForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        kw['outfile'] = self.temporary_path()
        reo = re.search(r'([\w\s\-,.=]+)', kw.pop('useropts') or '')
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
        with execution(None) as ex:
            output = eval(all_tools[kw.pop('tool')])(ex, **kw)
        self.new_file(output, 'bedtools_result')
        return 1
