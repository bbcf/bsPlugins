from bsPlugins import *
from bbcflib.gfminer.figure import heatmap
from bbcflib.track import track
from numpy import log2, median, array, vstack
from os.path import basename, splitext

nb_colors_def = 10

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'track'},
                 {'id': 'log', 'type': 'boolean'},
                 {'id': 'rowids', 'type': 'boolean'},
                 {'id': 'nb_colors', 'type': 'int'}]
out_parameters = [{'id': 'Heatmap', 'type': 'pdf'},
                  {'id': 'List', 'type': 'file'}]


class HeatmapForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
                            help_text='Enter a table in which the first column contains the IDs and the others the values',
                            validator=twb.BsFileFieldValidator(required=True))
    log = twf.CheckBox(label='Log: ',
                       value=False,
                       help_text='Take the log2(1+x) of each value x')
    rowids = twf.CheckBox(label='List: ',
                          value=False,
                          help_text='Number the rows in the heatmap and make a file with the corresponding IDs')
    nb_colors = twf.TextField(label='Number of colors: ',
                              validator=twc.IntValidator(required=False),
                              value=nb_colors_def,
                              help_text='Number of colors between blue and red (default: 10)')
    submit = twf.SubmitButton(id="submit", value="Plot")


class HeatmapPlugin(BasePlugin):
    """Creates a heatmap of the table using *rows* as row labels and *columns* as column labels.
The values are assumed to be equal to log2(raw data). If not, you have to select the option *Log*.
Selecting the option *List* will print numbers beside the rows in the heatmap and make a file with the corresponding IDs."""
    info = {
        'title': 'Make a heatmap of a numeric table',
        'description': __doc__,
        'path': ['Graphics', 'Heatmap'],
        'output': HeatmapForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        table = open(kw.get('table'))
        title = splitext(basename(table))[0]
        nb_colors = int(kw.get('nb_colors') or nb_colors_def)
        logscale = kw.get('log',False)
        if isinstance(logscale, basestring):
            logscale = (linscale.lower() in ['1', 'true', 't','on'])
        make_list = kw.get('rowids',False)
        if isinstance(make_list, basestring):
            make_list = (make_list.lower() in ['1', 'true', 't','on'])
        names = []
        values = None
        for line in table:
            newline = line.strip("\n\r").split("\t")
            names.append(newline.pop(0))
            if values is None: values = array(newline,dtype=float)
            else: values = vstack([values,array(newline,dtype=float)])
        if logscale: values = log2(1+values)
        values -= median(values,axis=0)
        values /= median(abs(values),axis=0)
        pdf = self.temporary_path(fname='%s.pdf' %title)
        pdf, roword = heatmap(values, output=pdf, rows=names, main=title, nb_colors=nb_colors, return_rowInd=True)
        if make_list:
            List_ID = self.temporary_path(fname='%s_row_order.txt' %title)
            with open(List_ID,"w") as L:
                L.write("\n".join(["%i\t%s" %(n,names[n-1]) for n in roword]))
            self.new_file(List_ID, 'List')
        self.new_file(pdf, 'Heatmap')
        return self.display_time()
