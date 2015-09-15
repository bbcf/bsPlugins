from bsPlugins import *
from bbcflib.gfminer.figure import heatmap
from bbcflib.track import track
from numpy import log2, median, array, vstack, nan
from os.path import basename, splitext

nb_colors_def = 10

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'file', 'required': True, 'label': 'Table: ', 'help_text': 'Enter a table in which the first column contains the IDs and the others the values'},
                 {'id': 'log', 'type': 'boolean', 'label': 'Log: ', 'help_text': 'Take the log2(1+x) of each value x', 'value': False},
                 {'id': 'cor', 'type': 'boolean', 'label': 'Correlation: ', 'help_text': 'Use the correlation as the distance function to build the clusters', 'value': False},
                 {'id': 'rowids', 'type': 'boolean', 'label': 'List: ', 'help_text': 'Number the rows in the heatmap and make a file with the corresponding IDs', 'value': False},
                 {'id': 'nb_colors', 'type': 'int', 'label': 'Number of colors: ', 'help_text': 'Number of colors between blue and red (default: 10)', 'value': nb_colors_def}]
out_parameters = [{'id': 'Heatmap', 'type': 'pdf'},
                  {'id': 'List', 'type': 'file'}]


class HeatmapForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
                            help_text='Enter a table in which the first column contains the IDs and the others the values',
                            validator=twb.BsFileFieldValidator(required=True))
    log = twf.CheckBox(label='Log: ',
                       value=False,
                       help_text='Take the log2(1+x) of each value x')
    cor = twf.CheckBox(label='Correlation: ',
                       value=False,
                       help_text='Use the correlation as the distance function to build the clusters')
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
#       'output': HeatmapForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        table = kw.get('table')
        title = splitext(basename(table))[0]
        nb_colors = int(kw.get('nb_colors') or nb_colors_def)
        logscale = kw.get('log',False)
        if isinstance(logscale, basestring):
            logscale = (logscale.lower() in ['1', 'true', 't','on'])
        cor = kw.get('cor',False)
        if isinstance(cor, basestring):
            cor = (cor.lower() in ['1', 'true', 't','on'])
        if cor: cor = True
        else: cor = False
        make_list = kw.get('rowids',False)
        if isinstance(make_list, basestring):
            make_list = (make_list.lower() in ['1', 'true', 't','on'])
        names = []
        values = None
        with open(table) as _tabl:
            col_labels = _tabl.readline().strip("\n\r").split("\t")[1:]
            for line in _tabl:
                newline = line.strip("\n\r").split("\t")
                names.append(newline.pop(0))
                if values is None: values = array(newline,dtype=float)
                else: values = vstack([values,array(newline,dtype=float)])
        if logscale: values = log2(1+values)
        values -= median(values,axis=0)
        mads = median(abs(values),axis=0)
        mads[mads < 1e-6] = nan
        values /= mads
        pdf = self.temporary_path(fname='%s.pdf' %title)
        if make_list:
            pdf, roword = heatmap(values, output=pdf, rows=names, columns=col_labels, main=title, nb_colors=nb_colors, return_rowInd=True, cor=cor)
            List_ID = self.temporary_path(fname='%s_row_order.txt' %title)
            with open(List_ID,"w") as L:
                L.write("Row"+"\t"+"ID.number"+"\t"+"ID.name"+"\n")
                L.write("\n".join(["%i\t%i\t%s" %(i+1,n,names[n-1]) for i, n in enumerate(roword)]))
            self.new_file(List_ID, 'List')
        else:
            pdf = heatmap(values, output=pdf, rows=names, columns=col_labels, main=title, nb_colors=nb_colors, return_rowInd=False, cor=cor)
        self.new_file(pdf, 'Heatmap')
        return self.display_time()
