from bsPlugins import *
from bbcflib.gfminer.figure import heatmap
from bbcflib.track import track
from numpy import array, log2, median, abs

nb_colors_def = 10

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'track'},
                 {'id': 'log', 'type': 'boolean'},
                 {'id': 'list', 'type': 'boolean'},
                 {'id': 'nb_colors', 'type': 'int'}]
out_parameters = [{'id': 'Heatmap', 'type': 'pdf'},
                  {'id': 'List', 'type': 'txt'}]


class HeatmapForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
                                help_text='Enter a table in which the first column contains the IDs and the others the values',
                                validator=twb.BsFileFieldValidator(required=True))
    log = twf.CheckBox(label='Log: ',
                                value=False,
                                help_text='Take the log2(1+x) of each value x')
    list = twf.CheckBox(label='List: ',
                                value=False,
                                help_text='Number the rows in the heatmap and make a file with the corresponding IDs')
    nb_colors = twf.TextField(label='Number of colors: ',
                                validator=twc.IntValidator(required=False),
                                value=nb_colors_def,
                                help_text='Number of colors between blue and red (default: 10)')
    submit = twf.SubmitButton(id="submit", value="Plot")


class HeatmapPlugin(BasePlugin):
    """Creates a heatmap of the table using *rows* as row labels and *columns* as column labels. The values are assumed to be equal to log2(raw data). If not, you have to select the option *Log*. Selecting the option *List* will print numbers beside the rows in the heatmap and make a file with the corresponding IDs."""
    info = {
        'title': 'Make heatmap of the table',
        'description': __doc__,
        'path': ['Graphics', 'Heatmap'],
        'output': HeatmapForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        table = open(kw.get('table'))
        names = []; values = []
        for line in table:
            newline = line.strip("\n").strip("\r").split("\t")
            name = newline[0]
            names.append(name)
            value = [float(newline[i]) for i in range(1,len(newline))]
            values.append(value)
        logscale = kw.get('log',False)
        if isinstance(logscale, basestring):
            logscale = (linscale.lower() in ['1', 'true', 't','on'])
        if logscale:
            for i in range(0,len(values)):
                values[i] = log2([x+1 for x in values[i]])
        values = array(values)
        med = median(values,axis=0)
        mad = median(abs(values-med),axis=0)
        values = ((values-med)/mad)
        nb_colors = int(kw.get('nb_colors') or nb_colors_def)
        #names = ['']*len(data)
        title = "List of genes"
        pdf = self.temporary_path(fname='Heatmap.pdf')
        heatmap(array(values), output=pdf, rows=names, main=title, nb_colors=nb_colors)
        self.new_file(pdf, 'Heatmap')
        return self.display_time()
