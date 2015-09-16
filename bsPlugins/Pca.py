from bsPlugins import *
import os
import subprocess

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'txt', 'required': True, 'label': 'Table: ', 'help_text':'Select scores table'},
                 {'id': 'columns', 'type': 'text', 'required': True, 'label': 'Columns selection: ', 'help_text':'Which columns? "all"/"1,3,5,6"', 'value': 'all'}]
out_parameters = [{'id': 'pca_biplot', 'type': 'pdf'}]

class PcaForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
        help_text='Select scores table',
        validator=twb.BsFileFieldValidator(required=True))
    columns = twf.TextField(label='Columns selection: ',
        value="all",
        help_text='Which columns? "all"/"1,3,5,6"')
    submit = twf.SubmitButton(id="submit", value="Plot")


class PcaPlugin(BasePlugin):
    """Generates a PCA biplot for diagnostic.

  The intput is a tab-delimited table with feature names in the first column,
  then one column of respective scores per sample.
  The first line is a header of the type "ID  sample1  sample2 ...".
  Ideal relevant column names are of the form
  "[prefix.]<group_name>.<#replicate>", e.g. "rpkm.Control.2" or "TRF2_KO.1",
  so that the plugin will recognize the group name and replicate number:
  all samples with a different group name will get a different color.

  Example of input:

    ID    rpkm.KO.1    rpkm.KO.2    rpkm.WT.1
    AAA   9.0          8.1          0.0
    BBB   3.2          9.0          7.4

    """
    info = {'title': 'PCA',
            'description': __doc__,
            'path': ['Graphics', 'PCA'],
#            'output': PcaForm,
            'in': in_parameters,
            'out': out_parameters,
            'meta': meta}

    def __call__(self, **kw):
        filename = kw.get('table')
        assert os.path.exists(str(filename)), "File not found: '%s'" % filename
        subprocess.call(['pca.R', filename, "pca_biplot", kw['columns']])
        self.new_file("pca_biplot.pdf", "pca_biplot")
        return self.display_time()
