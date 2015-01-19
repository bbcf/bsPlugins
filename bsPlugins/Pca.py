from bsPlugins import *
import os
import subprocess

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'table', 'type': 'txt', 'required': True},
                 {'id': 'columns', 'type': 'text', 'required': True},]
out_parameters = [{'id': 'pca_biplot', 'type': 'pdf'}]

class PcaForm(BaseForm):
    table = twb.BsFileField(label='Table: ',
        help_text='Select scores table',
        validator=twb.BsFileFieldValidator(required=True))
    columns = twf.TextField(label='Columns selection: ',
        value="all",
        help_text='Which columns? "all"/"rpkm" (name contains rpkm)/"1,3,5,6,..."')
    submit = twf.SubmitButton(id="submit", value="Plot")


class PcaPlugin(BasePlugin):
    """Generates a PCA biplot for diagnostic."""
    info = {'title': 'PCA',
            'description': __doc__,
            'path': ['Graphics', 'PCA'],
            'output': PcaForm,
            'in': in_parameters,
            'out': out_parameters,
            'meta': meta}

    def __call__(self, **kw):
        filename = kw.get('table')
        assert os.path.exists(str(filename)), "File not found: '%s'" % filename
        subprocess.call(['pca.R', filename, "pca_biplot", kw['columns']])
        self.new_file("pca_biplot.pdf", "pca_biplot")
        return self.display_time()
