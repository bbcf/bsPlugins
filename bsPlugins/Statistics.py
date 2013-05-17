from bsPlugins import *
from bbcflib.btrack import track,stats
from bbcflib import genrep
import os


class StatisticsForm(BaseForm):
    child = twd.HidingTableLayout()
    sample = twb.BsFileField(label='Input file: ',
        help_text='Select the file to examine',
        validator=twb.BsFileFieldValidator(required=True))
    #assembly = twf.SingleSelectField(label='Assembly: ',
    #    prompt_text=None,
    #    options=genrep.GenRep().assemblies_available(),
    #    help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id':'sample', 'type':'track', 'required':True},
        #{'id':'assembly', 'type':'assembly'},
]
out_parameters = [{'id':'stats', 'type':'file'},
                  {'id':'density_plot', 'type':'file'}]


class StatisticsPlugin(BasePlugin):
    description = """Calculates diverse statistics from a track file,
    such as a distribution of scores and feature lengths, and prints them to the output file.
    """
    info = {
        'title': 'Statistics',
        'description': description,
        'path': ['Files', 'Statistics'],
        'output': StatisticsForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        self.debug(**kw)
        #assembly = genrep.Assembly(kw['assembly'])
        sample = kw['sample']
        fname = os.path.splitext(os.path.basename(sample))[0]
        output = self.temporary_path(fname=fname+'_stats.txt')
        with open(output,'wb') as out:
            stats(sample,out=out)
        self.new_file(output, 'stats')
        return self.display_time()

