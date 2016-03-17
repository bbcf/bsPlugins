from bsPlugins import *
from bein import execution
from bbcflib.track import track, convert
from bbcflib import genrep, c4seq
import os, tarfile, re


default_path = "/mnt/common/epfl/share"
size_def = 50 #wmax_BRICKS=50
height_def = 500 #wmaxDomainograms=500

# path = '/scratch/cluster/monthly/mleleu/tmp/'
#kw={'sample':path+'segToFrag_FB_HoxD13_all.bedGraph', 'name':'test_domainogram', 'region':'chr2:74450000-74500000', 'wmax_domainograms':'50', 'wmax_BRICKS':'500'}
#segToFrag_KCTD13_CT_all_part.bedGraph chr16:29927651-29928325

class DomainogramsForm(BaseForm):
    sample = twb.BsFileField(label='Input file: ',
        help_text='Select the signal file',
        validator=twb.BsFileFieldValidator(required=True))
    name = twf.TextField(label='Name: ',
        placeholder="(No space or special character)",
        help_text="Prefix name" )
    region = twf.TextField(label='Region: ',
        placeholder="(chr2 or chr2:75450000-75500000)",
        help_text="Chromosome to treat or region to exclude" )
    wmax_BRICKS = twf.TextField(label='window size: ',
        validator=twc.IntValidator(required=True),
        value=size_def,
        help_text='maximum window size (for BRICKS)')
    wmax_domainograms = twf.TextField(label='plot height: ',
        validator=twc.IntValidator(required=True),
        value=height_def,
        help_text='maximum window size (for plots)')

    submit = twf.SubmitButton(id="submit", value="Submit")


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [
        {'id':'sample', 'type':'track', 'required':True, 'label': 'Input file: ', 'help_text': 'Select the signal file' },
        {'id':'name', 'type':'text', 'label': 'Name: ', 'help_text': 'Prefix name', 'placeholder': '(No space or special character)'},
        {'id':'region', 'type':'text', 'label': 'Region: ', 'help_text': 'Chromosome to treat or region to exclude', 'placeholder': '(chr2 or chr2:75450000-75500000)'},
        {'id': 'wmax_BRICKS', 'type': 'int', 'required': True, 'label': 'Window size: ', 'help_text': 'maximum window size (for BRICKS)', 'value': size_def},
        {'id': 'wmax_domainograms', 'type': 'int', 'required': True, 'label': 'Plot height: ', 'help_text': 'maximum window size (for plots)', 'value': height_def}]

out_parameters = [{'id': 'domainograms_tar', 'type': 'file'}]


class DomainogramsPlugin(BasePlugin):
    """Run a domainogram analysis (ref article) on a set of fragments scores (e.g., such as the one obtained with the 4C-seq pipeline)
    """
    info = {
        'title': 'Domainogram analysis',
        'description': __doc__,
        'path': ['Analysis', 'Domainograms'],
#        'output': DomainogramsForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        filename = kw.get('sample')
        assert os.path.exists(str(filename)), "File not found: '%s'" %filename
        script_path = kw.get("script_path",default_path)

        tarname = kw.get('name')+"_domainogram.tar.gz"
        domainograms_tar = tarfile.open(tarname, "w:gz")

        if re.search(r'sql$',str(filename)):
            convert((str(filename),'sql'),(str(filename)+'.bedGraph','bedGraph'))
            filename=str(filename)+'.bedGraph'

        with execution(None) as ex:
            res = c4seq.runDomainogram(ex,infile=filename,name=kw.get('name'),prefix=None,regCoord=kw.get('region'),wmaxDomainogram=str(kw.get('wmaxDomainogram')),wmax_BRICKS=str(kw.get('wmax_BRICKS')),script_path=script_path)

            start = False
            with open(res) as f:
                for s in f:
                    s = s.strip()
                    if re.search('####resfiles####',s):
                        start = True
                    elif start and not re.search("RData",s):
                        domainograms_tar.add(s)
        domainograms_tar.close()

        self.new_file(tarname, 'domainograms_tar')
        return self.display_time()
