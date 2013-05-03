from bsPlugins import *
from bein import execution
from bbcflib import genrep
from bbcflib.common import fasta_composition
from bbcflib.motif import save_motif_profile
from bbcflib.btrack import track, FeatureStream
import re, os


g = genrep.GenRep()
input_types = [(0, 'Fasta upload'), (1, 'Select regions from genome')]
input_map = {0: ['fastafile'], 1: ['assembly', 'regions']}

class MotifScanForm(BaseForm):

    motif_list = g.motifs_available()
    assembly_list = g.assemblies_available()

    child = twd.HidingTableLayout()
    input_type = twd.HidingRadioButtonList(label='Sequence source: ',
                                           options=input_types,
                                           value=0,
                                           mapping=input_map)
    s1 = twf.Spacer()
    fastafile = twb.BsFileField(label='Fasta file: ', help_text='Sequences to scan')
    assembly = twf.SingleSelectField(label='Assembly: ', options=assembly_list,
                                     prompt_text=None,
                                     help_text='Assembly to fetch sequences from')
    regions = twb.BsFileField(label='Regions: ', help_text='Genomic regions to scan (e.g. bed)')
    s2 = twf.Spacer()
    background = twb.BsFileField(label='Background: ',
                               help_text='File of background frequencies (default: genome-wide frequencies)')
    motifs = twf.MultipleSelectField(label='Motifs: ', options=motif_list,
                                     help_text='Select motifs to be scanned')
    customMotif = twb.BsFileField(label='Custom motif: ',
                                help_text='An optional custom/additional motif to scan (.mat)')
    threshold = twf.TextField(label='Threshold: ', value='0.0')
    submit = twf.SubmitButton(id="submit", value='Scan sequences')


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'input_type', 'type': 'radio'},
                 {'id': 'fastafile', 'type': 'userfile'},
                 {'id': 'background', 'type': 'txt'},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'regions', 'type': 'track'},
                 {'id': 'motifs', 'type': 'list'},
                 {'id': 'customMotif', 'type': 'txt'},
                 {'id': 'threshold', 'type': 'float', 'required': True}]
out_parameters = [{'id': 'motif_track', 'type': 'track'}]

class MotifScanPlugin(BasePlugin):
    info = {
        'title': 'Motif scanner',
        'description': 'Scan motifs PWM on a set of a sequences',
        'path': ['Sequence analysis', 'Motif scanner'],
        'output': MotifScanForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }


    def __call__(self, **kw):
        fasta_file = kw.get('fastafile')
        background = kw.get('background') or None
        assembly_id = kw.get('assembly') or None
        regions_file = kw.get('regions')
        motifs_list = kw.get('motifs')
        motif_add = kw.get('customMotif')
        threshold = float(kw.get('threshold') or 0)

        if motifs_list is None: motifs_list = []
        if not isinstance(motifs_list, list): motifs_list = [motifs_list]

        if background is None and assembly_id is None:
            background = self.temporary_path(fname='background.txt')
            stats = {'A': 0.25,'C': 0.25, 'G': 0.25, 'T': 0.25}
            if fasta_file:
                with execution(None) as ex:
                    stats = fasta_composition(ex,fasta_file,frequency=True)
            with open(background,"w") as bgr:
                bgr.write(" ".join(["1"]+[str(stats[n]) for n in ['A','C','G','T']]))
        if assembly_id is not None:
            assembly = genrep.Assembly(assembly_id)
        else:
            if regions_file is not None:
                raise ValueError("Please specify an assembly if you specify regions.")
            assembly = None

        motifs = {}
        if motif_add is not None:
            mname = os.path.basename(os.path.splitext(motif_add)[0])
            if mname: motifs[mname] = motif_add
        for mot in motifs_list:
            gid, mname = mot.split(' ')
            pwmfile = self.temporary_path()
            _ = g.get_motif_PWM(int(gid), mname, output=pwmfile)
            motifs[mname] = pwmfile

        if len(motifs) == 0:
            raise ValueError("Please give at least one motif to scan for")

        track_output = self.temporary_path(fname='motif_scan', ext="sql")
        with execution(None) as ex:
            _ = save_motif_profile( ex, motifs, assembly, regions_file, fasta_file,
                                    background=background, threshold=threshold,
                                    output = track_output,
                                    description=None, via='local' )
        self.new_file(track_output, 'motif_track')
        return self.display_time()



