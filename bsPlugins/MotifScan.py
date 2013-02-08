from bsPlugins import *
from bein import execution
from bbcflib import genrep
from bbcflib.motif import save_motif_profile
from bbcflib.btrack import track, FeatureStream
import re


input_types = [(0, 'Custom upload'), (1, 'From assembly')]
input_map = {input_types[0][1]: ['fastafile', 'background'],
             input_types[0][1]: ['assembly', 'regions']}

class MotifForm(DynForm):
    genomes = g.get_genrep_objects('genomes', 'genome') ####????
    motif_list = []
    for genome in genomes:
        if genome.motif_matrix_url != None and genome.motif_matrix_url != 'null':
            curMotifs = g.get_genrep_objects('genomes/' + str(genome.id) + '/get_matrix', 'motif') ####??
            for motif in curMotifs:
                motif_list.append((str(genome.id) + " " + motif.name, genome.name + " - " + motif.name))
    motif_list.sort(key=lambda x: x[1])

    #jQuery's already on the page
    twj.jquery_js.no_inject = True

    sequenceSource = twd.HidingRadioButtonList(label='Sequence source', 
                                               options=input_types, mapping=input_map)
    _ = twf.Spacer()
    fastafile = twf.FileField(label='Fasta file: ', help_text='Sequences to scan')
#### TextField?
    background = twf.FileField(label='Background: ',
                               help_text='File of background frequencies (default: genome-wide frequencies)')
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=genrep.GenRep().assemblies_available(),
                                     help_text='Assembly to fetch sequences from')
    regions = twf.FileField(label='Regions: ', help_text='Genomic regions to scan (e.g. bed)')
    _ = twf.Spacer()
    motifs = twjqSelect.Select2MultipleSelectField(placeholder='Select the motifs to scan for', 
                                                   options=motif_list,
                                                   help_text='')
    customMotif = twf.FileField(label='Custom motif: ',
                                help_text='An optional custom/additional motif to scan (.mat)')
    threshold = twf.TextField(label='Threshold: ',
                              value='0.0',
                              validator=twc.RegexValidator(required=True, regex=re.compile('^\d+(\.\d+)?$')))
    submit = twf.SubmitButton(id="submit", value='Scan sequences')


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'sequenceSource', 'type': 'text', 'required': True},
                 {'id': 'fastafile', 'type': 'txt'},
                 {'id': 'background', 'type': 'txt'},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'regions', 'type': 'userfile'},
                 {'id': 'motifs', 'type': 'text'},
                 {'id': 'customMotif', 'type': 'txt'},
                 {'id': 'threshold', 'type': 'float', 'required': True}]
out_parameters = [{'id': 'motif_results', 'type': 'track'}]

class MotifScanPlugin(OperationPlugin):
    info = {
        'title': 'Motif scanner',
        'description': 'Scan motifs PWM on a set of a sequences',
        'path': ['Signal', 'Motif scanner'],
        'output': MotifForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }


    def __call__(self, **kw):
        sequence_source = kw.get('sequenceSource')
        fasta_file = kw.get('fastafile')
        background = kw.get('background') or None
        assembly_id = kw.get('assembly') or None
        regions_file = kw.get('regions')
        motifs_list = kw.get('motifs')
        motif_add = kw.get('customMotif')
        threshold = float(kw.get('threshold') or 0)
 
        if motifs_list is None: motifs_list = []
        if not isinstance(motifs_list, list): motifs_list = [motifs_list]
 
        if background is None and assembly is None:
            background = self.temporary_path(fname='background_')
            with open(background,"w") as bgr: bgr.write('1 0.25 0.25 0.25 0.25')

        if assembly_id is not None:
            assembly = genrep.Assembly(assembly_id)
        else:
            if regions_file is not None:
                raise ValueError("Please specify an assembly if you specify regions.")
            assembly = None

        motifs = []
        if motif_add is not None:
            motifs.append({"name": "Custom Motif", "file": motif_add})
        for mot in motifs_list:
            gid, name = mot.split(' ')
            motif = g.get_genrep_objects('genomes/'+str(gid)+'/get_matrix', 'motif',
                                         params={"gene_name": name})
#### get file path directly?
            motifs.append({"name": motif.name, "file": self.saveMotifToFile(motif)})

        if len(motifs) == 0:
            raise ValueError("Please give at least one motif to scan for")

        with execution(None) as ex:
            track_output = save_motif_profile( ex, motifs, assembly, regions_file, fasta_file, 
                                               background=background, threshold=threshold, 
                                               description=None, via='local' )
        if track_output is not None:
            self.new_file(track_output, 'motif_track')
        return 1


