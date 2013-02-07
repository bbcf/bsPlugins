from . import *
from bein import execution
from bbcflib import genrep
from bbcflib.motif import save_motif_profile
from bbcflib.btrack import track, FeatureStream
import re, os, urllib


input_types = [(0, 'Custom upload'), (1, 'From assembly')]
input_map = {input_types[0][1]: ['fastafile', 'background'],
             input_types[0][1]: ['assembly', 'regions']}

class MotifForm(DynForm):
    genomes = g.get_genrep_objects('genomes', 'genome') ####????
    motif_array = []
    for genome in genomes:
        if genome.motif_matrix_url != None and genome.motif_matrix_url != 'null':
            curMotifs = g.get_genrep_objects('genomes/' + str(genome.id) + '/get_matrix', 'motif') ####??
            for motif in curMotifs:
                motif_array.append((str(genome.id) + " " + motif.name, genome.name + " - " + motif.name))
    motif_array.sort(key=lambda motif: motif[1])

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
    motifs = twjqSelect.Select2MultipleSelectField(placeholder='Select the motifs to scan for', options=motif_array,
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
        try:
            sequence_source = kw.get('sequenceSource')
            fasta_file = kw.get('fastafile')
            background_file = kw.get('background')
            assembly_id = kw.get('assembly')
            regions_file = kw.get('regions')
            motifs_list = kw.get('motifs')

            if motifs_list == None:
                motifs_list = []
            if not isinstance(motifs_list, list):
                motifs_list = [motifs_list]

            motif_add = kw.get('customMotif',)
            threshold = int(kw.get('threshold') or 0)

            isAssembly = assembly_id is not None and regions_file is not None
            isSequence = fasta_file is not None

            if sequence_source == 'Custom input' or ((sequence_source == None or len(sequence_source) == 0) and isSequence and not isAssembly):
                sequence = seq.path
                background = background_file
                if background == None:
                    background = self.temp(fname='background', ext='mat',
                                           data='1 0.25 0.25 0.25 0.25')

                chrmeta = seq.chrmeta
            else:
                if not isAssembly:
                    raise ValueError('Invalid assembly supplied (check if an assembly was selected and that a region file was set)')
                assembly = genrep.Assembly(assembly_id)
                chrmeta = assembly.chrmeta

                stats = assembly.statistics(output='background.mat', frequency=True, matrix_format=True)
                sequence = self.temporary_path(fname="sequence", ext='fasta')
                (sequence, retrieved_length) = assembly.fasta_from_regions(str(regions_file), out=sequence)
                print("Got " + str(retrieved_length) + " nucleotides -> " + str(os.path.getsize(sequence)) + " bytes")

                if os.path.getsize(sequence) < retrieved_length or retrieved_length == 0:
                    self.dump_file(sequence)
                    raise ValueError("The given region wasn't retrieved correctly. Check if you gave the correct chromosome names (chr1 isn't the same as chrI!).")
            motifs = []
            if motif_add is not None:
                motifs.append({"name": "Custom Motif", "file": motif_add})

            for motif_entry in motifs_list:
                genome_id, motif_name = motif_entry.split(' ')
                [motif] = g.get_genrep_objects('genomes/' + str(genome_id) + '/get_matrix', 'motif', params={"gene_name": urllib.quote(str(motif_name))})
#### get file path directly?
                motifs.append({"name": motif.name, "file": self.saveMotifToFile(motif)})

            if len(motifs) == 0:
                raise ValueError("Please give at least one motif to scan for")

            with execution(None) as ex:
##### assembly?
                track_output = save_motif_profile( ex, motifs, assembly, regions, background=background,
                                                   threshold=threshold, description=None, via='local' )
            if track_output != None:
                self.new_file(track_output, 'motif_track')
            return 1
        except ValueError, e:
            return 'Error (' + self.display_time() + ').\n' + str(e)


