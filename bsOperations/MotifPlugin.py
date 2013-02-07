from . import *
from bein import execution
from bbcflib import genrep, motif
from bbcflib.motif import motif_scan
from bbcflib.btrack import track, FeatureStream
import re, os, urllib


class MotifForm(DynForm):
    #Build the motif array
    genomes = g.get_genrep_objects('genomes', 'genome')
    motif_genomes = []
    motif_array = []

    #jQuery's already on the page
    twj.jquery_js.no_inject = True

    for genome in genomes:
        if genome.motif_matrix_url != None and genome.motif_matrix_url != 'null':
            motif_genomes.append(genome)
            curMotifs = g.get_genrep_objects('genomes/' + str(genome.id) + '/get_matrix', 'motif')
            for motif in curMotifs:
                motif_array.append((str(genome.id) + " " + motif.name, genome.name + " - " + motif.name))

    motif_array.sort(key=lambda motif: motif[1])

    #Build form
    sequenceSource = twd.HidingRadioButtonList(
        label='Sequence source', options=('Custom input', 'From assembly'),
        mapping={
            'Custom input': ['sequence', 'background'],
            'From assembly': ['assembly', 'regions']
        })

    spacer12 = twf.Spacer()

    fastafile = twf.FileField(label='FASTA: ', help_text='FASTA file')
    background = twf.FileField(label='Background: ',
        help_text='A matrix file containing the background frequencies (default: uniform background)')

    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Assembly used as data source')
    regions = twf.FileField(label='Regions: ',
        help_text='Regions to scan (e.g. bed)')

    spacer23 = twf.Spacer()

    motifs = twjqSelect.Select2MultipleSelectField(
        placeholder='Select the motifs to scan for',
        options=motif_array,
        help_text='')

    customMotif = twf.FileField(label='Custom motif: ',
        help_text='An optional custom/additional motif to scan (.mat)')

    threshold = twf.TextField(label='Threshold: ',
        value='0.0',
        validator=twc.RegexValidator(required=True, regex=re.compile('^\d+(\.\d+)?$')))

    submit = twf.SubmitButton(id="submit", value='Scan for Motifs')


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


class MotifPlugin(OperationPlugin):
    info = {
        'title': 'Motif scanner',
        'description': 'Scan motifs PWM on a set of a sequences',
        'path': ['Signal', 'Motif scanner'],
        'output': MotifForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def motif_scan(self, fasta, motif, background, threshold):
        """Launch the command-line utility that scans for the given motif.

        Fasta, motif and background should be paths, threshold a double/float"""
        #FIXME: This should use bein (the two following commented lines), but it wouldn't work on the
        #local machine. Change if possible.
#        with execution(None) as ex:
#            return motif_scan(ex, fasta, motif, background, threshold)
        try:
            #Assumes S1K is in the PATH. See last comment for a better solution
            return subprocess.check_output(["S1K", motif, background, threshold, fasta])
        except subprocess.CalledProcessError, e:
            raise e

    def dump_file(self, file):
        """Dump a file to the command-line (for debugging)"""
        f = open(file, 'r')
        for line in f:
            if len(line) > 80:
                print(str(line[0:77] + "..."))
            else:
                print(str(line))
        f.close()

    def motif_scan_to_track(self, fasta, motifName, motif, background, threshold, chrmeta, output=None):
        """Perform a motif scan and write the results to a track.
        It executes motif_scan(fasta, motif, background, threshold) and inserts all results back into a track. If
        the output track is None, a SQL track is created.
        Returns the track"""
        #The buffer size (used to speed up insertion into SQL tracks althoug it will probably help for most formats)
        COLLECT_SIZE = 1000
        results = self.motif_scan(fasta, motif, background, threshold)
        if output == None:
            output = self.temporary_path(fname='motif_finder_results', ext='sql')
        track_output = track(output, fields=['start', 'end', 'score', 'name', 'strand'], chrmeta=chrmeta, info={'datatype': 'features'})
        #Sample: "chr1|chr1:1-230207" -> ["chr1", "chr1", "1", "230207"]
        parse_name = re.compile("^(.*)\|(.*):(.*)-(.*)$")

        lines = results.splitlines()

        features = []
        for line in lines:
            # name: Name of the FASTA part
            # seq: Matched sequence
            # score: Score
            # pos: Starting position (1 -> first nucleotide)
            # strand: +/- -> Watson/Crick
            [name, seq, score, pos, strand] = line.split("\t")

            score = float(score)
            pos = int(pos) - 1
            length = len(seq)
            regionFrom = 0
            regionTo = length

            #Name parsing is a bit more complicated as we need handle more different cases. If the name is in the assembly
            #format, it can be parsed by the parse_name regex. If not, the name is taken as-is and the positions
            #(regionFrom & regionTo) are assumed to be simple (0 -> length).
            fullName = motifName
            #Parse name
            if parse_name.match(name) != None:
                #Sample: ">chr1|chr1:1-230207"
                [(name, _, regionFrom, regionTo)] = parse_name.findall(name)
                regionFrom = int(regionFrom)
                regionTo = int(regionTo)

            #Generate a more explicit name
            if chrmeta[name] != None:
                if hasattr(chrmeta[name], 'real_name') and chrmeta[name]['real_name'] != None:
                    fullName = chrmeta[name]['real_name'] + " - " + motifName

            #Most track formats doesn't handle the case where to < from -> flip to correct
            if regionTo < regionFrom:
                strand = "+" if strand.strip() == "-" else "+"
                [regionFrom, regionTo] = [regionTo, regionFrom]

            features.append((name, regionFrom + pos, regionFrom + pos + length, score, fullName, strand))

            if len(features) >= COLLECT_SIZE:
                #Buffer full -> flush
                stream = FeatureStream(features, fields=['chr', 'start', 'end', 'score', 'name', 'strand'])
                track_output.write(stream)
                features = []

        if len(features) > 0:
            #Finished -> flush
            stream = FeatureStream(features, fields=['chr', 'start', 'end', 'score', 'name', 'strand'])
            track_output.write(stream)
        track_output.close()

        return output

    def temp(self, data, fname='temp', ext='txt'):
        """Create a temporary file and write data into it"""
        file = self.temporary_path(fname=fname, ext=ext)
        f = open(file, 'w')
        f.write(data)
        f.close()
        return file

    def parseMotifFile(self, motifFile):
        """Parse a motif matrix file and return an equivalent GenrepObject (the file you\'d receive if the motif existed on the server)."""
        result = {"name": "Custom",
                  "n": 0,
                  "alphabet": ["A", "C", "G", "T"],
                  "motif": []}
        f = open(motifFile, 'r')

        for line in f:
            line = line.strip()
            if len(line) > 0 and not line.startswith(">"):
                [one, a, c, g, t] = re.findall(r"[^\s]+", line)
                result["motif"].append([a, c, g, t])
        f.close()
        result["n"] = len(result["motif"])
        return GenrepObject({"motif": result}, "motif")

    def saveMotifToFile(self, motif, output_file=None):
        """Save a GenrepObject motif to a file. Returns the path to the written file."""
        if output_file == None:
            output_file = self.temporary_path(fname='motif', ext='mat')
        f = open(output_file, 'w')
        #Transposition table
        t = {}
        for i in range(len(motif.alphabet)):
            t[motif.alphabet[i]] = i

        first = True
        for row in motif.motif:
            if first:
                first = False
            else:
                f.write("\n")
            f.write("1 " + str(row[t['A']]) + " " + str(row[t['C']]) + " " + str(row[t['G']]) + " " + str(row[t['T']]))

        f.close()

        return output_file

    def __call__(self, **kw):
        try:
            #Get the input
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

            motif_add = kw.get('customMotif')
            threshold = kw.get('threshold')

            isAssembly = not (assembly_id == None or len(assembly_id) == 0 or regions_file == None or len(regions_file) == 0)
            isSequence = not (fasta_file == None or len(fasta_file) == 0)

            if sequence_source == 'Custom input' or ((sequence_source == None or len(sequence_source) == 0) and isSequence and not isAssembly):
                #The user is supplying his own sequence
                if not isSequence:
                    raise ValueError('Invalid sequence supplied (check if a sequence file was set)')

                #Rewrite to the expected format (and simultaneously check whether it's valid)
                seq = CustomSequenceFile(self, fasta_file)
                if not seq.validate():
                    raise ValueError("Invalid sequence file given")

                #Get the corrected sequence file
                sequence = seq.path

                background = background_file

                #No background -> uniform
                if background == None:
                    background = self.temp(fname='background', ext='mat',
                                           data='1 0.25 0.25 0.25 0.25')

                chrmeta = seq.chrmeta
            else:
                #The user is using an assembly sequence
                if not isAssembly:
                    raise ValueError('Invalid assembly supplied (check if an assembly was selected and that a region file was set)')
                assembly = genrep.Assembly(assembly_id)
                chrmeta = assembly.chrmeta

                #The statistics method doesn't work completely as advertised, construct background file manually
                stats = assembly.statistics()
                tot = float(stats['A'] + stats['C'] + stats['G'] + stats['T'])
                background = self.temp(fname='background', ext='mat',
                                       data='1 ' + str(stats['A'] / tot) + " " + str(stats['C'] / tot) + " " + str(stats['G'] / tot) + " " + str(stats['T'] / tot))

                if regions_file.endswith(".bed"):
                    #BED file handling is a bit less flexible than it should be. Replace spaces by tabs:
                    f = open(str(regions_file), 'r')
                    regions_file = self.temporary_path(fname="region", ext='bed')
                    fOut = open(regions_file, 'w')
                    for line in f:
                        fOut.write("\t".join(re.findall("[^\s]+", line)) + "\n")
                    f.close()
                    fOut.close()
                #Get the asked sequence
                sequence = self.temporary_path(fname="sequence", ext='fasta')
                (sequence, retrieved_length) = assembly.fasta_from_regions(str(regions_file), out=sequence)
                print("Got " + str(retrieved_length) + " nucleotides -> " + str(os.path.getsize(sequence)) + " bytes")

                if os.path.getsize(sequence) < retrieved_length or retrieved_length == 0:
                    self.dump_file(sequence)
                    raise ValueError("The given region wasn't retrieved correctly. Check if you gave the correct chromosome names (chr1 isn't the same as chrI!).")

                #Pass the file through the CustomSequenceFile class -> prevents certain bugs in S1K
                sequence = CustomSequenceFile(self, sequence).path

            motifs = []

            if motif_add != None and len(motif_add) > 0:
                #Make sure the motif file has the correct format (this parser is
                #a bit more liberal than the one in S1K, so it'll correct some
                #of the mistakes
                motif = self.parseMotifFile(motif_add)
                motif_add = self.saveMotifToFile(motif)
                motifs.append({"name": "Custom Motif", "file": motif_add})

            g = genrep.GenRep('127.0.0.1:3000', '')

            for motif_entry in motifs_list:
                [genome_id, motif_name] = re.match("^([0-9]*) (.*)$", motif_entry).groups()
                [motif] = g.get_genrep_objects('genomes/' + str(genome_id) + '/get_matrix', 'motif', params={"gene_name": urllib.quote(str(motif_name))})

                motifs.append({
                    "name": motif.name,
                    "file": self.saveMotifToFile(motif)
                })

            if len(motifs) == 0:
                raise ValueError("Please give at least one motif to scan for")

            track_output = None  # self.temporary_path(fname='motif_scan_track', ext='sql')

            for motif in motifs:
                track_output = self.motif_scan_to_track(sequence, motif["name"], motif["file"], background, threshold, chrmeta, track_output)

            if track_output != None:
                self.new_file(track_output, 'motif_results')
            return 'Success (' + self.display_time() + ")."
        except ValueError, e:
            #ValueError generated locally -> Show a nicer error message
            return 'Error (' + self.display_time() + ').\n' + str(e)


class CustomSequenceFile:
    """Class used to handle specific problems with custom sequence files and S1K.
    When creating the class, a corrected version of the sequence is immediately
    created (the path is saved in the "path" property)."""
    def __init__(self, motif_plugin, path, chr=None, regionFrom=None, regionTo=None):
        self.motif_plugin = motif_plugin
        self.path = path
        self.chr = chr if chr != None else 'chr1'
        self.elements = []

        self.transformFile()

    def validate(self):
        "Validate the sequence. Doesn't do anything."
        #TODO: Implement ;D
        return True

    def getName(self, chr):
        """Get the real name for a certain chromosome.

        chr: The chromosome identifier (i.e. pass 1 to get the name of chr1)"""
        if isinstance(chr, int):
            chr = 'chr' + str(chr)

        for el in self.elements:
            if chr == 'chr' + str(el['chr']):
                return el["name"]

        return None

    @property
    def chrmeta(self):
        """Get the chrmeta of the current sequence"""
        meta = {}

        for el in self.elements:
            meta['chr' + str(el['chr'])] = {"length": el["length"], "real_name": el['name']}

        return meta

    def transformFile(self):
        """Private. Transform the file so S1K will accept it"""
        MAX_LINE_SIZE = 1000000

        f = open(self.path, 'r')
        tmpPath = self.motif_plugin.temporary_path(fname='sequence', ext='fasta')
        fOut = open(tmpPath, 'w')
        i = 1
        copiedLine = False
        hasName = False
        for line in f:
            line = str(line)
            line = line.strip()

            if len(line) > 0:  # Remove all empty lines
                if line.startswith(">"):
                    hasName = True
                    fOut.write('>chr' + str(i) + "\n")
                    self.elements.append({"name": line[1:].strip(), "chr": i, "length": 0})
                    i += 1
                else:
                    if not hasName:
                        #EVERY part needs a name, if there's none in the file -> create one!
                        fOut.write('>chr' + str(i) + "\n")
                        self.elements.append({"name": '', "chr": i, "length": 0})
                        i += 1
                    self.elements[-1]['length'] += len(line)

                    #S1K can't handle lines of arbitrary length -> cut them into blocks
                    for j in xrange(0, len(line), MAX_LINE_SIZE):
                        fOut.write(line[j:j + MAX_LINE_SIZE] + "\n")

        f.close()
        fOut.close()

        self.path = tmpPath
