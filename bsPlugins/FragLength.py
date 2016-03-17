from bsPlugins import *
import tarfile, os, sys, pysam

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'bamfiles', 'type': 'bam', 'required': True, 'multiple': 'BamMulti', 'label': 'Paired-ended BAM files: ' },
                 {'id': 'minlength', 'type': 'int', 'label': 'Minimum fragment length: '},
                 {'id': 'maxlength', 'type': 'int', 'label': 'Maximum fragment length: '}]
out_parameters = [{'id': 'fragment_track', 'type': 'track'},
                  {'id': 'fragment_track_tar', 'type': 'file'}]

class FragLengthForm(BaseForm):
    class BamMulti(twb.BsMultiple):
        label = 'Paired-end BAM files: '
        bamfiles = twb.BsFileField(label=' ',
                                 validator=twb.BsFileFieldValidator(required=True))
    minlength = twf.TextField(label='Minimum fragment length: ',
                              validator=twc.IntValidator(required=False))
    maxlength = twf.TextField(label='Maximum fragment length: ',
                              validator=twc.IntValidator(required=False))
    submit = twf.SubmitButton(id="submit", value="Filter")

class FragLengthPlugin(BasePlugin):
    """Computes BAM files with fragment lengths between minlength and maxlength."""
    info = {
        'title': 'Filtering of paired-end BAM files by fragment length',
        'description': __doc__,
        'path': ['Files', 'Fragment length filtering'],
#        'output': FragLengthForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta}

    def __call__(self, **kw):
        #bamfiles = kw.get('BamMulti',{}).get('bamfiles',[])
        bamfiles = kw.get('bamfiles',[])
        if not isinstance(bamfiles, (tuple,list)): bamfiles = [bamfiles]
        bamfiles = [pysam.Samfile(bam) for bam in bamfiles]
        minlength = kw.get('minlength')
        maxlength = kw.get('maxlength')
        if minlength or maxlength:
            minlength = int(minlength or 0)
            maxlength = int(maxlength or sys.maxint)
            if minlength > maxlength:
                raise ValueError("Empty range: %i:%i" %(minlength,maxlength))
        all_tracks = []
        for bam in bamfiles:
            tname = bam.filename.split(".")[0]+"_minlength"+str(minlength)+"_maxlength"+str(maxlength)+".bam"
            outname = self.temporary_path(fname=tname)
            all_tracks.append(outname)
            trout = pysam.Samfile(outname, "wb", template=bam)
            for read in bam:
                if not read.is_proper_pair: continue
                if read.is_reverse:
                    size = -read.isize
                else:
                    size = read.isize
                if (size >= minlength and size <= maxlength):
                    trout.write(read)
        if len(all_tracks) > 1:
            tarname = self.temporary_path(fname='BAM_filtered_by_fragment_length.tgz')
            tar_tracks = tarfile.open(tarname, "w:gz")
            [tar_tracks.add(f,arcname=os.path.basename(f)) for f in all_tracks]
            tar_tracks.close()
            self.new_file(tarname, 'fragment_track_tar')
        else:
            self.new_file(all_tracks[0], 'fragment_track')
        return self.display_time()
