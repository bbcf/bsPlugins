from bsPlugins import *
from bein import execution
from bbcflib.track import track, convert
from bbcflib.mapseq import bam_to_density
from bbcflib.gfminer.stream import merge_scores
import os, sys, pysam

__requires__ = ["pysam"]
output_opts=["sql", "bedGraph", "bigWig"]


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'sample', 'type': 'bam', 'required': True, 'multiple': True, 'label': ' Test BAMs: ', 'help_text': 'Select main bam file(s)'},
                 {'id': 'control', 'type': 'bam', 'label': 'Control BAM: ', 'help_text': 'Select control bam file to compute enrichment' },
                 {'id': 'output', 'type': 'listing', 'label': 'Output format: ', 'help_text': 'Format of the output file', 'options': output_opts, 'prompt_text': None},
                 {'id': 'normalization', 'type': 'int', 'label': 'Normalization: ', 'help_text': 'Normalization factor, default is total number of reads'},
                 {'id': 'merge_strands', 'type': 'int', 'label': 'Shift and merge strands: ', 'help_text': 'Shift value (in bp) if you want to merge strand-specific densities (will not merge if negative)', 'value': -1},
                 {'id': 'read_extension', 'type': 'int', 'label': 'Read extension: ','help_text': 'Read extension (in bp) to be applied when constructing densities (will use read length if negative)', 'value': -1 },
                 {'id': 'no_nh_flag', 'type':'boolean', 'required':True, 'label': 'Do not use NH flag: ', 'help_text': 'Do not use NH (multiple mapping counts) as weights', 'value': False},
                 {'id': 'single_end', 'type':'boolean', 'required':True, 'label': 'As single end: ', 'help_text': 'Considered a paired-end bam as single-end (default: False, namely whole-fragment densities instead of read densities)', 'value': False},
                 {'id': 'stranded', 'type':'boolean', 'required':True, 'label': 'As strand-specific: ', 'help_text': 'If the sequencing protocol was paired-end strand-specific, generate plus and minus densities (default: False)', 'value': False}
]
out_parameters = [{'id': 'density_merged', 'type': 'track'},
                  {'id': 'density_fwd', 'type': 'track'},
                  {'id': 'density_rev', 'type': 'track'},
                  {'id': 'density_plus_merged', 'type': 'track'},
                  {'id': 'density_plus_fwd', 'type': 'track'},
                  {'id': 'density_plus_rev', 'type': 'track'},
                  {'id': 'density_minus_merged', 'type': 'track'},
                  {'id': 'density_minus_fwd', 'type': 'track'},
                  {'id': 'density_minus_rev', 'type': 'track'}]


class Bam2DensityForm(BaseForm):
    class BamMulti(twb.BsMultiple):
        label='Test BAMs: '
        sample = twb.BsFileField(label=' ',
                                 help_text='Select main bam file(s)',
                                 validator=twb.BsFileFieldValidator(required=True))
    control = twb.BsFileField(label='Control BAM: ',
                              help_text='Select control bam file to compute enrichment',
                              validator=twb.BsFileFieldValidator(required=False))
    format = twf.SingleSelectField(label='Output format: ',
                                   options=["sql", "bedGraph", "bigWig"],
                                   prompt_text=None,
                                   help_text='Format of the output file')
    normalization = twf.TextField(label='Normalization: ',
                                  validator=twc.IntValidator(),
                                  help_text='Normalization factor, default is total number of reads')
    merge_strands = twf.TextField(label='Shift and merge strands: ',
                                  validator=twc.IntValidator(),
                                  value=-1,
                                  help_text='Shift value (in bp) if you want to merge strand-specific densities (will not merge if negative)')
    read_extension = twf.TextField(label='Read extension: ',
                                   validator=twc.IntValidator(),
                                   value=-1,
                                   help_text='Read extension (in bp) to be applied when constructing densities (will use read length if negative)')
    single_end = twf.CheckBox(label='As single end: ',
                              value=False,
                              help_text='Considered a paired-end bam as single-end (default: False, namely whole-fragment densities instead of read densities)')
    no_nh_flag = twf.CheckBox(label='Do not use NH flag: ',
                              value=False,
                              help_text='Do not use NH (multiple mapping counts) as weights')
    stranded = twf.CheckBox(label='As strand-specific: ',
                              value=False,
                              help_text='If the sequencing protocol was paired-end strand-specific, generate plus and minus densities (default: False)')
    submit = twf.SubmitButton(id="submit", value='bam2density')


class Bam2DensityPlugin(BasePlugin):
    """From a BAM file, creates a track file of the read count/density along the whole genome,
in the chosen format.

Read counts are divided by 10^-7 times the normalization factor (which is total number of reads by default).
Positive and negative strand densities are generated and optionally merged (averaged) if a
shift value >=0 is given. The read extension is the number of basepairs a read will cover,
starting from its most 5' position (e.g. with a read extension of 1, only the starting position of
each alignment will be considered, default is read length).
"""
    info = {
        'title': 'Genome-wide reads density from BAM',
        'description': __doc__,
        'path': ['Files', 'Bam2density'],
#        'output': Bam2DensityForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        b2wargs = []
        control = None
        #samples = kw.get('BamMulti',{}).get('sample', [])
        samples = kw.get('sample', [])
        if not isinstance(samples, list): samples = [samples]
        samples = {'_': [os.path.abspath(s) for s in samples if os.path.exists(s)]}
        if kw.get('control'):
            control = kw['control']
            b2wargs = ["-c", str(control)]
            assert os.path.exists(str(control)), "Control file not found: '%s'." % control
            control = os.path.abspath(control)
        try:
            nreads = int(kw.get('normalization'))
        except (ValueError, TypeError):
            nreads = -1
        bamfiles = [track(s, format='bam') for s in samples]
        if nreads < 0:
            _nreads = [0]*len(samples)
            if control is not None:
                b2wargs += ["-r"]
        else:
            _nreads = [nreads for s in samples]
        try:
            merge_strands = int(kw.get('merge_strands'))
        except (ValueError, TypeError):
            merge_strands = -1
        try:
            read_extension = int(kw.get('read_extension'))
        except (ValueError, TypeError):
            read_extension = -1
        single_end = kw.get('single_end',False)
        if isinstance(single_end, basestring):
            single_end = (single_end.lower() in ['1', 'true', 't','on'])
        no_nh = kw.get('no_nh_flag',False)
        if isinstance(no_nh, basestring):
            no_nh = (no_nh.lower() in ['1', 'true', 't','on'])
        if no_nh:  b2wargs += ["--no_nh"]
        output = {'_': [self.temporary_path(fname=b.name+'_density_') for b in bamfiles]}
        stranded = kw.get('stranded',False)
        if isinstance(stranded, basestring):
            stranded = (stranded.lower() in ['1', 'true', 't','on'])
        if stranded:
            if single_end: sys.exit("Error: the option stranded only works with paired-end data")
            output = {'_plus_': [], '_minus_': []}
            samples = {'_plus_': [], '_minus_': []}
            trout = {}
            for bam in bamfiles:
                for orient in output:
                    bamname = self.temporary_path(fname= bam.name+orient+".bam")
                    trout[orient] = pysam.Samfile(bamname, "wb", template=bam.filehandle)
                    samples[orient].append(os.path.abspath(bamname))
                    outname = self.temporary_path(fname=bam.name+orient)
                    output[orient].append(os.path.abspath(outname))
                bam.open()
                for read in bam.filehandle:
                    if not (read.is_paired and read.is_proper_pair): continue
                    if (read.is_read1 and read.is_reverse) or (read.is_read2 and read.mate_is_reverse):
                        trout['_plus_'].write(read)
                    elif (read.is_read2 and read.is_reverse) or (read.is_read1 and read.mate_is_reverse):
                        trout['_minus_'].write(read)
                (t.close() for t in trout.values())
        format = kw.get('output', 'sql')
        info = {'datatype': 'quantitative', 'read_extension': read_extension}
        if merge_strands >= 0:
            suffixes = ["merged"]
            info['shift'] = merge_strands
        else:
            suffixes = ["fwd", "rev"]
        chrmeta = bamfiles[0].chrmeta
        files = {}
        with execution(None) as ex:
            files = dict((o,[bam_to_density( ex, s, output[o][n], nreads=_nreads[n],
                                             merge=merge_strands,
                                             read_extension=read_extension,
                                             sql=True, se=single_end, args=b2wargs )
                             for n,s in enumerate(s)]) for o,s in samples.items())
        for suf in suffixes:
            all_s_files = dict((o,[x for y in files for x in y if x.endswith(suf+".sql")]) for o,f in files.items())
            for orient, sfiles in all_s_files.iteritems():
                if len(sfiles) > 1:
                    x = self.temporary_path(fname="Density_average"+orient+suf+".sql")
                    tsql = track( x, fields=['start', 'end', 'score'], chrmeta=chrmeta, info=info )
                    insql = []
                    for f in sfiles:
                        t = track( f, format='sql', fields=['start', 'end', 'score'],
                                   chrmeta=chrmeta, info=info )
                        t.save()
                        insql.append(t)
                    for c in tsql.chrmeta:
                        tsql.write(merge_scores([t.read(c) for t in insql]),chrom=c)
                else:
                    x = sfiles[0]
                    tsql = track( x, format='sql', fields=['start', 'end', 'score'],
                                  chrmeta=chrmeta, info=info )
                    tsql.save()
                if format in [None,"sql"]:
                    outname = x
                else:
                    outname = os.path.splitext(x)[0]+"."+format
                    convert(x, outname, mode="overwrite")
                self.new_file(outname, 'density'+orient+suf)
        return self.display_time()
