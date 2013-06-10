from bsPlugins import *
from bbcflib.bFlatMajor.stream import merge_scores
from bbcflib.bFlatMajor.numeric import correlation
from bbcflib import track as track
from bbcflib import genrep


class MergeTracksForm(BaseForm):
    forward = twb.BsFileField(label='Forward: ',
                              help_text='Select forward density file',
                              validator=twb.BsFileFieldValidator(required=True))
    reverse = twb.BsFileField(label='Reverse: ',
                              help_text='Select reverse density file',
                              validator=twb.BsFileFieldValidator(required=True))
    assembly = twf.SingleSelectField(label='Assembly: ',
                                     options=genrep.GenRep().assemblies_available(),
                                     help_text='Reference genome')
    shift = twf.TextField(label='Shift: ',
                          validator=twc.IntValidator(required=True),
                          value=0,
                          help_text='Enter positive downstream shift ([fragment_size-read_length]/2), \
                                     or a negative value to estimate shift by cross-correlation')
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bed",'bedGraph','wig','sga'],
        prompt_text=None,
        help_text='Format of the output file', )
    submit = twf.SubmitButton(id="submit", value='Merge tracks')


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'forward', 'type': 'track', 'required': True},
                 {'id': 'reverse', 'type': 'track', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'shift', 'type': 'int', 'required': True}]
out_parameters = [{'id': 'density_merged', 'type': 'track'}]


class MergeTracksPlugin(BasePlugin):
    """Shift and average scores from forward and reverse strand densities. <br /><br />
Typically built to merge ChIP-seq signals from both DNA strands, it can also be used to add (average)
several numeric genomic tracks, replicates for instance.<br />
The output is the average of all the input signals, position by position.
    """
    info = {
        'title': 'Merge strands',
        'description': __doc__,
        'path': ['Signal', 'Merge strands'],
        'output': MergeTracksForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }

    def __call__(self, **kw):
        def _shift(stream, shift):
            istart = stream.fields.index('start')
            iend = stream.fields.index('end')
            i1 = min(istart, iend)
            i2 = max(istart, iend)

            def _apply_shift(x):
                return x[:i1] + (x[i1] + shift,) + x[i1 + 1:i2] + (x[i2] + shift,) + x[i2 + 1:]
            return track.FeatureStream((_apply_shift(x) for x in stream),
                                       fields=stream.fields)

        assembly = kw.get('assembly') or 'guess'
        tfwd = track.track(kw.get('forward'), chrmeta=assembly)
        trev = track.track(kw.get('reverse'), chrmeta=assembly)
        chrmeta = tfwd.chrmeta

        shiftval = int(kw.get('shift', 0))
        if shiftval < 0:  # Determine shift automatically
            shiftval = None
            xcor_lim = 300
            for chrom, v in chrmeta.iteritems():
                chrsize = v['length']
                xcor_lim = min(xcor_lim, 0.01 * chrsize)
                xcor = correlation([tfwd.read(chrom), trev.read(chrom)], regions=(1, chrsize),
                                   limits=(-xcor_lim, xcor_lim))
                max_xcor_idx = xcor.argmax()
                if xcor[max_xcor_idx] > 0.2:
                    shiftval = (max_xcor_idx - xcor_lim - 1)/2
                    break
            if not shiftval:
                raise ValueError("Unable to detect shift automatically. Must specify a shift value.")

        output = self.temporary_path(fname=tfwd.name+'-'+trev.name+'_merged', ext=kw['format'])
        tout = track.track(output, chrmeta=chrmeta,
                           info={'datatype': 'quantitative', 'shift': shiftval})
        mode = 'write'
        for chrom in chrmeta.keys():
            tout.write(merge_scores([_shift(tfwd.read(selection=chrom), shiftval),
                                     _shift(trev.read(selection=chrom), -shiftval)]),
                       chrom=chrom, mode=mode, clip=True)
            mode = 'append'
        tout.close()
        trev.close()
        tfwd.close()
        self.new_file(output, 'density_merged')
        return self.display_time()

# nosetests --logging-filter=-tw2 test_MergeTracks.py
