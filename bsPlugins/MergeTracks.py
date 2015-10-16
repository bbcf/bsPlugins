from bsPlugins import *
from bbcflib.gfminer.stream import merge_scores
from bbcflib.gfminer.numeric import correlation
from bbcflib.track import track, FeatureStream
from bbcflib import genrep

output_opts = ['sql','bed','bedGraph','wig','bigWig','sga']

method_opts = ['mean','min','max','geometric','median','sum']

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'forward', 'type': 'track', 'required': True, 'label': 'Forward: ', 'help_text': 'Select forward density file' },
                 {'id': 'reverse', 'type': 'track', 'required': True, 'label': 'Reverse: ', 'help_text': 'Select reverse density file' },
                 {'id': 'assembly', 'type': 'assembly', 'label': 'Assembly: ', 'help_text': 'Reference genome', 'options': genrep.GenRep().assemblies_available()},
                 {'id': 'shift', 'type': 'int', 'required': True, 'label': 'Shift: ', 'help_text': 'Enter positive downstream shift ([fragment_size-read_length]/2), \nor a negative value to estimate shift by cross-correlation', 'value': 0},
                 {'id': 'output', 'type': 'listing', 'label': 'Output format: ', 'help_text': 'Format of the output file', 'options': output_opts, 'prompt_text': None},
                 {'id': 'method', 'type': 'radio', 'label': 'Method: ', 'help_text': 'Select the score combination method', 'options': ['mean','min','max','geometric','median','sum'], 'value': 'mean'}]
out_parameters = [{'id': 'density_merged', 'type': 'track'}]



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
                                   options=['sql','bed','bedGraph','wig','bigWig','sga'],
                                   prompt_text=None,
                                   help_text='Format of the output file', )
    method = twf.RadioButtonList(label='Method: ',
                                 options=['mean','min','max','geometric','median','sum'],
                                 value='mean',
                                 help_text='Select the score combination method')
    submit = twf.SubmitButton(id="submit", value='Merge tracks')


class MergeTracksPlugin(BasePlugin):
    """Shift and average scores from forward and reverse strand densities.

Typically built to merge ChIP-seq signals from both DNA strands, it can also be used to add (average)
several numeric genomic tracks, replicates for instance.

The output is the average of all the input signals, position by position.
    """
    info = {
        'title': 'Merge strands',
        'description': __doc__,
        'path': ['Signal', 'Merge strands'],
#        'output': MergeTracksForm,
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
            return FeatureStream((_apply_shift(x) for x in stream),
                                       fields=stream.fields)

        assembly = kw.get('assembly') or 'guess'
        tfwd = track(kw.get('forward'), chrmeta=assembly)
        trev = track(kw.get('reverse'), chrmeta=assembly)
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

        output = self.temporary_path(fname=tfwd.name+'-'+trev.name+'_merged', 
                                     ext=kw.get('format',tfwd.format))
        tout = track(output, chrmeta=chrmeta,
                     info={'datatype': 'quantitative', 'shift': shiftval})
        mode = 'write'
        method = kw.get("method","mean")
        for chrom in chrmeta.keys():
            tout.write(merge_scores([_shift(tfwd.read(selection=chrom), shiftval),
                                     _shift(trev.read(selection=chrom), -shiftval)],
                                    method=method),
                       chrom=chrom, mode=mode, clip=True)
            mode = 'append'
        tout.close()
        trev.close()
        tfwd.close()
        self.new_file(output, 'density_merged')
        return self.display_time()

# nosetests --logging-filter=-tw2 test_MergeTracks.py
