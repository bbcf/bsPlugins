from bsPlugins import *
from bbcflib.bFlatMajor.stream import merge_scores
from bbcflib.bFlatMajor.numeric import correlation
from bbcflib import btrack as track
from bbcflib import genrep


class MergeTracksForm(BaseForm):
    forward = twf.FileField(label_text='Forward: ',
        help_text='Select forward density file',
        validator=twf.FileValidator(required=True))
    reverse = twf.FileField(label_text='Reverse: ',
        help_text='Select reverse density file',
        validator=twf.FileValidator(required=True))
    assembly = twf.SingleSelectField(label_text='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    shift = twf.TextField(label_text='Shift: ',
        validator=twc.IntValidator(required=True),
        value=0,
        help_text='Enter positive downstream shift ([fragment_size-read_length]/2), or a negative value to estimate shift by cross-correlation')
    submit = twf.SubmitButton(id="submit", value='Merge tracks')


meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'forward', 'type': 'track', 'required': True},
                 {'id': 'reverse', 'type': 'track', 'required': True},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'shift', 'type': 'int', 'required': True}]
out_parameters = [{'id': 'density_merged', 'type': 'track'}]


class MergeTracksPlugin(OperationPlugin):

    info = {
        'title': 'Merge Tracks',
        'description': 'Shift and average scores from forward and reverse strand densities',
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

        tfwd = track.track(kw.get('forward'), chrmeta=kw.get('assembly') or None)
        trev = track.track(kw.get('reverse'), chrmeta=kw.get('assembly') or None)
        if not kw.get('assembly'):  # btrack does the job, take the max of both chromosome lengths
            chrmeta = tfwd.chrmeta
            for k, v in trev.chrmeta.iteritems():
                chrmeta.setdefault(k, {})['length'] = max(v['length'], chrmeta.get(k, {}).get('length', 0))
        elif tfwd.chrmeta:
            chrmeta = tfwd.chrmeta  # For sql files, btrack doesn't make it,
        elif trev.chrmeta:
            chrmeta = trev.chrmeta  # so one can contain the info while the second does not.
        else:
            raise ValueError("Must specify an assembly.")  # In case nothing works - should not happen

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
                    shiftval = (max_xcor_idx - xcor_lim - 1) / 2
                    #print "Autocorrelation shift=%i, correlation is %f at index %d for chromosome %s." \
                    #       % (shiftval,xcor[max_xcor_idx],max_xcor_idx,chrom)
                    break
            if not shiftval:
                raise ValueError("Unable to detect shift automatically. Must specify a shift value.")

        output = self.temporary_path(fname='density_merged', ext='sql')
        fields = ['chr', 'start', 'end', 'score']
        tout = track.track(output, format='sql', fields=fields, chrmeta=chrmeta,
                           info={'datatype': 'quantitative'})
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
        return 1
