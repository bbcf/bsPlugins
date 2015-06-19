from bsPlugins import *
from bein import execution
from bbcflib.common import fasta_length
from bbcflib.motif import meme
from bbcflib import genrep
from bbcflib.track import track, FeatureStream
import os, tarfile

input_types = [(0, 'Fasta upload'), (1, 'Select regions from genome')]
input_map = {0: ['fastafile'], 1: ['assembly', 'regions']}
_nm = 4

assembly_list = genrep.GenRep().assemblies_available()

meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}

in_parameters = [{'id': 'input_type', 'type': 'radio'},
                 {'id': 'fastafile', 'type': 'userfile'},
                 {'id': 'assembly', 'type': 'assembly'},
                 {'id': 'regions', 'type': 'track'},
                 {'id': 'nmotifs', 'type': 'int'}]
out_parameters = [{'id': 'meme_archive', 'type': 'file'}]


class MotifSearchForm(BaseForm):

    child = twd.HidingTableLayout()
    input_type = twd.HidingRadioButtonList(label='Sequence source: ',
                                           options=input_types,
                                           value=0,
                                           mapping=input_map)

    fastafile = twb.BsFileField(label='Fasta file: ', help_text='Sequences to scan')
    assembly = twf.SingleSelectField(label='Assembly: ', options=assembly_list,
                                     prompt_text=None,
                                     help_text='Assembly to fetch sequences from')
    regions = twb.BsFileField(label='Regions: ', help_text='Genomic regions to search (e.g. bed)',
                              validator=twb.BsFileFieldValidator())
    nmotifs = twf.TextField(label='Number of motifs: ',
                            validator=twc.IntValidator(),
                            value=_nm,
                            help_text='Number of motifs to search')
    submit = twf.SubmitButton(id="submit", value='Search motifs')


class MotifSearchPlugin(BasePlugin):
    """Search over-represented motifs in a set of a genomic regions using MEME"""
    info = {
        'title': 'Motif search by MEME',
        'description': __doc__,
        'path': ['Sequence analysis', 'Motif search'],
        'output': MotifSearchForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }


    def __call__(self, **kw):
        input_type = kw.get('input_type', 0)
        ass = kw.get('assembly','')
        regions_file = kw.get('regions') or ''
        if regions_file: regions_file = os.path.abspath(regions_file)
        fasta = kw.get('fastafile') or ''
        if fasta: fasta = os.path.abspath(fasta)
        else: fasta = os.path.abspath(self.temporary_path(fname=regions.name+'.fa'))
        output = self.temporary_path(fname=name+"_meme.tgz")
        outdir = os.path.abspath(os.path.join(os.path.split(fasta)[0],name+"_meme"))
        bfile = self.temporary_path(fname="background")
        with execution(None) as ex:
            if str(input_type) in [str(x[0]) for x in input_types]:
                input_type = int(input_type)
            if input_type in input_types[0]: #fasta
                name = os.path.splitext(os.path.basename(fasta))[0]
                if ass in [x[0] for x in genrep.GenRep().assemblies_available()]:
                    assembly = genrep.Assembly(ass)
                else:
                    assembly = genrep.Assembly(ex=ex,fasta=fasta)
                size = None
            elif input_type in input_types[1]: #regions
                assembly = genrep.Assembly(ass)
                if not os.path.exists(regions_file):
                    raise ValueError("File not found: %s" %regions_file)
                regions = track(regions_file,chrmeta=assembly.chrmeta)
                name = regions.name
                gRef = assembly.fasta_by_chrom
                (fasta, size) = assembly.fasta_from_regions( list(regions.read(fields=['chr','start','end'])), 
                                                             out=fasta, path_to_ref=gRef )
            else:
                raise ValueError("Input type not implemented: %s" %input_type)
            background = assembly.statistics(bfile, frequency=True)
            meme_args = kw.get("meme_args",[])
            nmotifs = kw.get('nmotifs') or _nm
            if not '-nmotifs' in meme_args:
                meme_args += ['-nmotifs',"%i" %int(nmotifs)]
            if size is None: size = sum(fasta_length(ex,fasta).values())
            meme_out = meme( ex, fasta, outdir, background, maxsize=(size*3)/2, 
                             args=meme_args )
        tarf = tarfile.open(output, "w:gz")
        tarf.add(outdir,arcname=os.path.basename(outdir))
        tarf.add(fasta,arcname=os.path.basename(fasta))
        tarf.close()
        self.new_file(output, 'meme_archive')
        return self.display_time()



