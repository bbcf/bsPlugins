from bsPlugins import *
from bbcflib.bFlatMajor import stream as gm_stream
from bbcflib import btrack as track
from bbcflib import genrep
from bbcflib.btrack import track
from bbcflib.bFlatMajor import common
from bbcflib.bFlatMajor.common import split_field,apply
import re
import math #log2, log10, sqrt, ...
__requires__ = ["numpy"]
# import toscawidget2 modules in order to build forms
import tw2.forms as twf
class loga2Form(BaseForm): #form
    class SigMulti(Multi):
        label='Signals: '
        track = twf.FileField(label=' ',
        help_text='Select files (e.g. bedgraph)',
        validator=twf.FileValidator(required=True))
    function = twf.RadioButtonList(label='operation: ',
        options=["log2","log10","sqrt"],
        validator=twc.Validator(required=True))
    format = twf.SingleSelectField(label='Output format: ',
        options=["sql","bedgraph","bigwig","wig","gff"],
        validator=twc.Validator(required=False),
        help_text='Format of the output file')
    assembly = twf.SingleSelectField(label='Assembly: ',
        options=genrep.GenRep().assemblies_available(),
        help_text='Reference genome')
    submit = twf.SubmitButton(id="submit", value="Log2")
meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [{'id': 'track', 'type': 'track', 'required': True}, # input parameters
                {'id': 'assembly', 'type': 'assembly', 'required': True},
                {'id': 'function', 'type': 'radio', 'required': True},
                {'id': 'format', 'type': 'format'}]
out_parameters = [{'id': 'track', 'type': 'track'}] # output parameters
class loga2Plugin(OperationPlugin):
    description = """ Operation on a track."""
    info = {
        'title': 'log2 or log10 or square root',
        'description': description,
        'path': ['track', 'Log2', "log10","sqrt"],
        'track': loga2Form,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        def method(x): # the function that is applied to the track score
            if kw.get('function')=="log2":
                return math.log( x,2) ;
            if kw.get('function')=="log10":
                return math.log10( x) ;
            if kw.get('function')=="sqrt":
                return math.sqrt(x)
        assembly_id = kw.get('assembly') #or None
        assembly = genrep.Assembly(assembly_id)
        chrmeta = assembly.chrmeta
        length=len( kw.get('track') ) # number of tracks
        for i in range( length  ) :
            tname = kw['track'][i]
            tinput = track(tname, chrmeta=kw.get('assembly'))
            (filepath, filename) = os.path.split(tname) # ( path, name of one track)
            (shortname, extension) = os.path.splitext(filename) # (name of one track (without path) , extension)
            modif= kw.get('function')
            if not isinstance(tinput, (tuple,list)):
                tinput = [tinput]
            for trackFILE in tinput:
                if  "score" in trackFILE.fields:
                    if len(str(kw.get('format')))==0:
                        #create a temporary directory, name a new file in this directory
                        output1_name = self.temporary_path(str(shortname)+'_'+str(kw.get('assembly'))+'_'+str(modif)+str(extension))
                        t1_name = track(output1_name,chrmeta=assembly.chrmeta) #define an output track file in the temporary directory
                        #apply method function to the score field trackFILE
                        #create and write in the new file
                        t1_name.write(common.apply(trackFILE.read(),'score',method ) ,mode='overwrite')
                        t1_name.close() #close the output track file
                        trackFILE.close() #close the input track file
                    if len(str(kw.get('format')))>0:
                        output2_name = self.temporary_path(str(shortname)+'_'+str(kw.get('assembly'))+'_'+str(modif)+'.'+ str(kw.get('format') ))
                        t2_name = track(output2_name,chrmeta=assembly.chrmeta)
                        t2_name.write(common.apply(trackFILE.read(),'score',method ) ,mode='overwrite')
                        t2_name.close()
                        trackFILE.close()
