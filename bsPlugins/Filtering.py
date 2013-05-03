from bsPlugins import *
from bbcflib.bFlatMajor import stream as gm_stream
from bbcflib import btrack as track
from bbcflib import genrep
from bbcflib.btrack import track
from bbcflib.bFlatMajor import common
from bbcflib.btrack import FeatureStream
from bbcflib.bFlatMajor.common import split_field,apply,score_threshold
import re
import sys
# import toscawidget2 modules in order to build forms
import tw2.forms as twf
dicoAssemblyChr={ }
for el in genrep.GenRep().assemblies_available():
    listeChr=[]
    for elchr in genrep.Assembly(el[0]).chrmeta.keys():
        listeChr.append(elchr)
    dicoAssemblyChr.update({el[0]:listeChr})
class FilteringForm(BaseForm):
    class SigMulti(twb.BsMultiple):
        label='Signals: '
        track = twb.BsFileField(label=' ',
        help_text='Select files (e.g. bedgraph)',
        validator=twb.BsFileFieldValidator(required=True))
    minimum_score = twf.TextField(label='Minimum score: ',
        validator=twc.IntValidator(required=False))
    maximum_score = twf.TextField(label='Maximum score: ',
        validator=twc.IntValidator(required=False))
    minimum_length = twf.TextField(label='Minimum length: ',
        validator=twc.IntValidator(required=False))
    maximum_length = twf.TextField(label='Maximum length: ',
        validator=twc.IntValidator(required=False))
    class HidingSingleSelectField(twd.HidingTableLayout):
        label='Select: '
        assembly = twd.HidingSingleSelectField(label='Assembly: ' ,prompt_text=None, options=dicoAssemblyChr.keys(),
            mapping={
                'TAIR10': ['TAIR10'],
                'ce6': ['ce6'],
                'NA1000': ['NA1000'],
                'danRer7': ['danRer7'],
                'dm3': ['dm3'],
                'galGal4': ['galGal4'],
                'hg19': ['hg19'],
                'TB40BAC4': ['TB40BAC4'],
                'rheMac2': ['rheMac2'],
                'monDom5': ['monDom5'],
                'mm10': ['mm10'],
                'mm9': ['mm9'],
                'MLeprae_TN': ['MLeprae_TN'],
                'mycoSmeg_MC2_155': ['mycoSmeg_MC2_155'],
                'mycoTube_H37RV': ['mycoTube_H37RV'],
                'tbR25': ['tbR25'],
                'Snake': ['Snake'],
                'ornAna1': ['ornAna1'],
                'panTro3': ['panTro3'],
                'sacCer2': ['sacCer2'],
                'pombe': ['pombe'],
                'tetNig2': ['tetNig2'],
                'vibrChol1': ['vibrChol1'],
            })
        TAIR10 = twf.MultipleSelectField(label='Chromosome(s): ', options=dicoAssemblyChr['TAIR10'])
        ce6 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['ce6'])
        mm9 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['mm9'])
        mycoSmeg_MC2_155 = twf.MultipleSelectField(label='Chromosome(s): ', options=dicoAssemblyChr['mycoSmeg_MC2_155'])
        rheMac2= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['rheMac2'])
        NA1000 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['NA1000'])
        vibrChol1 = twf.MultipleSelectField(label='Chromosome(s): ', options=dicoAssemblyChr['vibrChol1'])
        panTro3= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['panTro3'])
        NA1000 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['NA1000'])
        galGal4 = twf.MultipleSelectField(label='Chromosome(s): ', options=dicoAssemblyChr['galGal4'])
        pombe= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['pombe'])
        mm10 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['mm10'])
        sacCer2 = twf.MultipleSelectField(label='Chromosome(s): ', options=dicoAssemblyChr['sacCer2'])
        MLeprae_TN= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['MLeprae_TN'])
        danRer7 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['danRer7'])
        tbR25 = twf.MultipleSelectField(label='Chromosome(s): ', options=dicoAssemblyChr['tbR25'])
        ornAna1= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['ornAna1'])
        tetNig2 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['tetNig2'])
        dm3= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['dm3'])
        TB40BAC4= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['TB40-BAC4'])
        Snake= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['Snake'])
        hg19 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['hg19'])
        mycoTube_H37RV= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['mycoTube_H37RV'])
        tetNig2 = twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['tetNig2'])
        monDom5= twf.MultipleSelectField(label='Chromosome(s): ',options=dicoAssemblyChr['monDom5'])
        chromosome_name=[TAIR10,ce6,NA1000,danRer7,dm3,galGal4,hg19,TB40BAC4,rheMac2,monDom5,mm10,mm9,MLeprae_TN,
        mycoSmeg_MC2_155,mycoTube_H37RV,tbR25,Snake,ornAna1,panTro3,sacCer2,pombe,tetNig2,vibrChol1]
        submit = twf.SubmitButton(id="submit", value="Submit")
meta = {'version': "1.0.0",
        'author': "BBCF",
        'contact': "webmaster-bbcf@epfl.ch"}
in_parameters = [{'id': 'track', 'type': 'track', 'required': True, 'multiple':'SigMulti'},
                {'id': 'assembly', 'type': 'assembly', 'required': True, 'single':'HidingSingleSelectField'},
                {'id': 'minimum_score', 'type': 'minimum_score'},
                {'id': 'maximum_score', 'type': 'maximum_score'},
                {'id': 'minimum_length', 'type': 'minimum_length'},
                {'id': 'maximum_length', 'type': 'maximum_length'},
                {'id': 'chromosome_name', 'type': 'chromosome_name', 'multiple': True}]
out_parameters = [{'id': 'output', 'type': 'track'}]
class FilteringPlugin(BasePlugin):
    description = """ Filter a track."""
    info = {
        'title': 'Apply a filter on a given track',
        'description': ' Apply a filter on a given track (by score, by name, by length, ...) ',
        'path': ['Files', 'Filtering'],
        'output': FilteringForm,
        'in': in_parameters,
        'out': out_parameters,
        'meta': meta,
        }
    def __call__(self, **kw):
        def length_threshold(FilinTrack, thresholdUp, thresholdDown):
            """
            Filter the features of a track which a length=abs(start-end) above or below a certain threshold.
            :param stream: FeatureStream, or list of FeatureStream objects.
            :param threshold: (float) threshold above/below which features are retained
            :rtype: filtered row of a track
            """
            thresholdUp = int(thresholdUp or sys.maxint)
            thresholdDown = int(thresholdDown or -sys.maxint)

            for element2 in FilinTrack.read():
                flelement_start = int(element2[1])
                flelement_end = int(element2[2])
                if abs(flelement_start-flelement_end) > int(thresholdDown) :
                    if abs(flelement_start-flelement_end) < int(thresholdUp) :
                        yield element2
        def filter_score(file1):
            """
            Filter the features of a track which a score is above or below a certain threshold.
            :param stream: FeatureStream.
            :rtype: FeatureStream
            """
            fl_maximum_score =  float(kw.get('maximum_score') or sys.maxint)
            fl_minimum_score =  float(kw.get('minimum_score') or -sys.maxint)
            file1Threshold=score_threshold(file1 ,threshold= fl_maximum_score , lower=True, strict=True, fields='score'  )# < maximum score
            return  score_threshold(file1Threshold ,threshold= fl_minimum_score , lower=False, strict=True, fields='score'  ) # > minimum score
        def filter_name(Stream2 ):
            """
            Filter the features of a track by chromosome name.
            :param stream: FeatureStream.
            :rtype: FeatureStream
            """
            liste2=[]
            chrname = kw['chromosome_name']
            if chrname != "":
                for el in Stream2 :
                    for chromosome_i in chrname :
                        if chromosome_i in el :
                            liste2.append(el)
                return FeatureStream(liste2,Stream2.fields)
            else:
                return (Stream2)
        assembly_id = kw['assembly']
        assembly = genrep.Assembly(assembly_id)
        chrmeta = assembly.chrmeta
        length=len( kw['track'])    # number of tracks
        for i in range( length  ) :
            tname = kw['track'][i]
            tinput = track(tname, chrmeta=kw.get('assembly'))
            (filepath, filename) = os.path.split(tname)    # ( path, name of one track)
            (shortname, extension) = os.path.splitext(filename)    # (name of one track (without path) , extension)
            if not isinstance(tinput, (tuple,list)):
                tinput = [tinput]
            for trackFILE in tinput:
                if  "score" in trackFILE.fields:
                    output2_name = self.temporary_path( shortname+'_'+kw.get('assembly')+'_filter'+str(extension))
                    out_track = track(output2_name,chrmeta=assembly.chrmeta)
                    #Filter by ' length= |start-end| '
                    track_length = FeatureStream(length_threshold(trackFILE,kw.get('maximum_length'), kw.get('minimum_length')),trackFILE.fields)
                    #Filter by score
                    track_length_score=filter_score(track_length)
                    #Filter by name
                    track_length_score_name=filter_name(track_length_score )
                    out_track.write( track_length_score_name, mode='write')
                    out_track.close()
                    trackFILE.close()
