#!/usr/bin/env python

import sys
import ROOT as r
from start import systematics
from start.ensembles import calibration_specs
from start.ensembles import ensemble_specs
from start.measurement import measurement
from start.options import opts
from start.inputs import channel_data
import os
from lib import batch

def chunk(L,n):
    return [L[i:i+n] for i in range(0,len(L), n)] if L else []

def chunktuple(T,n):
    return [(i,min(i+n,T[1])) for i in range(T[0],T[1], n)] if T[0]!=T[1] else []

if __name__ == '__main__':
    r.gROOT.SetBatch(1)
    options = opts()

    measure = 0
    if options.partitions==True:
        print ','.join(systematics.partitions)
        exit()
    partitions = options.partitions.split(',')

    allSystematics = [s['label'] for s in systematics.systematics()]
    if options.systematics == True:
        print ','.join(allSystematics)
        exit()
    systs = ([] if not options.systematics else 
             allSystematics if options.systematics=='all' else 
             options.systematics.split(','))

    allEnsembles = [e['label'] for e in ensemble_specs()]
    if options.ensembles == True:
        print ','.join(allEnsembles)
        exit()
    ensembles = ([] if not options.ensembles else
                 allEnsembles if options.ensembles=='all' else
                 options.ensembles.split(','))

    allCalibrations = [c['which'] for c in calibration_specs()]
    if options.calibrations == True:
        print ','.join(allCalibrations)
        exit()
    calibrations = ([] if not options.calibrations else
                    allCalibrations if options.calibrations=='all' else
                    options.calibrations.split(','))

    ensSlice = ((None,None) if not options.ensSlice else tuple(int(i) for i in options.ensSlice.split(':')) if ':' in options.ensSlice else eval("[%s]"%options.ensSlice) if ',' in options.ensSlice else (int(options.ensSlice),1+int(options.ensSlice)))
    calSlice = ((None,None) if not options.calSlice else tuple(int(i) for i in options.calSlice.split(':')) if ':' in options.calSlice else eval("[%s]"%options.calSlice) if ',' in options.calSlice else (int(options.calSlice),1+int(options.calSlice)))
    templates = ((0,0) if not options.templates else tuple(int(i) for i in options.templates.split(':')) if ':' in options.templates else (int(options.templates),1+int(options.templates)))

    if options.batch:
        stack = []
        for part in partitions:
            syschunks = chunk(systs, options.chunk)
            enschunks = chunk(ensSlice, options.chunk) if type(ensSlice)==list else chunktuple(ensSlice, options.chunk)
            calchunks = chunk(calSlice, options.chunk) if type(calSlice)==list else chunktuple(calSlice, options.chunk)
            tmpchunks = chunktuple(templates, options.chunk)
            letter = 'L' if options.XL else 'T'
            if syschunks: stack.extend(["./start.py --X%s --partition %s --systematics %s"%(letter, part, ','.join(s)) for s in syschunks])
            if templates: stack.extend(["./start.py --X%s --partition %s --templates %d:%d"%((letter, part,)+t) for t in tmpchunks])
            if enschunks: stack.extend(["./start.py --X%s --partition %s --ensembles %s --ensSlice "%(letter, part, e) + (','.join(str(s_) for s_ in s) if type(s)==list else "%d:%d"%s) for s in enschunks for e in ensembles])
            if calchunks: stack.extend(["./start.py --X%s --partition %s --calibrations %s --calSlice "%(letter, part, e) + (','.join(str(s_) for s_ in s) if type(s)==list else "%d:%d"%s) for s in calchunks for e in calibrations])
            if not any([syschunks,tmpchunks,enschunks,calchunks]): stack.append("./start.py --X%s --partition %s"%(letter,part))
        for item in ['onlyel','onlymu','nobg','rebin','no3D','twoStage','alttt','sepchan']:
            if getattr(options,item):
                stack = [' --'.join([c,item]) for c in stack]
        print '\n'.join(stack)
        batch.batch(stack, site=options.site)
    else:
        for part in partitions:
            for tID in ([None] if templates[0]==templates[1] else 
                        range(channel_data.nTemplates)[slice(*templates)]):
                mp = systematics.measurement_pars(partition=part, var = 'XL' if options.XL else 'XT')
                mp.update({'doVis':options.visualize,
                           'evalSystematics':systs if tID == None else [],
                           'ensembles':ensembles if tID == None else [],
                           'ensSlice':ensSlice,
                           'calibrations':calibrations if tID == None else [],
                           'calSlice':calSlice,
                           'templateID':tID,
                           'only':'_mu' if options.onlymu else '_el' if options.onlyel else '',
                           'nobg':'nobg' if options.nobg else '',
                           'rebin':options.rebin,
                           'no3D':options.no3D,
                           'twoStage':options.twoStage,
                           'alttt':options.alttt,
                           'sepchan':options.sepchan
                       })
                print mp
                if '_CONDOR_SCRATCH_DIR' in os.environ:
                    mp['outDir'] = os.environ['_CONDOR_SCRATCH_DIR'] + '/output/'
                measurement(**mp)
