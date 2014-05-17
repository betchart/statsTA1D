import ROOT as r
from lib import roo
import lib
import systematics
from ensembles import ensemble_specs
from ensembles import calibration_specs
from fitroutine import fit
from asymmNames import genNames
import os
import inputs

class measurement(object):
    def __init__(self, label, signal, R0_,
                 doVis=False, evalSystematics=[],
                 ensembles=None, ensSlice=(None,None),
                 calibrations=None, calSlice=(None,None),
                 outDir='output/', templateID=None, only=""):
        os.system('mkdir -p %s' % outDir)
        self.doVis = doVis
        self.outNameBase = (outDir + 
                            '_'.join(label.split(',')) + only +
                            ('_t%03d'%templateID if templateID!=None else ''))

        with open(self.outNameBase + 'SM.log', 'w') as log:
            pars = systematics.central()
            pars['only'] = only
            pars['label'] += 'SM'
            self.SM = fit(signal=signal, R0_=R0_, log=log,
                          fixSM=True, **pars)
            self.SM.model.print_n(log)

        with open(self.outNameBase + '.log', 'w') as log:
            pars = systematics.central()
            pars['log'] = log
            pars['only'] = only
            if templateID!=None: pars['label'] = 'T%03d'%templateID
            self.central = fit(signal=signal, R0_=R0_, templateID=templateID, **pars)
            self.central.model.print_n(log)
            self.central.ttreeWrite(self.outNameBase+'.root')

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in 
                         ['alpha', 'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']])

        if doVis: self.central.model.visualize(printName=self.outNameBase+'.pdf')

        if ensembles: 
            pars = systematics.central()
            pars.update({'signal':signal, 'R0_':R0_, 'log':log})
            for ensPars in ensemble_specs():
                if ensPars['label'] not in ensembles: continue
                ensPars.update({'ensSlice':ensSlice})
                self.ensembles(pars, **ensPars)

        if calibrations:
            pars = systematics.central()
            pars.update({'signal':signal, 'R0_':R0_, 'log':log})
            for calPars in calibration_specs():
                if calPars['which'] not in calibrations: continue
                calPars.update({'calSlice':calSlice})
                self.calibrations(pars, **calPars)

        for sys in systematics.systematics():
            if sys['label'] not in evalSystematics: continue
            pars = systematics.central()
            pars.update(sys)
            fname = self.outNameBase +'_sys_'+ sys['label'] + '.log'
            with open(fname, 'w') as log:
                f = fit(signal=signal, R0_=R0_,
                        quiet=False, defaults=defaults, log=log, **pars)
                f.ttreeWrite(fname.replace('.log','.root'))

    @roo.quiet
    def ensembles(self, pars, alpha=1.0, lumiFactor=1.0, ensSlice=(None,None), Nens=1000, label=''):
        wGen = self.central.model.w
        truth = {'fit': alpha*self.central.scale}
        wGen.arg('alpha').setVal(alpha)
        wGen.arg('lumi_factor').setVal(lumiFactor)
        pars['lumiFactor'] = lumiFactor

        for item in fit.modelItems(): truth[item] = wGen.arg(item).getVal()

        mcstudy = r.RooMCStudy(wGen.pdf('model'),
                               wGen.argSet(','.join(self.central.model.observables+['channel'])),
                               r.RooFit.Binned(True),
                               r.RooFit.Extended(True)
                           )
        mcstudy.generate(Nens,0,True)
        for i in range(Nens)[slice(*ensSlice)]:
            alt = mcstudy.genData(i)
            pars['label'] = '%s_ens%03d'%(label,i)
            pars['quiet'] = True
            with open(self.outNameBase + pars['label'] + '.log', 'w') as log:
                pars['log']=log
                f = fit(altData=alt, **pars)
            f.ttreeWrite(self.outNameBase + pars['label'] + '.root', truth)
            if self.doVis: f.model.visualize(self.outNameBase + pars['label'] + '.pdf')

    @roo.quiet
    def calibrations(self, pars, which='mn', calSlice=(None,None), N=1000, label='', **kwargs):

        sampleList = [c['sample'] for c in calibration_specs() if c['which']==which]
        prePre = pars['dirIncrement'] in [0,4,5]

        args = {
            'signal':pars['signal'],
            'sigPrefix':pars['sigPre'],
            'dirPrefix':"R%02d" % (pars['R0_'] + pars['dirIncrement']),
            'genDirPre':pars['genDirPre'],
            'prePre':prePre,
            'templateID':None,
            'sampleList': sampleList
            }
        alt_channels = dict( [ ((lep,part), inputs.channel_data(lep, part, **args))
                               for lep in ['el','mu'] for part in ['top','QCD']] )

        if 'xsfactor' in kwargs:
            alt_channels[('el','top')].samples['ttalt'].xs *= kwargs['xsfactor']
            alt_channels[('mu','top')].samples['ttalt'].xs *= kwargs['xsfactor']

        # get Ac_phi_ttalt and Ac_y_ttalt
        filePattern = 'data/stats_top_mu_%s.root'
        tag = 'ph_sn_jn_20'
        tfile = r.TFile.Open(filePattern%tag)
        h = lib.get(tfile,'genTopTanhDeltaAbsY_genTopDPtDPhi/'+ sampleList[0])
        Ac_y_ttalt = lib.asymmetry(h.ProjectionX())[0]
        Ac_phi_ttalt = lib.asymmetry(h.ProjectionY())[0]
        tfile.Close()

        model = self.central.model
        model.import_alt_model(alt_channels)
        wGen = model.w

        # bring xs_ttalt to the most consistant value possible
        wGen.arg('d_xs_ttalt').setVal((wGen.arg('expect_mu_tt').getVal() + wGen.arg('expect_el_tt').getVal()) /
                                      (wGen.arg('expect_mu_ttalt').getVal() + wGen.arg('expect_el_ttalt').getVal()) - 1)
        if not (-0.5 < wGen.arg('d_xs_ttalt').getVal() < 1.5):
            print 'ttalt xs invalid! Adjust calibration_specs!'
            exit()
        # not clear how to do the same for factor_*_qcd (equivalent bg representations)

        truth = {'Ac': Ac_y_ttalt if genNames['XL'][3:] in pars['signal'] else Ac_phi_ttalt}
        altItems = ['expect_%s_ttalt'%s for s in ['el','mu']]
        for item in (set(fit.modelItems()+altItems)-set()): truth[item] = wGen.arg(item).getVal()
        truth.update({'Ac_raw_el_model':model.Ac_raw('el','altmodel'),
                      'Ac_raw_mu_model':model.Ac_raw('mu','altmodel')})

        mcstudy = r.RooMCStudy(wGen.pdf('altmodel'),
                               wGen.argSet(','.join(model.observables+['channel'])),
                               r.RooFit.Binned(True),
                               r.RooFit.Extended(True)
                           )
        mcstudy.generate(N,0,True)
        for i in range(N)[slice(*calSlice)]:
            alt = mcstudy.genData(i)
            pars['label'] = '%s_cal%s%03d'%(label,which,i)
            pars['quiet'] = True
            with open(self.outNameBase + pars['label'] + '.log', 'w') as log:
                pars['log']=log
                f = fit(altData=alt, **pars)
            f.ttreeWrite(self.outNameBase + pars['label'] + '.root', truth)
            if self.doVis: f.model.visualize(self.outNameBase + pars['label'] + '.pdf')
