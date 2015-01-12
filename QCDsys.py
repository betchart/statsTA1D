#!/usr/bin/env python

import lib
from start import inputs
from start.systematics import partitions, measurement_pars
import ROOT as r
r.gROOT.SetBatch(1)

fracs = {'el':{'tt':0.635,
               'st':0.043,
               'wj':0.151},
         'mu':{'tt':0.712,
               'st':0.049,
               'wj':0.173}}

class QCDsys(object):

    def __init__(self):
        i,j = 0,0
        
        pars = measurement_pars(partitions[j])
        R0_,diffR0_ = pars['R0_'] if type(pars['R0_'])==tuple else (pars['R0_'],None)
    
        channels = dict([(lep,
                          inputs.channel_data(lep,
                                              'QCD',
                                              signal=pars['signal'], rebin=True,
                                              dirPrefix="R%02d" % R0_)) for lep in
                         ['el','mu'] ])


        fname = 'qcdtemplates.pdf'
        c = r.TCanvas()
        c.Divide(2,1)
        c.Print(fname + '[')
        for lep,chan in channels.items():
            template = chan.samples['data'].datasX[0].Clone(lep+'template')
            template.Sumw2()
            Ndata = template.Integral()
            for s,sd in chan.samples.items():
                c.cd(1)
                sd.datasX[1].SetMinimum(0)
                sd.datasX[1].Draw()
                c.cd(2)
                sd.datasX[2].Draw()
                c.Print(fname)
                if s!='data':
                    template.Add(sd.datasX[0], -fracs[lep][s]*Ndata/sd.datasX[0].Integral())
            tsymm,tanti = lib.symmAnti(template)
            c.cd(1)
            tsymm.SetMinimum(min(1.1*tsymm.GetMinimum(),0))
            tsymm.Draw()
            c.cd(2)
            tanti.Draw()
            c.Print(fname)
            
        c.Print(fname + ']')
        
if __name__=="__main__":
    QCDsys()
