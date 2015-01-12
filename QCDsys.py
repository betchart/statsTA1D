#!/usr/bin/env python

import lib
from start import inputs
from start.systematics import partitions, measurement_pars
import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)

fracs = {'el':{'tt':0.170,
               'st':0.011,
               'wj':0.039},
         'mu':{'tt':0.402,
               'st':0.025,
               'wj':0.084}}

colors = {'tt':r.kViolet,
          'wj':r.kGreen,
          'st':r.kGray,
          'data':r.kBlack}

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
        y = 0.05
        c = r.TCanvas()
        c.Divide(2,1)
        c.Print(fname + '[')
        for lep,chan in channels.items():
            template = chan.samples['data'].datasX[0].Clone(lep+'template')
            template.Sumw2()
            Ndata = template.Integral()
            print Ndata
            for s,sd in chan.samples.items():
                if s!='data':
                    template.Add(sd.datasX[0], -fracs[lep][s]*Ndata/sd.datasX[0].Integral())
                c.cd(1)
                scale = 1./sd.datasX[1].Integral()
                sd.datasX[1].SetTitle(sd.datasX[0].GetName().replace('_px',''))
                sd.datasX[1].Scale(scale)
                sd.datasX[2].Scale(scale)
                sd.datasX[1].SetMinimum(0)
                sd.datasX[2].SetMinimum(-y)
                sd.datasX[2].SetMaximum(y)
                sd.datasX[1].SetLineColor(colors[s])
                sd.datasX[2].SetLineColor(colors[s])
                sd.datasX[1].Draw()
                c.cd(2)
                sd.datasX[2].Draw()
                c.Print(fname)
            template.Scale(1./template.Integral())
            tsymm,tanti = lib.symmAnti(template)
            c.cd(1)
            tsymm.SetTitle(lep+'_qcd')
            tsymm.SetLineColor(r.kRed)
            tanti.SetLineColor(r.kRed)
            tsymm.SetMinimum(0)
            tanti.SetMinimum(-y)
            tanti.SetMaximum(y)
            tsymm.Draw()
            c.cd(2)
            tanti.Draw()
            c.Print(fname)
            
        c.Print(fname + ']')
        
if __name__=="__main__":
    QCDsys()
