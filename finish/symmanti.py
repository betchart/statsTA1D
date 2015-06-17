#!/usr/bin/env python

import ROOT as r
import lib
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L lib/tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(0.5);

Xl = "#varUpsilon_{t#bar{t}}"
labels = {'tt':"pp #rightarrow t#bar{t}",
          "ttqq":"q#bar{q} #rightarrow t#bar{t}",
          "ttgg":"gg #rightarrow t#bar{t}",
          "ttqg":"qg #rightarrow t#bar{t}",
          "ttag":"#bar{q}g #rightarrow t#bar{t}"
      }


class symmanti(object):

    def __init__(self, props, outName, fName = 'data/stats_top_mu_ph_sn_jn_20.root', dName = 'genTopTanhDeltaAbsY', propOrder = []):
        tFile = r.TFile.Open(fName)
        hists = {}
        smax = 0
        amax = 0
        leg = r.TLegend(0.7,0.2,0.9,0.45)
        leg.SetBorderSize(0)
        leg.SetFillColor(r.kWhite)
        leg.SetTextFont(42)
        for key in (propOrder if propOrder else props):
            h = tFile.Get(dName + '/' + key)
            h.UseCurrentStyle()
            h.Rebin(2)
            h.SetTitle(';%s;(1/#sigma)(#partial#sigma/#partial%s)'%(Xl,Xl))
            h.Scale(1./h.Integral(),'width')
            col,width,style = props[key]
            h.SetLineColor(col)
            h.SetLineWidth(width)
            h.SetLineStyle(style)
            h.GetXaxis().SetNdivisions(5,4,0,False)
            leg.AddEntry(h, labels[key] if key in labels else key, 'L')
            hists[key] = lib.symmAnti(h)
            smax = max(smax, hists[key][0].GetMaximum())
            amax = max(amax, hists[key][1].GetMaximum())
        
        c = r.TCanvas("canvas","", 800,800)
        c.Print(outName + '[')
        for i,(n,(h,_)) in enumerate(hists.items()):
            h.SetMaximum(1.1*smax)
            h.SetMinimum(0)
            h.Draw('hist' + ('' if not i else 'same'))
            leg.Draw()
        c.Print(outName)
        for i,(n,(_,h)) in enumerate(hists.items()):
            h.SetMaximum(1.1*amax)
            h.SetMinimum(-1.1*amax)
            h.Draw('hist' + ('' if not i else 'same'))
            leg.Draw()
        c.Print(outName)
        c.Print(outName + ']')
        tFile.Close()


outNames = ['output/symmAnti.pdf','output/symmAntiAlt.pdf','output/symmAntiCalib.pdf']
props_s = [#{'tt':(r.kBlack, 3,r.kSolid), 'ttqq':(r.kRed,1,r.kSolid), 'ttgg':(r.kBlue,1,r.kSolid), 'ttqg':(r.kViolet,1,r.kSolid), 'ttag':(r.kViolet,1,r.kDashed)},
           {'ttqq':(r.kRed,3,r.kSolid), 'ttgg':(r.kBlack,1,r.kSolid), 'ttqg':(r.kBlue,3,7), 'ttag':(r.kBlue,1,r.kDashed)},
           {'tt':(r.kBlack, 3,r.kSolid), 'calib_mg.pu.sf':(r.kRed, 2, r.kSolid), 'calib_mn.pu.sf':(r.kBlue, 2, r.kSolid)},
           {'tt':(r.kBlack, 3,r.kSolid),
            'calib_R200.pu.sf':(r.kRed, 2, r.kSolid),
            'calib_R2K.pu.sf':(r.kRed, 2, r.kDashed),
            'calib_A200.pu.sf':(r.kBlue, 2, r.kSolid),
            'calib_A2K.pu.sf':(r.kBlue, 2, r.kDashed),
            'calib_L200.pu.sf':(r.kViolet, 2, r.kSolid),
            'calib_ZP.pu.sf':(r.kGray, 2, r.kSolid)}]

orders = [['ttgg','ttqq','ttqg','ttag'],
         [],
         []]

if __name__=="__main__":
    for outName, props, order in zip(outNames, props_s, orders):
        symmanti(props, outName, propOrder = order)
