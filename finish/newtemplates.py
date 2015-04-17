#!/usr/bin/env python

from inputs import channel_data
import lib
import ROOT as r
import math
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L lib/tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(r.TStyle().GetErrorX())
#r.tdrStyle.SetPadTopMargin(0.065)
r.TGaxis.SetMaxDigits(3)
r.tdrStyle.SetEndErrorSize(4)
#r.tdrStyle.SetPadRightMargin(0.06)

rebin = True
import systematics

pars = systematics.measurement_pars()
pars.update(systematics.central())
print pars

channels = dict([(lep, channel_data(lep, 'top', tag=pars['tag'],
                                    signal=pars['signal'],
                                    sigPrefix=pars['sigPre'],
                                    dirPrefix="R%02d" % (pars['R0_'] + pars['dirIncrement']),
                                    genDirPre=pars['genDirPre'],
                                    rebin = rebin
                                )) for lep in ['el','mu']])
linetypes = [1, 2]


print channels['el']
print channels['mu']


comps = ['tt']

projections = {}

for lt,(lep,ch) in zip(linetypes,channels.items()):
    for comp in comps: 
        symm = ch.samples[comp].datas[1]
        anti = ch.samples[comp].datas[2]
        [h.SetLineStyle(lt) for h in [symm,anti]]
        
        sf = 1./symm.Integral()
        [h.Scale(sf) for h in [symm,anti]]

        projections[(lep,comp)] = (symm.ProjectionX('symmx'+lep+comp), anti.ProjectionX('antix'+lep+comp))

fn = 'output/template.pdf'
xlabel = '#varUpsilon_{t#bar{t}}^{rec}'
colors = [r.kBlack, r.kRed]

amax = 0.002
c = r.TCanvas()
c.Print(fn+'[')
cmsstamp = r.TText(0.2, 0.87, "CMS Simulation")
cmsstamp.SetNDC()
cmsstamp.SetTextSize(0.9 * cmsstamp.GetTextSize())

text = r.TText()
text.SetTextFont(42)
text.SetTextSize(0.9 * text.GetTextSize())

for j,sublabel in enumerate(['symmetric','antisymmetric']):
    init = False
    leg = r.TLegend(0.7,0.25,0.9,0.4)
    leg.SetBorderSize(0)
    leg.SetFillColor(r.kWhite)
    leg.SetTextFont(42)
    for lep,color in reversed(zip(channels,colors)):
        h = projections[(lep,comp)][j]
        h.UseCurrentStyle()
        h.SetLineColor(color)
        h.SetMarkerColor(color)
        if j:
            h.SetBinError(3,0.00001)
            h.SetMinimum(-amax)
            h.SetMaximum(amax)
        else:
            h.SetMinimum(0)
            h.SetMaximum(0.28)
        h.GetXaxis().SetNdivisions(5,4,0,False)
        h.SetMarkerSize(1.5)
        if lep=='mu': h.SetMarkerStyle(4)
        h.SetLineWidth(4 if lep=='el' else 2)
        h.GetYaxis().SetTitle('probability / bin')
        h.GetXaxis().SetTitle(xlabel)
        h.Draw('same e1' if init else 'e1')
        init = True
        leg.AddEntry(h, "%s+jets" % ('#mu' if lep=='mu' else 'e'), 'lp')
        text.DrawTextNDC(0.2,0.82,sublabel)
        cmsstamp.Draw()
    leg.Draw()
    c.Print(fn)
c.Print(fn+']')

