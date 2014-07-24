import ROOT as r
import lib
r.gROOT.SetBatch(True)

r.gROOT.ProcessLine(".L lib/tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(0.5);

outName = 'output/symmAnti.pdf'
fName = 'data/stats_top_mu_ph_sn_jn_20.root'
tFile = r.TFile.Open(fName)
dName = 'genTopTanhDeltaAbsY'

props = {'tt':(r.kBlack, 3,r.kSolid), 'ttqq':(r.kRed,1,r.kSolid), 'ttgg':(r.kBlue,1,r.kSolid), 'ttqg':(r.kViolet,1,r.kSolid), 'ttag':(r.kViolet,1,r.kDashed)}

hists = {}
smax = 0
amax = 0
for key in tFile.Get(dName).GetListOfKeys():
    h = tFile.Get(dName + '/' + key.GetName())
    h.UseCurrentStyle()
    h.Rebin(2)
    h.SetTitle(';X_{L};(1/#sigma)(#partial#sigma/#partial X_{L})')
    h.Scale(1./h.Integral(),'width')
    col,width,style = props[key.GetName()]
    h.SetLineColor(col)
    h.SetLineWidth(width)
    h.SetLineStyle(style)
    hists[key.GetName()] = lib.symmAnti(h)
    smax = max(smax, hists[key.GetName()][0].GetMaximum())
    amax = max(amax, hists[key.GetName()][1].GetMaximum())
print hists

c = r.TCanvas("canvas","", 800,800)
c.Print(outName + '[')
for i,(n,(h,_)) in enumerate(hists.items()):
    h.SetMaximum(1.1*smax)
    h.SetMinimum(0)
    h.Draw('hist' + ('' if not i else 'same'))
c.Print(outName)
for i,(n,(_,h)) in enumerate(hists.items()):
    h.SetMaximum(1.1*amax)
    h.SetMinimum(-1.1*amax)
    h.Draw('hist' + ('' if not i else 'same'))
c.Print(outName)
c.Print(outName + ']')


tFile.Close()
