#!/usr/bin/env python

from model import simpleModel
from fivebin import fivebin
import lib
import ROOT as r
import math
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetPalette(55)
r.gStyle.SetPaintTextFormat(".2f")

class slope(object):
    def __init__(self, ctemp, cdata, a=0.1):
        h = fivebin(100, a, ctemp)
        m = simpleModel(h)

        datas = [fivebin(10000, i*a/2, cdata) for i in [-8,-1,1,8]]
    
        points = []
        for i,d in enumerate(datas):
            m.import_data(d, str(i))
            m.minimize(str(i))
            points.append((m.val('Ac_data%d'%i), m.val('Ac')))

        self.points = points

    def __call__(self):
        return sum(meas/true for true,meas in self.points) / len(self.points)


class slope_hist(object):
    def __init__(self, a):
        nbins = 10
        h = r.TH2D('biases', 'log( A^{meas} / A^{true} ); Template  N_{c} / N;Data  N_{c} / N', nbins, 0, 1, nbins, 0, 1)
        
        for iX in range(1, 1+nbins):
            x = h.GetXaxis().GetBinCenter(iX)
            for iY in range(1, 1+nbins):
                y = h.GetYaxis().GetBinCenter(iY)
                h.SetBinContent(iX, iY, math.log( slope(x, y, a=a)()))
        
        h.SetContour(100)
        self.h = h

    def __call__(self):
        return self.h

if __name__=='__main__':
    c = r.TCanvas('c','c',600,600)
    #c.SetLogz(1)
    outName = 'slope.pdf'
    c.Print(outName + '[')
    for a in [0.1]:
        h = slope_hist(a)()
        h.Draw('colz text')
        c.Print(outName)
    c.Print(outName+']')
