#!/usr/bin/env python

import ROOT as r
from lib.__autoBook__ import autoBook
import lib
import array,math
from itertools import izip

r.gROOT.SetBatch(1)
r.gStyle.SetOptFit(1)
r.gROOT.ProcessLine(".L lib/tdrstyle.C")
r.setTDRStyle()
r.gStyle.SetPalette(1)


class ensemble1D(object):
    def __init__(self, tree, tree2=None, datatrees=None):
        if datatrees:
            for d in datatrees.values(): d.GetEntry(0)
        self.datatrees = datatrees
        meanbook = autoBook("means")
        pullbook = autoBook("pulls")
        sigmbook = autoBook("sigma")
        nllsbook = autoBook("nlls")

        for e,m in izip(tree, tree2 if tree2 else tree):
            fit, sigma = lib.combined_result([(e.fit,e.sigma),(m.fit,m.sigma)])
            alpha = fit / e.fit * e.alpha
            gen_fit = e.gen_alpha *  e.fit / e.alpha

            label = "%d"%(100*e.lumi_factor)
            permil = '#circ#kern[-0.2]{#/}#kern[-0.6]{#lower[0.4]{#circ#circ}}'
            meanbook.fill(alpha,label, 30, e.gen_alpha - 1.5, e.gen_alpha+1.5, title=';#alpha')
            pullbook.fill( (fit-gen_fit)/sigma, label, 30, -5, 5, title = ';#Delta/#sigma')
            sigmbook.fill( sigma, label, 200, 0, 6./1000, title = ';#sigma')

        for item in ['mean','pull','sigm']:
            setattr(self,item+'book',eval(item+'book'))

    def plot(self, outName = 'ensembleFits_lumi.pdf'):

        c = r.TCanvas('','',800,400)
        c.Divide(2,1)
        c.Print(outName+'[')
        mpointsXYe = []
        ppointsXYe = []
        order = sorted(int(label) for label in self.meanbook)
        for iLab in order:
            lab = "%d"%iLab

            c.cd(1)
            m = self.meanbook[lab]
            m.Fit('gaus','QEML')
            mgaus = m.GetFunction('gaus')
            mpointsXYe.append((iLab/100., mgaus.GetParameter(1) - iLab/100., mgaus.GetParError(1)))

            c.cd(2)
            p = self.pullbook[lab]
            p.Fit('gaus','QEML')
            pgaus = p.GetFunction('gaus')
            ppointsXYe.append((iLab/100., pgaus.GetParameter(2), pgaus.GetParError(2)))

            c.Print(outName)

        mxye = [array.array('d', i) for i in zip(*mpointsXYe)]
        mgraph = r.TGraphErrors(len(mpointsXYe), mxye[0], mxye[1], array.array('d',[0]*len(mpointsXYe)), mxye[2])
        mgraph.SetTitle(';#alpha_{gen};#LT#alpha#GT - #alpha_{gen}')
        mgraph.Fit("pol1")
        delta = max(abs(mgraph.GetYaxis().GetXmin()), abs(mgraph.GetYaxis().GetXmax()))
        mgraph.SetMinimum(-delta)
        mgraph.SetMaximum(delta)
        c.cd(1)
        mgraph.Draw("AP")
        x1 = r.TF1('x1','0',-4,4)
        x1.SetLineColor(r.kBlack)
        x1.SetLineWidth(1)
        x1.Draw("same")

        pxye = [array.array('d', i) for i in zip(*ppointsXYe)]
        pgraph = r.TGraphErrors(len(ppointsXYe), pxye[0], pxye[1], array.array('d',[0]*len(ppointsXYe)), pxye[2])
        pgraph.SetTitle(';#alpha_{gen};#LT#Delta/#sigma#GT')
        pgraph.Fit("pol0")
        delta = max(abs(1-pgraph.GetYaxis().GetXmin()), abs(1-pgraph.GetYaxis().GetXmax()))
        pgraph.SetMinimum(1-delta)
        pgraph.SetMaximum(1+delta)
        c.cd(2)
        pgraph.Draw("AP")
        one = r.TF1('one','1', -4, 4)
        one.SetLineColor(r.kBlack)
        one.SetLineWidth(1)
        one.Draw('same')

        c.Print(outName)

        c.Print(outName+']')


    def plotSigmas(self, outName = 'ensembleSigmas_lumi.pdf'):
        r.TGaxis.SetMaxDigits(3)
        c = r.TCanvas()
        c.SetRightMargin(0.11)
        c.SetLeftMargin(0.09)
        c.Print(outName+'[')
        mpointsXYe = []
        order = sorted(int(label) for label in self.sigmbook)
        for iLab in order:
            lab = "%d"%iLab
            m = self.sigmbook[lab]
            m.Fit('gaus','QEML')
            mgaus = m.GetFunction('gaus')
            mpointsXYe.append((iLab/100., mgaus.GetParameter(1), mgaus.GetParError(1)))
            c.Print(outName)

        c.SetLogx(1)
        c.SetLogy(1)
        c.SetLeftMargin(0.15)
        c.SetRightMargin(0.05)
        scale = next(y for x,y,e in mpointsXYe if x==1)
        ryePoints = [(x,y/scale,e/scale) for x,y,e in mpointsXYe]
        rye = [array.array('d', i) for i in zip(*ryePoints)]
        mgraph = r.TGraphErrors(len(mpointsXYe), rye[0], rye[1], array.array('d',[0]*len(mpointsXYe)), rye[2])
        mgraph.GetXaxis().SetLimits(0.1,10)
        mgraph.SetMinimum(0.1)
        mgraph.SetMaximum(10)
        mgraph.SetTitle(';luminosity / 19.59fb^{-1}  ;#sigma / #hat{#sigma}')
        mgraph.Draw("AP")
        expected = r.TF1('invsqrt','1./sqrt(x)', 0.1,10)
        expected.Draw('same')
        c.Print(outName)
        c.Print(outName+']')

if __name__=='__main__':
    import sys
    if len(sys.argv) < 2:
        print 'Usage: ensemble1D.py <ensemblefile.root> [<datafile.root>]'
        exit()

    tf = r.TFile.Open(sys.argv[1])
    tree = tf.Get('fitresult')
    etree = tf.Get('elfitresult')
    mtree = tf.Get('mufitresult')
    tfdata = None if len(sys.argv)<3 else r.TFile.Open(sys.argv[2])
    datatrees = {'':tfdata.Get('fitresult'), 'el':tfdata.Get('elfitresult'), 'mu':tfdata.Get('mufitresult') } if tfdata else None
    e1D = ensemble1D(etree,mtree,datatrees)
    #e1D.plot()
    e1D.plotSigmas()
    if tfdata: tfdata.Close()
    tf.Close()
    
