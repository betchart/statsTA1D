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

            label = "%+d"%(100*e.gen_alpha)
            permil = '#circ#kern[-0.2]{#/}#kern[-0.6]{#lower[0.4]{#circ#circ}}'
            meanbook.fill(alpha,label, 30, e.gen_alpha - 1.5, e.gen_alpha+1.5, title=';#alpha')
            pullbook.fill( (fit-gen_fit)/sigma, label, 30, -5, 5, title = ';#Delta/#sigma')
            sigmbook.fill( 1000*sigma, label, 30, 2.15, 2.85, title = ';#sigma (%s)'%permil)
            nlle = -3365000
            nllm = -3525000
            nlldelta = 30000
            nllsbook.fill( (e.NLL, m.NLL), label, (30, 30), (nlle - nlldelta,nllm - nlldelta), (nlle+nlldelta,nllm+nlldelta), title = ';NLL (e);NLL (#mu)'  )

        for item in ['mean','pull','sigm','nlls']:
            setattr(self,item+'book',eval(item+'book'))

    def plot(self, outName = 'ensembleFits.pdf'):

        c = r.TCanvas('','',800,400)
        c.Divide(2,1)
        c.Print(outName+'[')
        mpointsXYe = []
        ppointsXYe = []
        order = sorted(int(label) for label in self.meanbook)
        for iLab in order:
            lab = "%+d"%iLab

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
        #x1 = r.TF1('x1','x',-4,4)
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


    def plotSigmas(self, outName = 'ensembleSigmas.pdf'):
        c = r.TCanvas()
        c.Print(outName+'[')
        order = sorted(int(label) for label in self.sigmbook)
        total = self.sigmbook['-300'].Clone('combined_sigmas')
        total.Reset()
        for iLab in order:
            lab = "%+d"%iLab
            m = self.sigmbook[lab]
            total.Add(m)
            m.Fit('gaus','QEML')
            c.Print(outName)
        total.Fit('gaus','QEML')
        
        if self.datatrees:
            a = r.TArrow()
            a.SetLineColor(r.kBlue)
            a.SetLineWidth(4)
            x = 1000 * lib.combined_error([self.datatrees['el'].sigma,self.datatrees['mu'].sigma])
            a.DrawArrow(x,3,x,0.4*total.GetMaximum(),0.05,"<")
        c.Update()

        c.Print(outName)
        c.Print(outName+']')

    def plotNlls(self, outName = 'ensembleNlls.pdf'):
        c = r.TCanvas()
        c.Print(outName+'[')
        order = sorted(int(label) for label in self.nllsbook)
        total = self.nllsbook['-300'].Clone('combined_nlls')
        total.Reset()
        for iLab in order:
            lab = "%+d"%iLab
            m = self.nllsbook[lab]
            total.Add(m)
            #m.Fit('gaus','QEML')
            m.Draw('colz')
            c.Print(outName)
        #total.Fit('gaus','QEML')
        total.Draw('colz')

        if self.datatrees:
            m = r.TMarker()
            m.SetMarkerColor(r.kBlack)
            m.SetMarkerStyle(5)
            m.SetMarkerSize(2)
            x = self.datatrees['el'].NLL
            y = self.datatrees['mu'].NLL
            m.DrawMarker(x,y)
        c.Update()

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
    e1D.plot()
    e1D.plotSigmas()
    e1D.plotNlls()
    if tfdata: tfdata.Close()
    tf.Close()
    
