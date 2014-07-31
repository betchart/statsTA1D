#!/usr/bin/env python

import ROOT as r
from lib.__autoBook__ import autoBook
import array,math
from itertools import izip

r.gROOT.SetBatch(1)
r.gStyle.SetOptFit(1)

class ensemble1D(object):
    def __init__(self, tree, tree2=None):
        meanbook = autoBook("means")
        pullbook = autoBook("pulls")

        for e,m in izip(tree, tree2 if tree2 else tree):
            esig2inv = 1./e.sigma**2
            msig2inv = 1./m.sigma**2 if tree2 else 0
            sig2inv = esig2inv + msig2inv
            alpha = (e.alpha * esig2inv + m.alpha * msig2inv ) / sig2inv
            fit = (e.fit * esig2inv + m.fit * msig2inv ) / sig2inv
            gen_fit = e.gen_alpha *  e.fit / e.alpha
            sigma = math.sqrt(1./sig2inv)
            label = "%+d"%(100*e.gen_alpha)
            meanbook.fill(alpha,label, 30, e.gen_alpha - 1.5, e.gen_alpha+1.5, title='mean;alpha')
            pullbook.fill( (fit-gen_fit)/sigma, label, 30, -5, 5, title = 'pull;(fit-gen_fit)/sigma')

        c = r.TCanvas()
        c.Divide(2,1)
        outName = 'ensembleFits.pdf'
        c.Print(outName+'[')
        mpointsXYe = []
        ppointsXYe = []
        order = sorted(int(label) for label in meanbook)
        for iLab in order:
            lab = "%+d"%iLab

            c.cd(1)
            m = meanbook[lab]
            m.Fit('gaus','QEM')
            mgaus = m.GetFunction('gaus')
            mpointsXYe.append((iLab/100., mgaus.GetParameter(1) - iLab/100., mgaus.GetParError(1)))

            c.cd(2)
            p = pullbook[lab]
            p.Fit('gaus','QEM')
            pgaus = p.GetFunction('gaus')
            ppointsXYe.append((iLab/100., pgaus.GetParameter(2), pgaus.GetParError(2)))

            c.Print(outName)

        mxye = [array.array('d', i) for i in zip(*mpointsXYe)]
        mgraph = r.TGraphErrors(len(mpointsXYe), mxye[0], mxye[1], array.array('d',[0]*len(mpointsXYe)), mxye[2])
        c.cd(1)
        mgraph.Draw("AP")
        #x1 = r.TF1('x1','x',-4,4)
        x1 = r.TF1('x1','0',-4,4)
        x1.SetLineColor(r.kRed)
        x1.SetLineWidth(1)
        x1.Draw("same")

        pxye = [array.array('d', i) for i in zip(*ppointsXYe)]
        pgraph = r.TGraphErrors(len(ppointsXYe), pxye[0], pxye[1], array.array('d',[0]*len(ppointsXYe)), pxye[2])
        c.cd(2)
        pgraph.Draw("AP")
        one = r.TF1('one','1', -4, 4)
        one.SetLineColor(r.kRed)
        one.SetLineWidth(1)
        one.Draw('same')

        c.Print(outName)

        c.Print(outName+']')






if __name__=='__main__':
    import sys
    if len(sys.argv) < 2:
        print 'Usage: ensemble1D.py <ensemblefile.root>'
        exit()
    
    tf = r.TFile.Open(sys.argv[1])
    tree = tf.Get('fitresult')
    etree = tf.Get('elfitresult')
    mtree = tf.Get('mufitresult')
    ensemble1D(etree,mtree)
    tf.Close()
    
