#!/usr/bin/env python
import numpy as np
import lib
from lib.__autoBook__ import autoBook
from itertools import izip
import ROOT as r

import matplotlib
class bias_plot(object):
    def __init__(self, trees, treesA, treesC):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        fs = 14
        lw = 1.3
        fig = plt.figure(figsize=(6.5,6.5))
        ax = fig.add_subplot(111)
        ax.set_ylim(-2,2)
        ax.set_xlim(-2,2)
        ax.set_ylabel(r'$A_c^y (\%)$ : Measured', fontsize=fs)
        ax.set_xlabel(r'$A_c^y (\%)$ : Calculated', fontsize=fs)
        ax.set_aspect('equal')
        
        t = np.arange(-2,2,0.01)
        ax.plot(t,t, lw=0.5, color='k')[0].set_dashes([1,1])

        bookA = autoBook('A')
        for e,m in izip(*treesA):
            fit,sigma = lib.combined_result([(e.fit,e.sigma),(m.fit,m.sigma)])
            label = e.label[1:-9]
            #print float(label)*e.scale, fit
            bookA.fill( fit*100, "mean_"+label, 100, (float(label)-1)*e.scale*100, (float(label)+1)*e.scale*100)

        gen = []
        fit = []
        err = []
        for k,v in sorted(bookA.items()):
            gen.append(float(k[5:])*e.scale*100)
            fit.append(v.GetMean())
            err.append(v.GetMeanError())
        ax.errorbar(gen,fit,yerr=err,fmt='o', color='k', mec=(0,0.9,0), mfc='none', label='Extended POWHEG', mew=1)

        bookC = autoBook('C')
        for e,m in izip(*treesC):
            fit,sigma = lib.combined_result([(e.fit,e.sigma),(m.fit,m.sigma)])
            label = e.label[4:6]
            bookC.fill( e.gen_Ac*100, 'gen_'+label, 1000, -2, 2)
            bookC.fill( fit*100, "mean_"+label, 100, (e.gen_Ac-e.scale)*100, (e.gen_Ac+e.scale)*100)

        cgen = []
        cfit = []
        cerr = []
        clab = []
        for k,v in bookC.items():
            if 'gen' in k: continue
            cgen.append(bookC[k.replace('mean','gen')].GetMean())
            cfit.append(v.GetMean())
            cerr.append(v.GetMeanError())
            clab.append(k)
        ax.errorbar(cgen,cfit,yerr=cerr,fmt='.', color=(0,0,0.85), mec='k', label=r'Alternative $\mathrm{t\bar{t}}$ models')
        for k,g,f,e in sorted(zip(clab,cgen,cfit,cerr), key=lambda x: x[1]):
            print k, g, f, e

        fit,sigma = lib.combined_result([(tree.fit,tree.sigma) for tree in trees])
        ax.axhspan( -100, -99, alpha=0.3, fc='k', hatch='', label=r'$(e\oplus\mu)\pm\sigma_{stat}$')
        ax.axhspan( 100*(fit-sigma), 100*(fit+sigma), alpha=0.2, fc='k', hatch='')
        ax.axhspan( 100*(fit-0.0039), 100*(fit+0.0039), alpha=0.15, fc='k', hatch='', label=r'$(e\oplus\mu)\pm\sigma_{stat}\pm\sigma_{sys}$')

        PH = (tree.scale*100, 0.0009*100)
        KR = (0.0102*100, 0.0005*100)
        BS = (0.0111*100, 0.0004*100)
        predictions = zip([PH, KR, BS],[(0.2,0.8,0),(0.75,0,0),(0.5,0,0)],['POWHEG','K&R','B&S'])
        for (f,s),c,L in predictions:
            ax.axvspan( f-s, f+s, alpha=0.6, fc=c, ec=c, label=L)

        ax.legend(loc='upper left', prop={'size':10}).draw_frame(False)

        labelsfonts = {'fontsize':8}
        ax.text(-0.4, -0.15, 'madgraph', labelsfonts, ha='right')
        ax.text(0.15, 0.5, r"$Z'$", labelsfonts, ha='right')
        ax.annotate('right', xy=(0.479041039944,0.418250670293), xytext=(0.4,0.1), arrowprops={'fc':'k', 'width':0.05, 'shrink':0.2, 'headwidth':2}, fontsize=8)
        ax.text(0.55, 0.3, 'mc@nlo', labelsfonts)
        ax.text(0.65, 0.45, 'RIGHT', labelsfonts)
        ax.text(0.5, 0.7, 'left', labelsfonts, ha='right')
        ax.text(1.15, 0.9, 'AXIAL', labelsfonts)
        ax.text(1.65, 1.3, 'axial', labelsfonts)

        output = 'output/bias_plot.pdf'
        pp = PdfPages(output)
        print 'Wrote:', output
        pp.savefig(fig)
        pp.close()


if __name__=="__main__":
    tfileA = r.TFile.Open("XL_full_rebin_twoStage_sepchanA_ens.root")
    treesA = [tfileA.Get(lep+'fitresult') for lep in ['el','mu']]
    tfileC = r.TFile.Open("XL_full_rebin_twoStage_sepchan_cal.root")
    treesC = [tfileC.Get(lep+'fitresult') for lep in ['el','mu']]
    tfile = r.TFile.Open("XL_full_rebin_twoStage_sepchan.root")
    etree = tfile.Get('elfitresult');    etree.GetEntry(0)
    mtree = tfile.Get('mufitresult');    mtree.GetEntry(0)

    bias_plot([etree,mtree], treesA, treesC)
    tfile.Close()
    tfileA.Close()
    tfileC.Close()
