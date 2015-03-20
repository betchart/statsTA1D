#!/usr/bin/env python
import numpy as np
import lib
from lib.__autoBook__ import autoBook
from itertools import izip
import ROOT as r
import math

unfold = False

import matplotlib
class bias_plot(object):
    def __init__(self, trees, treesA, treesC):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        just = -0.6

        fs = 14
        lw = 1.3
        fig = plt.figure(figsize=(6.5,6.5))
        ax = fig.add_subplot(111)
        ax.set_ylim(-4.5,0.5 +[0,1][unfold])
        ax.set_xlim(-2,2)
        ax.set_xlabel(r'$A_c^y (\%)$', fontsize=fs)
        #ax.set_xlabel(r'$A_c^y (\%)$ : Calculated', fontsize=fs)
        #ax.set_aspect('equal')
        ax.get_yaxis().set_visible(False)
        
        bookC = autoBook('C')
        for e,m in izip(*treesC):
            fit,sigma = lib.combined_result([(e.fit,e.sigma),(m.fit,m.sigma)])
            label = e.label[4:6]
            if label in ['R.','L.','A.','ZP','A2','R2','L2','mg']: continue
            bookC.fill( e.gen_Ac*100, 'gen_'+label, 1000, -2, 2)
            bookC.fill( fit*100, "mean_"+label, 100, (e.gen_Ac-e.scale)*100, (e.gen_Ac+e.scale)*100)

        lw = 2
        if unfold:
            ax.errorbar( 0.5, 0, xerr=math.sqrt(0.7**2+0.6**2), marker='.', markersize=15, mfc='k', mec='k', color='r', linewidth=lw, )
            ax.text(just, 0, 'CMS 8TeV (unfold)', ha='right')

        fit,sigma = lib.combined_result([(tree.fit,tree.sigma) for tree in trees])
        print fit,sigma
        sigmaboth = 0.0042
        ax.axvspan( -100, -99, alpha=0.2, fc='k', hatch='', label=r'$68\%\ \mathrm{CI}$')
        ax.axvspan( 100*(fit-sigmaboth), 100*(fit+sigmaboth), alpha=0.1, fc='k', hatch='')
        ax.axvspan( 100*(fit-2*sigmaboth), 100*(fit+2*sigmaboth), alpha=0.1, fc='k', ec='k', hatch='', label=r'$95\%\ \mathrm{CI}$')
        ax.errorbar( 100*fit, [0,1][unfold], xerr=100*sigmaboth, color='r', marker='.', markersize=15, mfc='k', mec='k', linewidth=lw)
        ax.text(just, [0,1][unfold], r'$\mathrm{CMS\ 8TeV\ (template)}$', ha='right')

        PHerr = 0.0009*100
        PH = (tree.scale*100, PHerr)
        KR = (0.0102*100, 0.0005*100)
        BS = (0.0111*100, 0.0004*100)
        predictions = zip([KR, BS, PH],[(0.75,0,0),(0.5,0,0),(0.2,0.8,0)],[r'$\mathrm{K\"{u}hn\ &\ Rodrigo}$',r'$\mathrm{Bernreuther\ &\ Si}$',r'$\mathrm{POWHEG}$'])
        for i,((f,s),c,L) in enumerate(predictions):
            ax.errorbar( f, -1-i, xerr=s, color='r', linewidth=lw)
            ax.text(just, -1-i, L, ha='right')

        names = {'mn':r'$\mathrm{MC@NLO}$'}
        order = [0]
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
        for i,(g,f,e,l) in enumerate(sorted(zip(cgen,cfit,cerr,clab))):
            ax.errorbar([g], [-4-order[i]], xerr=PHerr, color='r', linewidth=lw )
            #ax.plot([g], [-4-order[i]], 'ob', color=(0,0,0.85), )
            #ax.arrow(f, -4-order[i], g-f, 0, color=(0,0,0.85), head_length=0.1, head_width=0.2)
            ax.text(just, -4-order[i], names[l[-2:]], ha='right')
        for k,g,f,e in sorted(zip(clab,cgen,cfit,cerr), key=lambda x: x[1]):
            print k, g, f, e


        ax.legend(loc='upper right', prop={'size':10}).draw_frame(False)

        outName = 'output/result_plot.pdf'
        pp = PdfPages(outName)
        pp.savefig(fig)
        pp.close()
        print "Wrote:", outName


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