#!/usr/bin/env python
import numpy as np
import lib
from lib.__autoBook__ import autoBook
from itertools import izip
import math
import matplotlib

import ROOT as r
r.gROOT.SetBatch(1)

class bias_plot(object):
    def __init__(self, trees, treesA, treesC):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        from matplotlib.ticker import MultipleLocator

        fs = 16
        lw = 1.3
        fig = plt.figure(figsize=(6.5,6.5))
        ax = fig.add_subplot(111)
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(fs)

        ylimabs = 0.8
        ax.set_ylim(-ylimabs,ylimabs)
        ax.set_xlim(-2,2)
        ax.set_ylabel(r'Bias $(\%)$', fontsize=fs+2)
        ax.set_xlabel(r'$A_c^y$ $(\%)$', fontsize=fs+2)
        #ax.set_aspect('equal')
        ax.xaxis.set_minor_locator(MultipleLocator(0.2))
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(0.2))
        ax.tick_params('both', length=10, width=1, which='major')
        ax.tick_params('both', length=5, width=1, which='minor')
        
        t = np.arange(-2,2,0.01)
        ax.plot(t,np.zeros(len(t)), lw=0.5, color='k')[0].set_dashes([1,1])

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
        ax.errorbar(gen,fit-np.array(gen),yerr=err,fmt='o', color='k', mec=(0,0.9,0), mfc='none', label='Extended POWHEG', mew=1, capsize=0, ms=7)

        bookC = autoBook('C')
        for e,m in izip(*treesC):
            fit,sigma = lib.combined_result([(e.fit,e.sigma),(m.fit,m.sigma)])
            label = e.label[4:6]
            bookC.fill( e.gen_Ac*100, 'gen_'+label, 1000, -2, 2)
            bookC.fill( fit*100, "mean_"+label, 100, (e.gen_Ac-2*e.scale)*100, (e.gen_Ac+2*e.scale)*100)

        cgen = []
        cfit = []
        cerr = []
        clab = []
        iLow = []
        iHi = []
        iOther = []
        iSM = []
        fitfile = "output/bias_plot2-fits.pdf"
        canvas = r.TCanvas()
        canvas.Print(fitfile+'[')
        for i, (k,v) in enumerate(x for x in bookC.items() if 'gen' not in x[0]):
            if 'gen' in k: continue
            cgen.append(bookC[k.replace('mean','gen')].GetMean())
            #cfit.append(v.GetMean())
            #cerr.append(v.GetMeanError())
            v.Scale(1./v.Integral())
            v.Fit("gaus","QEML")
            canvas.Print(fitfile)
            cfit.append(v.GetFunction("gaus").GetParameter('Mean'))
            cerr.append(v.GetFunction("gaus").GetParameter('Sigma'))
            clab.append(k)
            print k
            if '.' in k:
                iLow.append(i)
            elif '2' in k:
                iHi.append(i)
            elif k in ['mean_mg','mean_mn']:
                iSM.append(i)
            else: iOther.append(i)
        canvas.Print(fitfile+']')

        #ax.errorbar(cgen,cfit-np.array(cgen),yerr=cerr,fmt='.', color=(0,0,0.85), mec='k', label=r'Alternative $\mathsf{t\bar{t}}$ models')
        ax.errorbar(np.array(cgen)[iOther],(cfit-np.array(cgen))[iOther],yerr=np.array(cerr)[iOther],fmt='s', color=(0,0,0.85), mec='k', capsize=0, label=r"$\mathrm{Z}'$  model", ms=5)
        ax.errorbar(np.array(cgen)[iLow],(cfit-np.array(cgen))[iLow],yerr=np.array(cerr)[iLow],fmt='^', color=(0,0,0.85), mec='b', mfc='none',label=r'200$\mathsf{\,GeV}$ axigluon models', capsize=0, ms=8)
        ax.errorbar(np.array(cgen)[iHi],(cfit-np.array(cgen))[iHi],yerr=np.array(cerr)[iHi],fmt='v', color=(0,0,0.85), mec='k', label=r'2$\mathsf{\,TeV}$ axigluon models', capsize=0, ms=6)
        ax.errorbar(np.array(cgen)[iSM],(cfit-np.array(cgen))[iSM],yerr=np.array(cerr)[iSM],fmt='.', color=(0,0,0.85), mec='k', label=r'Standard model simulations ', capsize=0)

        for k,g,f,e in sorted(zip(clab,cgen,cfit,cerr), key=lambda x: x[1]):
            print k, g, f, e

        fit,sigma = lib.combined_result([(tree.fit,tree.sigma) for tree in trees])
        #ax.axhspan( -100, -99, alpha=0.3, fc='k', hatch='', label=r'$(e\oplus\mu)\pm\sigma_{stat}$')
        #ax.axhspan( 100*(fit-sigma), 100*(fit+sigma), alpha=0.2, fc='k', hatch='')
        #ax.axhspan( 100*(fit-0.0039), 100*(fit+0.0039), alpha=0.15, fc='k', hatch='', label=r'$(e\oplus\mu)\pm\sigma_{stat}\pm\sigma_{sys}$')

        sys = {'mcstat':0.153, 'modeling': 0.017, 'pdf': 0.018, 'scale': 0.136}
        sys_th = math.sqrt(sum(s*s for s in sys.values()))
        ax.axhspan( -sys_th, sys_th, alpha=0.4, fc='white', hatch='//', edgecolor='k', label=r'Modeling systematic uncertainties')


        handles_, labels_ = [],[]
        handles, labels = ax.get_legend_handles_labels()
        for frag in ['Ext','Stand','GeV','TeV','Z','Mod']:
            i = next(i for i,j in enumerate(labels) if frag in j)
            handles_.append(handles[i])
            labels_.append(labels[i])
        ax.legend(handles_, labels_, loc = 'lower left', prop={'size':13}, numpoints=1).draw_frame(False)
        
        #ax.legend(loc='lower left', prop={'size':13}, numpoints=1).draw_frame(False)

        labelsfonts = {'fontsize':13}
        ax.text(-0.4, 0.15, 'MadGraph', labelsfonts, ha='right')
        #ax.text(0.15, 0.4, r"$Z'$", labelsfonts, ha='right')
        #ax.annotate('right', xy=(0.479041039944,0.418250670293), xytext=(0.4,0.1), arrowprops={'fc':'k', 'width':0.05, 'shrink':0.2, 'headwidth':2}, fontsize=8)
        ax.text(0.4, -0.19, 'MC@NLO', labelsfonts, ha='right')
        #ax.text(0.65, 0.45, 'RIGHT', labelsfonts)
        #ax.text(0.5, 0.7, 'left', labelsfonts, ha='right')
        #ax.text(1.15, 0.9, 'AXIAL', labelsfonts)
        #ax.text(1.65, 1.3, 'axial', labelsfonts)

        ax.text( -1.85, 0.7, "CMS", fontsize=20, weight='heavy')
        ax.text( 0.6, 0.82, "19.6$\,\mathsf{fb^{-1}}$ (8 TeV)", fontsize=16)

        plt.subplots_adjust(top=0.95,right=0.95,left=0.15)
        
        output = 'output/bias_plot2.pdf'
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
