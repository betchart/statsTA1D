#!/usr/bin/env python
import ROOT as r
import numpy as np
import lib

class nll_plot_root(object):
    
    def __init__(self,trees):
        canvas = r.TCanvas()
        graphs = []
        parbs = []
        boxes = []
        markers = [34,28]
        lines = [r.kSolid, r.kDashed]
        top = 2.01
        for i,tree in enumerate(trees):
            graphs.append( r.TGraphErrors(tree.points.size(), np.array('d',tree.points), np.array('d',tree.pllPoints)) )
            graphs[i].SetMarkerStyle(markers[i])
            graphs[i].SetMarkerSize(1.5)
            graphs[i].GetXaxis().SetLimits(-3,3)
            graphs[i].SetMaximum(top)
            graphs[i].SetTitle(';#alpha;-log L')
            if not i: graphs[i].Draw('pAX+')
            boxes.append(r.TBox((tree.fit-tree.sigma)/tree.scale,0,(tree.fit+tree.sigma)/tree.scale,top))
            boxes[i].Draw()
            graphs[i].Draw('pX+')
            parbs.append( r.TF1("parb","%f * x**2 + 2 * %f * x + %f"%tuple(tree.parbABC), tree.points[0], tree.points[-1]))
            parbs[i].SetLineColor(r.kBlue)
            parbs[i].SetLineStyle(lines[i])
            parbs[i].Draw("same")
        gaxis = r.TGaxis(graphs[0].GetXaxis().GetXmin(), 0, graphs[0].GetXaxis().GetXmax(), 0, 100*trees[0].scale*graphs[0].GetXaxis().GetXmin(), 100*trees[0].scale*graphs[0].GetXaxis().GetXmax())
        gaxis.SetTitle("A_{c}^{y} (%)")
        gaxis.Draw()
        canvas.Update()
        raw_input()


import matplotlib
class nll_plot(object):
    def __init__(self,trees):
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages

        sigtot = 0.0042
        top = 2.01
        fs = 20
        lw = 1.3
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax2 = ax.twiny()
        ax.set_ylim(0,top)
        ax.set_xlim(-3,3)
        ax.set_ylabel(r'-log$L$', fontsize=fs)
        ax.set_xlabel(r'$\alpha$', fontsize=fs)
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax2.xaxis.tick_bottom()
        ax2.xaxis.set_label_position('bottom')
        ax2.set_xlabel(r'$A_c^y (\%)$', fontsize=fs)

        tree = trees[0]
        ax2.set_xlim(-300*tree.scale, 300*tree.scale)

        for i,tree in enumerate(trees):
            print tree.fit, tree.sigma
            fit = tree.fit / tree.scale
            sig = tree.sigma / tree.scale
            kwargs = {'length_includes_head':True, 'fc':['k','none'][i],'ec':'k','linestyle':['solid','dashed'][i], 'head_width':0.06, 'head_length':0.03}
            ax.arrow(fit-sig,0.5, 0, -0.5, **kwargs)
            ax.arrow(fit+sig,0.5, 0, -0.5, **kwargs)
            kwargs['linestyle']='solid'
            ax.arrow(fit-sig,0.05, 0, -0.05, **kwargs)
            ax.arrow(fit+sig,0.05, 0, -0.05, **kwargs)
        for i,tree in enumerate(trees):
            ax.plot(tree.points, tree.pllPoints, 'o', markersize=6, mfc=['black','none'][i], mec='black', mew=lw)
        for i,tree in enumerate(trees):
            t = np.arange(tree.points[0], tree.points[-1], 0.01)
            y = [eval("%f*x**2 + 2*%f*x + %f"%tuple(tree.parbABC)) for x in t]
            p5 = np.array([0.5]*len(y))
            ax.plot(t, y, linewidth=lw, linestyle='--'[:i+1], color='black')
            ax.fill_between(t, y, 0.5, where=y<p5, color=(1,1,1), alpha=0.5, hatch='\/'[i], edgecolor='k', linewidth=0.0, linestyle='--'[:i+1])

        ax.plot([0],[100], 'o', linewidth=lw, linestyle='-', color='k', markersize=6, mfc='k', mec='k', mew=lw, label=r'-log$L$, $e$')
        ax.plot([0],[100], 'o', linewidth=lw, linestyle='--',color='k', markersize=6, mfc='none', mec='k', mew=lw, label=r'-log$L$, $\mu$')

        fit,sigma = lib.combined_result([(tree.fit,tree.sigma) for tree in trees])
        ax.axvspan( -100, -99, alpha=0.3, fc='k', hatch='', label=r'$(e\oplus\mu)\pm\sigma_{stat}$')
        ax.axvspan( (fit-sigma)/tree.scale, (fit+sigma)/tree.scale, alpha=0.2, fc='k', hatch='')
        ax.axvspan( (fit-sigtot)/tree.scale, (fit+sigtot)/tree.scale, alpha=0.15, fc='k', hatch='', label=r'$(e\oplus\mu)\pm\sigma_{stat}\pm\sigma_{sys}$')

        PH = (1, 0.0009/tree.scale)
        KR = (0.0102/tree.scale, 0.0005/tree.scale)
        BS = (0.0111/tree.scale, 0.0004/tree.scale)
        predictions = zip([PH, KR, BS],[(0.2,0.8,0),(0.75,0,0),(0.5,0,0)],['POWHEG','K&R','B&S'])
        for (f,s),c,L in predictions:
            ax.axvspan( f-s, f+s, alpha=0.6, fc=c, ec=c, label=L)

        ax.legend(loc='lower left').draw_frame(False)
        #plt.show()

        output = 'output/nll_plot.pdf'
        pp = PdfPages(output)
        print 'Wrote:', output
        pp.savefig(fig)
        pp.close()


if __name__=="__main__":
    tfile = r.TFile.Open("XL_full_rebin_twoStage_sepchan.root")
    etree = tfile.Get('elfitresult')
    etree.GetEntry(0)
    mtree = tfile.Get('mufitresult')
    mtree.GetEntry(0)
    nll_plot([etree,mtree])
    tfile.Close()
