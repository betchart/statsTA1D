#!/usr/bin/env python

from model import simpleModel
import lib
import ROOT as r
r.gROOT.SetBatch(0)

class hists(dict):
    def __init__(self, lep="el"):
        tf = r.TFile.Open("../data/stats_top_%s_ph_sn_jn_20.root"%lep)
        tdir = "R03_genPdfWeights57_signalhists_/fitTopTanhDeltaAbsY_TridiscriminantWTopQCD/"

        samples = ["ph","mg",'mn','ZP','A2K','A200','L200','R2K','R200']
        for s in samples:
            self[s] = self.five1D(tf.Get(tdir + ("tt" if s=="ph" else "calib_%s.pu.sf"%s)))

    @staticmethod
    def five1D(h):
        htemp = h.ProjectionX()
        return htemp.Rebin(5)



def classify5(h):
    assert h.GetNbinsX()==5
    N = h.Integral()
    m = h.GetBinContent(3) / N
    Nmm = h.GetBinContent(1)
    Npp = h.GetBinContent(5)
    Nm = h.GetBinContent(2)
    Np = h.GetBinContent(4)
    c = 1 - (Nmm + Npp) / N
    a_out = (Npp - Nmm) / (Npp + Nmm)
    a_in = (Np - Nm) / (Np + Nm)

    return (m, c, a_in, a_out)

if __name__=="__main__":
    print '\t'.join(['', "m(%)",'c(%)','ain(%)','aout(%)','','Ac(%)', 'Ac_measure(%)', '', 'slope'])
    for lep in ['el','mu']:
        print lep
        lep_hists = hists(lep)
        m = simpleModel(lep_hists['ph'])
        for k,v in sorted(lep_hists.items()):
            m.import_data(v, k)
            m.minimize(k)
            Ac = lib.asymmetry(v)[0]
            Ac_meas = m.val('Ac')
            print k.rjust(6), (4*"% 0.2f  ") % tuple(100*c for c in classify5(v)), 
            print '|\t\t', '% .2f'%(100*Ac),
            print '\t', "% .2f" % (100*Ac_meas),
            print 3*'\t', "% .2f" % (Ac_meas / Ac)

