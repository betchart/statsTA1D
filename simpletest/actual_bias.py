#!/usr/bin/env python

from model import simpleModel
from fivebin import fivebin
import math
import lib
import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptStat(0)
r.gStyle.SetPalette(55)
r.gStyle.SetPaintTextFormat(".2f")

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


        pars = dict(zip(['middle','central','ain','aout'], classify5(lep_hists['ph'])))
        pars['N'] = 10000
        h = r.TH2D('ph_symm_mods','log (A_{meas} / A_{true});log(c/c_{0});log[(m/c) / (m_{0}/c_{0})]', 11, -0.2, 0.2, 11, -0.2, 0.2)
        for iC in range(h.GetNbinsX()):
            for iM in range(h.GetNbinsY()):
                local_pars = dict(pars)
                cfactor = math.exp(h.GetXaxis().GetBinCenter(1+iC))
                mfactor = math.exp(h.GetYaxis().GetBinCenter(1+iM))
                local_pars.update({'central':pars['central']*cfactor,
                                   'middle':pars['middle'] * mfactor/cfactor})
                tmp = fivebin(**local_pars)
                name = '%s_cm_%d_%d'%(lep, 1+iC, 1+iM)
                m.import_data(tmp, name)
                m.minimize(name)
                Ac = lib.asymmetry(tmp)[0]
                Ac_meas = m.val('Ac')
                h.SetBinContent(1+iC, 1+iM, math.log(Ac_meas / Ac))
                del tmp
        
        c = r.TCanvas('c','',600,600)
        h.Draw('colz text')
        c.Print('%s_iC_iM.pdf'%lep)
        del c

        h = r.TH1D('ph_anti_mods','A_{meas} / A_{true};a^{in} / a^{out}', 50, -1, 1)
        for iX in range(h.GetNbinsX()):
            local_pars = dict(pars)
            x = h.GetXaxis().GetBinCenter(1+iX)
            local_pars.update({'ain':pars['aout'] * x})
            tmp = fivebin(**local_pars)
            name = '%s_aa_%d'%(lep, 1+iX)
            m.import_data(tmp, name)
            m.minimize(name)
            Ac = lib.asymmetry(tmp)[0]
            Ac_meas = m.val('Ac')
            h.SetBinContent(1+iX, Ac_meas / Ac)
            del tmp
        
        c = r.TCanvas('c','',600,600)
        h.Draw('text')
        c.Print('%s_ain_aout.pdf'%lep)
        del c
