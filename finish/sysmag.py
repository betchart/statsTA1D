#!/usr/bin/env python

import sys
import math
import numpy as np
import ROOT as r
import os

from lib.__autoBook__ import autoBook
import lib
import array,math
from itertools import izip

r.gROOT.SetBatch(1)
r.gStyle.SetOptFit(1)
r.gROOT.ProcessLine(".L lib/tdrstyle.C")
r.setTDRStyle()
r.gStyle.SetPalette(1)


class fitresult(object):
    pairs = dict([(item,(item+'_dn',item+'_up')) for item in ['Q','JER','JES','PU','lumi','DY','ST','as','WBB','TT','WJ','QCDe','QCDm']])
    pairs.update(dict([('PD:%02d-%02d'%(i,i+1),('PD_%02d'%i,'PD_%02d'%(i+1))) for i in range(1,53,2)]))
    pairs.update(dict([('muid',('mu2','mu3')),
                       ('mutrig',('mu0','mu1')),
                       ('elid',('el2','el3')),
                       ('eltrig',('el0','el1')),
                       ('PT', ('PT','PT')),
                       ('Modeling', ('_calmn000','_calmn000')),
                       #('Modeling', ('_alttt-mn','_alttt-mn')), #should be this!!!
                       ('MC statistics',('','')),
                   ]))
    labels = sum(pairs.values(),())

    @staticmethod
    def getTrees(tfile):
        return {'el':tfile.Get('elfitresult'),
                'mu':tfile.Get('mufitresult')}

    def __init__(self,base):
        self.sfile = r.TFile.Open('%s_sys_with_model.root'%base)
        self.cfile = r.TFile.Open('%s.root'%base)
        self.tfile = r.TFile.Open('%s_t.root'%base)
        self.strees = self.getTrees(self.sfile)
        self.ctrees = self.getTrees(self.cfile)
        self.ttrees = self.getTrees(self.tfile)
        for t in self.ctrees.values(): t.GetEntry(0)
        self.extract()

    def central(self): return lib.combined_result([(self.ctrees['el'].fit,self.ctrees['el'].sigma),
                                                   (self.ctrees['mu'].fit,self.ctrees['mu'].sigma)])

    def delta(self,e,m): return (lib.combined_result([(e.fit,e.sigma),
                                                      (m.fit,m.sigma)])[0] - 
                                 self.central()[0])

    @property
    def sigma_mcStat(self):
        return math.sqrt( sum([self.delta(e,m)**2 for e,m in izip(self.ttrees['el'],self.ttrees['mu'])]) / self.ttrees['el'].GetEntries() )

    def extract(self):
        self.values = {}
        for e,m in izip(self.strees['el'],self.strees['mu']):
            label = max(os.path.commonprefix([e.label,l]) for l in self.labels)
            if label not in self.labels: continue
            self.values[label] = self.delta(e,m)

        self.order = sorted(self.values, key = lambda k: self.values[k], reverse=True)
        self.pvalues = dict([(name, math.sqrt(0.5*(self.values[a]**2 + self.values[b]**2))) for name,(a,b) in self.pairs.items() if a])
        self.pvalues['MC stat.'] = self.sigma_mcStat
        self.pvalues['PDF'] = math.sqrt(sum(x**2 for L,x in self.pvalues.items() if 'PD' in L))
        for key in list(self.pvalues):
            if 'PD:' in key: del self.pvalues[key]
        self.porder = sorted(self.pvalues, key = lambda k: self.pvalues[k], reverse=True)

        self.pvalues['Total'] = math.sqrt(sum(x**2 for x in self.pvalues.values()))
        self.porder.insert(0,'Total')

    def form(self,key):
        #top = key in self.order[:5]+self.porder[:5]
        top = False
        form = (r"\textbf{% .3f}" if top else "% .3f ").rjust(20)
        d = (self.values[key] if key in self.order else self.pvalues[key])*100
        return (form%d).ljust(20)


if __name__ == '__main__':

    result = fitresult('XL_full_rebin_twoStage_sepchan')

    def formkey(key) :
        subs = [('PD_',r'pdf'),
                ('mu0',r'$\mu$ trig_up'),
                ('mu1',r'$\mu$ trig_dn'),
                ('mu2',r'$\mu$ id_up'),
                ('mu3',r'$\mu$ id_dn'),
                ('el0',r'$e$ trig_up'),
                ('el1',r'$e$ trig_dn'),
                ('el2',r'$e$ id_up'),
                ('el3',r'$e$ id_dn'),
                ('_up',r'${}^{\uparrow}$'),
                ('_dn',r'${}_{\downarrow}$'),
                ('PU',r'pileup'),
                ('lumi',r'\lumi'),
                ('ST',r'$\sigma_{ST}$'),
                ('DY',r'$\sigma_{DY}$'),
                ('mutrig',r'$\mu$ trigger'),
                ('muid', r'$\mu$ ID'),
                ('eltrig',r'$e$ trigger'),
                ('elid',r'$e$ ID'),
                ('PD:',r'pdf '),
                ('as',r'$\alpha_s$'),
                ('WBB',r'$\mathrm{Wb\bar{b}}$'),
                ('PT',r'$p^\Pqt_{\mathrm{T}}$'),
                ('WJ',r'$\delta_{\W j}$'),
                ('TT',r'$\delta_{$\ttbar$}$'),
                ('Modeling', r'$\ttbar$ modeling'),
                ('QCDe',r'$F^e_{\MJ}$'),
                ('QCDm',r'$F^\mu_{\MJ}$'),
                ('Q',r'$Q^2$ scale')
                 ]
        def rep(key,ps):
            if not ps: return key
            return rep(key.replace(*ps[0]), ps[1:])
        return rep(key,subs)

    vspace=r'''
&\\
'''
    caption = r'''
\caption{\label{list_systematics} %s 
due to sources of systematic variations, ordered by decreasing
magnitude.
}
'''%('Uncertainty on $A_c^y$')


    print r'\begin{table}'
    print r'\begin{tabular}{cl}'
    print r'\hline'
    print r'(\%)                                & Systematic\\'
    print r'\hline'
    print r'\hline'
    print '\n'.join(result.form(key) + ' & ' + formkey(key) + r' \\' + (vspace if i%3==0 else '')
                    for i,key in enumerate(result.porder))
#    print '\n'.join(
#        formkey(key).ljust(28) +' & '+ ' & '.join(fr.form(key) for fr in [result]) + r'  \\' + (vspace if i%5==4 else '')
#        for i,key in enumerate(result.porder if True else result.order))
    print r'\hline'
    print r'\end{tabular}'
    print caption
    print r'\end{table}'
    print
    central = result.central()
    Ac = 100 * central[0]
    stat = 100 * central[1]
    sys = 100 * result.pvalues['Total']
    print Ac, stat, sys
    print Ac, math.sqrt(stat**2 + sys**2)
