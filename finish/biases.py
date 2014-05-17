#!/usr/bin/env python


import ROOT as r
r.gROOT.SetBatch(1)


def mean(tree,query):
    hist = r.TH1D('hist','',1000,-10,10)
    tree.Draw("(%s)>>hist"%query)
    mean = hist.GetMean()
    hist.Clear()
    del hist
    return mean


class biases(object):

    def __init__(self, tree):
        
        exprs = {"d_xs_tt":"(d_xs_tt-gen_d_xs_tt)",
                 "d_xs_wj":"(d_xs_wj-gen_d_xs_wj)",
                 "d_expect_mu_tt":"(expect_mu_tt-gen_expect_mu_tt)/gen_expect_mu_tt",
                 "d_expect_mu_wj":"(expect_mu_wj-gen_expect_mu_wj)/gen_expect_mu_wj",
                 "d_expect_mu_mj":"(expect_mu_mj-gen_expect_mu_mj)/gen_expect_mu_mj",
                 "d_expect_el_tt":"(expect_el_tt-gen_expect_el_tt)/gen_expect_el_tt",
                 "d_expect_el_wj":"(expect_el_wj-gen_expect_el_wj)/gen_expect_el_wj",
                 "d_expect_el_mj":"(expect_el_mj-gen_expect_el_mj)/gen_expect_el_mj",
                 "Ac_gen (%)":"100*gen_Ac",
                 "Ac (%) ":"100*fit",
                 "Ac_sigma (%)":"100*sigma"
             }
        
        self.values = dict([(label,mean(tree,query)) for label, query in exprs.items()])
    
    def __str__(self):
        return '\n'.join("%s:\t%+0.3f"%(lab,val) for lab,val in sorted(self.values.items()))
    

import sys

if __name__=='__main__':
    if len(sys.argv) < 2:
        print "Usage: biases.py <filename.root>"
        exit()

    tf = r.TFile.Open(sys.argv[1])
    print biases(tf.Get("fitresult"))
    tf.Close()
    exit()
