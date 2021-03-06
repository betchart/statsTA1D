#!/usr/bin/env python

import ROOT as r
import sys

class diffROOT(object):
    '''Check for histogram differences in all subdirectories of a pair of files.'''

    def __init__(self,fname1,fname2):
        self.f1 = r.TFile.Open(fname1)
        self.f2 = r.TFile.Open(fname2)
        self.checkdirs(self.f1,self.f2)

    def checkdirs(self,d1,d2):
        '''Recursive directory diff.'''
        keys,not1,not2 = self.compareListsOfKeys( d1, d2 )
        if not1: print "Missing from arg1:", not1
        if not2: print "Missing from arg2:", not2
    
        for item in sorted(keys):
            name = '/'.join([':'.join(d1.GetPath().split(':')[1:]),item])[1:]
            i1 = self.get(self.f1,name)
            i2 = self.get(self.f2,name)
            if not i1 :
                print "nil i1:", name
            elif not i2:
                print "nil i2:", name
            elif type(i1)==r.TDirectoryFile:
                self.checkdirs(i1,i2)
            elif not self.histosEqual(i1,i2):
                print 'diff:', name
            else: pass
        return

    @staticmethod
    def get(d,key) :
        '''Generalization of Get() which handles `;` in names.'''
        dirs = key.split('/')[:-1]
        key = key.split('/')[-1]
        try:
            for d_ in dirs: d = d.GetKey(d_).ReadObj()
            obj = d.GetKey(key).ReadObj()
        except: obj = d.Get(key)
        return obj
    
    @staticmethod
    def compareListsOfKeys(i1,i2):
        '''Intersection and differences of key lists.'''
        KL1 = i1.GetListOfKeys()
        KL2 = i2.GetListOfKeys()
        ks1 = set(k.GetName() for k in KL1)
        ks2 = set(k.GetName() for k in KL2)
        return ( set.intersection(ks1,ks2), 
                 (ks1 & ks2) - ks1, 
                 (ks1 & ks2) - ks2 )
    
    @staticmethod
    def histosEqual(h1,h2):
        '''Histogram equality test for TH1 or TH2.'''
        nbinsX = h1.GetNbinsX()
        nbinsY = h1.GetNbinsY()
        bins = nbinsX == h2.GetNbinsX() and nbinsY == h2.GetNbinsY()
        loX = h1.GetXaxis().GetBinLowEdge(0) == h2.GetXaxis().GetBinLowEdge(0)
        hiX = h1.GetXaxis().GetBinLowEdge(nbinsX+1) == h2.GetXaxis().GetBinLowEdge(nbinsX+1)
    
        pDiff = max([(abs(v1-v2)/(v1+v2) if v1+v2 else (v1+v2))
                     for v1,v2 in [( h1.GetBinContent(iX,iY),
                                     h2.GetBinContent(iX,iY) ) 
                                   for iX in range(nbinsX+2) for iY in range(nbinsY+2)]])
        same = pDiff < 1e-6
        return all([bins,loX,hiX,same])
    


if __name__=='__main__':
    if len(sys.argv)<3 or not all('.root' in v for v in sys.argv[1:]):
        print "Usage: compare_files.py <file1.root> <file2.root>"
        exit()
    diffROOT(sys.argv[1],sys.argv[2])
    print 'Done' # program can hang trying to close ROOT files
