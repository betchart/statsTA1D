import ROOT as r
import lib
import numpy as np
from asymmNames import genNameX,genNameY,genNames

class sample_data(object):
    def __init__(self, signalDistribution, xs=None, sigma=None,
                 selectionEfficiency=1.0, preselectionFraction=1.0):
        self.xs = xs
        self.xs_sigma = sigma
        self.eff = selectionEfficiency
        self.frac = preselectionFraction
        self.datas = self.format(signalDistribution)

        self.alphaMax = lib.alphaMax(*self.datas[1:3]) if self.datas[0].GetNbinsX() > 1 else 0

        self.datasX = tuple(d.ProjectionX() if d else None for d in self.datas)
        self.datasY = tuple(d.ProjectionY() if d else None for d in self.datas)
        for d in self.datasX + self.datasY + self.datas: d.SetDirectory(0)

    @staticmethod
    def format(sd):
        return ((sd,) + tuple(lib.symmAnti(sd)))

    def subtract(self,other):
        assert self.xs == other.xs
        assert self.xs_sigma == other.xs_sigma
        assert self.frac == other.frac
        assert self.eff >= other.eff

        self.eff -= other.eff
        for group in ['datas','datasX','datasY']:
            for d,od in zip(getattr(self,group),getattr(other,group)):
                d.Add(od,-1)

    def asymmStr(self):
        Ac = tuple([100*f for f in lib.asymmetry(self.datas[0])])
        return ("%. 2f(%.2f)" % Ac).rjust(8)


    def __str__(self):
        return ('data' if not self.xs else
                ';  '.join([('xs: % 8.2f' % self.xs),
                            ('eff: % .5f' % self.eff).rjust(7),
                            ('f: %.4f' % self.frac).rjust(8),
                            ('d: %.4f' % self.xs_sigma if self.xs_sigma else
                             'd: None   ').rjust(8),
                        ]))

    @property
    def key(self): return (self.xs, self.frac)


class channel_data(object):
    __samples__ = ['data', 'wj', 'dy', 'st', 'tt']
    __xs_uncertainty__ = {'tt': 1.0, 'wj': 2.0, 'st': 0.04, 'dy': 0.04}
    nTemplates = 1000

    def __init__(self, lepton, partition, tag = 'ph_sn_jn_20',
                 signal="", sigPrefix="", dirPrefix="R04", genDirPre="R01",
                 prePre = False, templateID=None, d_wbb=0, sampleList=[], rename=True, rebin=False, no3D=False, only3D=False, fullDir=None, Rst=None):
        filePattern="data/stats_%s_%s_%s.root"
        tfile = r.TFile.Open(filePattern % (partition, lepton, tag))
        self.templateID = templateID
        self.lepton = lepton
        self.lumi = tfile.Get('lumiHisto/data').GetBinContent(1)
        self.lumi_sigma = 0.05
        self.rebin = rebin
        self.no3D = no3D
        self.only3D = only3D
        self.Rst = Rst

        def full(pf) :
            return next((ky.GetName() + '/' for ky in tfile.GetListOfKeys()
                         if pf == ky.GetName().split('_')[0]),
                        '')
        fullDirName = fullDir if fullDir!=None else full(dirPrefix)

        paths = (fullDirName + sigPrefix + signal,
                 fullDirName + signal)

        prepaths = (full(genDirPre) + 
                    (sigPrefix if prePre else '') + 
                    '%s; %s/'%(genNameX,genNameY),
                    'meweighted/')

        self.samples = {}
        for s in (sampleList or self.__samples__):
            self.add(s, tfile, paths, prepaths, sname=('ttalt' if sampleList and rename else ''))
        if d_wbb:
            self.add( 'Wbb', tfile, paths, prepaths)
            for group in ['datas','datasX','datasY']:
                for d,od in zip(getattr(self.samples['wj'],group),
                                getattr(self.samples['Wbb'],group)):
                    d.Add(od,d_wbb)
            del self.samples['Wbb']
        tfile.Close()

    def add(self, s, tfile, paths, prepaths, sname=''):
        def get(s,ps):
            return next(iter(filter(None, [lib.get(tfile,p+s) for p in ps])), None)

        pre = get( s, prepaths)
        if not pre and not s == 'data': return

        data = get('/' + s,paths).Clone(self.lepton + '_' + s)
        data.SetDirectory(0)
        if self.Rst and s=='st':
            R = self.Rst
            T = get( 't', prepaths).Integral()
            T_ = get( 't_', prepaths).Integral()
            t = get('/'+ 't', paths)
            t_= get('/'+ 't_', paths)
            #print T, T_, (T+T_), R
            TT = R*(T+T_)/(R+1)
            TT_ = (T+T_) / (R+1)
            #print TT, TT_, (T+T_), TT/TT_
            #print data.Integral()
            data.Add(t_, (TT_-T_) / T_)
            data.Add(t, (TT-T) / T)
            #print data.Integral()
        if self.rebin: data.RebinX(5)
        if self.no3D: data.RebinY(5)
        elif self.only3D: data.RebinX(data.GetNbinsX())
        if s not in ['data'] or 'QCD' in tfile.GetName() : self.jiggle(data)
        #if s not in ['data']: self.jiggle(data)

        xs = tfile.Get('xsHisto/' + s).GetBinContent(1) if s != 'data' else None
        if s=='tt': xs/=4.0 # magic number hack fix for incorrect xs in data files.  4 for 4 sources of tt: {gg, qq_, gq, gq_}
        delta = (self.__xs_uncertainty__[s[:2]]
                 if s[:2] in self.__xs_uncertainty__ else None)

        named = \
            {'selectionEfficiency': (data.Integral() / pre.Integral() if pre else 0),
             'preselectionFraction': (1.0 if s[:2] != 'tt' else
                                      pre.Integral() / get( 'tt', prepaths ).Integral())
             }
        self.samples[sname or s] = sample_data(data, xs, delta, **named)

    def jiggle(self,hist):
        if self.templateID == None: return
        factor = hist.GetEffectiveEntries() / hist.Integral()
        for iZ in range(1,1+hist.GetNbinsZ()):
            for iX in range(1,1+hist.GetNbinsX()):
                for iY in range(1,1+hist.GetNbinsY()):
                    hist.SetBinContent(iX,iY,iZ,np.random.poisson(hist.GetBinContent(iX,iY,iZ)*factor, self.nTemplates)[self.templateID]/factor)
        return

    def subtract(self, other):
        for s,samp in self.samples.items():
            samp.subtract(other.samples[s])

    def asymmStr(self):
        return '\n'.join('%s:: %s' % (s.rjust(5), data.asymmStr())
                         for s, data in sorted(self.samples.items(),
                                               key=lambda (s, d): d.key,
                                               reverse=True))

    def __str__(self):
        return ('%s : %.2f/pb  (%.2f)\n' % (self.lepton, self.lumi, self.lumi_sigma) +
                '\n'.join('%s:: %s' % (s.rjust(5), str(data))
                          for s, data in sorted(self.samples.items(),
                                                key=lambda (s, d): d.key,
                                                reverse=True)
                          if data.xs))


if __name__ == '__main__':
    import sys
    from systematics import partitions, measurement_pars
    i,j = [int(k) for k in sys.argv[1:3]] if len(sys.argv)>2 else (0,0)
    print '#', partitions[j]
    
    pars = measurement_pars(partitions[j])
    R0_,diffR0_ = pars['R0_'] if type(pars['R0_'])==tuple else (pars['R0_'],None)

    channels = dict([(lep,
                      channel_data(lep,
                                   'top',
                                   signal=pars['signal'],
                                   dirPrefix="R%02d" % R0_)) for lep in
                     ['el','mu'] ])
    if diffR0_:
        for lep,chan in channels.items():
            chan.subtract(channel_data(lep,'top',signal=pars['signal'],
                                       dirPrefix="R%02d" % diffR0_))
    print
    print channels['el']
    print
    print channels['mu']
