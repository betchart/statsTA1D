import math
import array
import numpy as np
import ROOT as r
from lib import roo
import inputs
import model
from lib.parabola import parabola
from asymmNames import genNameX,genNameY


class fit(object):
    def __init__(self, label, signal, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, genDirPre, d_wbb,
                 quiet = False, templateID=None, defaults = {},
                 log=None, fixSM=False, altData=None, lumiFactor=1.0, only=""):

        np.random.seed(1981)
        for item in ['label','quiet','fixSM','only'] : setattr(self,item,eval(item))
        self.log = log if log else sys.stdout
        if type(R0_) == tuple:
            diffR0_ = R0_[1]
            R0_ = R0_[0]
        else: diffR0_ = None
        prePre = dirIncrement in [0,4,5]
        channels = dict([((lep,part),
                          inputs.channel_data(lep, part, tag, signal, sigPre,
                                              "R%02d" % (R0_ + dirIncrement),
                                              genDirPre, prePre=prePre, templateID=templateID,
                                              d_wbb = d_wbb))
                         for lep in ['el', 'mu']
                         for part in ['top', 'QCD']
                         ])
        channels['gen'] = inputs.channel_data('mu', 'top', tag,
                                              '%s; %s'%(genNameX,genNameY),
                                              sigPrefix = sigPre if dirIncrement in [0,4,5] else '',
                                              dirPrefix=genDirPre, genDirPre=genDirPre, 
                                              prePre = prePre, sampleList=['tt'], rename=False)

        if diffR0_ :
            for lepPart,chan in channels.items():
                if type(lepPart) != tuple: continue
                lep,part = lepPart
                chan.subtract(inputs.channel_data(lep,part,tag,signal,sigPre,
                                                  "R%02d" % (diffR0_ + dirIncrement),
                                                  genDirPre, prePre = prePre ))

        print "###", label
        print>>self.log, "###", label
        self.model = model.topModel(channels)
        self.model.w.arg('lumi_factor').setVal(lumiFactor)
        for k,v in defaults.items(): self.model.w.arg(k).setVal(v)
        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: self.model.w.arg(item).setVal(eval(item))

        self.fitArgs = [r.RooFit.Extended(True), r.RooFit.NumCPU(1),
                        r.RooFit.PrintLevel(-1)]
        self.model.import_data(altData)

        if fixSM: self.doSM()
        else: self.doFit()

    @roo.quiet
    def doSM(self):
        self.model.w.arg('alpha').setConstant()
        nll = self.model.w.pdf('model'+self.only).createNLL(self.model.w.data('data'+self.only), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        minu.setStrategy(2)
        minu.migrad()

    @roo.quiet
    def doFit(self):
        w = self.model.w
        (alpha,alphaE), pll = self.minimize()

        def pllEval(p):
            w.arg('alpha').setVal(p)
            return pll.getVal()
        step = 0.2 * alphaE
        points = [alpha+i*step for i in range(-10,11)]
        pllPoints = [pllEval(p) for p in points]
        iBounds = zip(*sorted((abs(1-v),i) for i,v in enumerate(pllPoints))[:2])[1]
        parb = parabola([(points[i],pllPoints[i]) for i in (10,)+iBounds])

        self.profVal = float(parb.xmin)
        self.profErr = parb.dx(0.5)
        self.scale = w.arg('Ac_tt').getVal()
        self.fit = self.profVal * self.scale
        self.sigma = self.profErr * self.scale

        self.profPLL = pllEval(self.profVal)
        for item in ['pll','points','pllPoints']: setattr(self,item,eval(item))
        self.parbABC = list(parb.ABC)

        if not self.quiet:
            print>>self.log
            print>>self.log, 'Ac: %f +/- %f'% (self.fit, self.sigma)
            print>>self.log
            print>>self.log, "fit alpha: %f +/- %f"%(alpha,alphaE)
            print>>self.log, "parb alpha: %f +/- %f"%(self.profVal, self.profErr)
            print>>self.log, "min PLL:", self.profPLL
            print>>self.log, parb
            for item in ['d_xs_tt','d_xs_wj','factor_elqcd','factor_muqcd']:
                print>>self.log, '\t', roo.str(w.arg(item))
            print>>self.log
            print>>self.log, 'iBounds', iBounds
            for i,(p,v) in enumerate(zip(points,pllPoints)):
                print>>self.log, i, p, v
        return

    @roo.quiet
    def minimize(self):
        w = self.model.w
        nll = w.pdf('model'+self.only).createNLL(w.data('data'+self.only), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        for j in range(10):
            minu.setStrategy(2)
            for i in range(10):
                self.fitstatus = minu.migrad()
                print>>self.log, i + 10*j,
                self.log.flush()
                if not self.fitstatus: break
            if not self.fitstatus: break
            minu.setStrategy(1)
            minu.migrad()
        print>>self.log
        print>>self.log, roo.str(w.arg('alpha'))
        self.NLL = nll.getVal()
        alpha = w.arg('alpha').getVal(), w.arg('alpha').getError()
        pll = nll.createProfile(w.argSet('alpha'))
        print>>self.log, roo.str(nll)
        print>>self.log, roo.str(pll)
        if hasattr(pll,'minimizer'):
            pll.minimizer().setStrategy(2)
        else: pll.minuit().setStrategy(2)
        return alpha,pll
        

    @staticmethod
    def modelItems():
        return ( ['d_xs_%s'%item for item in ['tt','wj','st','dy']] +
                 ['expect_%s_%s'%(lep,s) for lep in ['el','mu'] for s in ['tt','wj','mj','st','dy']] +
                 ['d_lumi','lumi_factor','alpha','factor_elqcd','factor_muqcd'] )

    def ttree(self, truth={}):
        # Note : ROOT and array.array use opposite conventions for upper/lowercase (un)signed
        #         name     array  ROOT  ROOT_typedef
        types = {int    : ("i", "I", "Int_t"),
                 long   : ("l", "L", "Long_t"),
                 float  : ("f", "F", "Float_t"),
                 bool   : ("B", "O", "Bool_t"),
                 str    : ("c", "C", "Char_t")
             }

        genvals = dict([(item,-99999999.) for item in (['fit']+self.modelItems())])
        genvals.update(truth)

        selfStuff = ['label','fit','sigma','NLL','fitstatus','fixSM','points','pllPoints','scale','parbABC']
        selfPairs = [(item,getattr(self,item)) for item in selfStuff]
        modelPairs = [(item,self.model.w.arg(item).getVal()) for item in self.modelItems()]
        
        address = {}
        tree = r.TTree('fitresult','')
        for name,value in selfPairs+modelPairs+[('gen_'+key,val) for key,val in genvals.items()]:
            if type(value) not in [list,tuple]:
                ar,ro,t = types[type(value)]
                address[name] = array.array(ar, value if type(value)==str else [value])
                tree.Branch(name, address[name], "%s/%s"%(name,ro))
            else:
                address[name] = r.std.vector('float')()
                for item in value: address[name].push_back(item)
                tree.Branch(name, address[name])
        tree.Fill()
        return tree

    def ttreeWrite(self, fname, truth={}):
        tfile = r.TFile.Open(fname,'RECREATE')
        tree = self.ttree(truth)
        tree.Write()
        tfile.Close()

