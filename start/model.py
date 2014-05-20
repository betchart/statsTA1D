import sys
import lib
from lib import roo
import math
import ROOT as r
from asymmNames import genNameX,genNameY,genNames


class topModel(object):
    @roo.quiet
    def __init__(self, channelDict, var='XL', w=None, quiet=True):

        self.fixFractions = False

        leptons = ['el', 'mu']
        observables = [var,'tridiscr']

        channels = dict((L, channelDict[(L,'top')]) for L in leptons)
        channels_qcd = dict((L + 'qcd', channelDict[(L, 'QCD')]) for L in leptons)
        gen = channelDict['gen']

        if not w: w = r.RooWorkspace('Workspace')

        for item in ['gen', 'channels', 'channels_qcd','quiet',
                     'observables', 'w']: setattr(self, item, eval(item))

        for item in ['xs_lumi', 'efficiencies', 'shapes', 'qcd', 'asymmetry',
                     'model', 'expressions']: getattr(self, 'import_' + item)(w)

        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: w.arg(item).setConstant()

    def import_xs_lumi(self, w):
        roo.factory(w, "d_lumi[0,-0.2,0.2]")
        roo.factory(w, 'lumi_factor[1.0]')
        w.arg('lumi_factor').setConstant()
        for L, channel in self.channels.items() + self.channels_qcd.items():
            roo.wimport_const(w, 'lumi_%s_hat' % L, channel.lumi)
            if 'qcd' in L: roo.factory(w, "expr::lumi_%s('(1+@0)*@1', {d_lumi, lumi_%s_hat})" % (L, L))
            else:          roo.factory(w, "expr::lumi_%s('(1+@0)*@1*@2', {d_lumi, lumi_%s_hat, lumi_factor})" % (L, L))

        xs_constraints = dict([(samp[:2], data.xs)
                               for samp, data in self.channels['el'].samples.items()
                               if data.xs > 0])

        for sample, xs in xs_constraints.items():
            roo.wimport_const(w, 'xs_%s_hat' % sample, xs)
            roo.factory(w, "d_xs_%s[0,-0.5,1.5]" % sample)
            roo.factory(w, "expr::xs_%s('(1+@0)*@1',{d_xs_%s, xs_%s_hat})" % (3 * (sample,)))

    def import_efficiencies(self, w, channels=None):
        if not channels: channels = self.channels
        [roo.wimport_const(w, 'eff_%s_%s' % (lepton, sample), data.eff)
         for lepton, channel in channels.items()
         for sample, data in channel.samples.items()
         if sample != 'data']

    def import_shapes(self, w, channels=None):
        if not channels:
            channels = self.channels
            [roo.factory(w, "%s[-1,1]" % obs) for obs in self.observables]
            roo.factory(w, "channel[%s]" % ','.join("%s=%d" % (s, i)
                                                    for i, s in enumerate(channels)))
            for v, X in zip(self.observables, 'XY'):
                w.var(v).setBins(getattr(self.channels['el'].samples['data'].datas[0],
                                         'GetNbins' + X)())

        [self.import_shape(w, lepton, sample, data)
         for lepton, channel in channels.items()
         for sample, data in channel.samples.items()
         if sample != 'data']

    @roo.quiet
    def import_shape(self, w, lepton, sample, data):
        name = '_'.join([lepton, sample])
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)

        for i, label in enumerate(['both', 'symm'][:None if sample in ['tt', 'dy'] else -1]):
            nL = (name, label)
            roo.wimport(w, r.RooDataHist('_sim_'.join(nL), '', arglist, data.datas[i]))
            roo.wimport(w, r.RooHistPdf('_'.join(nL), '', argset, w.data('_sim_'.join(nL))))

        roo.factory(w, "prod::expect_%s(%s)" %
                    (name, ','.join(['lumi_' + lepton,
                                     'xs_' + sample,
                                     'eff_' + name,
                                     '-1',
                                     'factor_%s' % lepton][:None if 'qcd' in lepton else -2])))

    def import_qcd(self, w):
        [roo.factory(w, "factor_%s[1,0,5]" % lepton) for lepton in self.channels_qcd]

        self.import_efficiencies(w, self.channels_qcd)
        self.import_shapes(w, self.channels_qcd)
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)
        for L, channel in self.channels_qcd.items():
            hist = channel.samples['data'].datas[0]
            dhist = '%s_data_sim_both' % L
            roo.wimport(w, r.RooDataHist(dhist, '', arglist, hist))
            roo.wimport(w, r.RooHistPdf('%s_data_both' % L, '', argset, w.data(dhist)))
            roo.factory(w, "expr::expect_%s_data('@0*%f',{factor_%s})" %
                        (L, hist.Integral(), L))

    def import_asymmetry(self, w):
        alpha_max = 0.99 * min(chan.samples['tt'].alphaMax for chan in
                                self.channels.values() +
                                self.channels_qcd.values())
        roo.factory(w, "alpha[1, -%f, %f]"%(alpha_max,alpha_max))

        [roo.factory(w, "SUM::%(n)s( alpha * %(n)s_both, %(n)s_symm )" % {'n': L+'_tt'})
         for L in (self.channels.keys() + self.channels_qcd.keys())]

        assert self.gen.samples['tt'].datas[0].GetXaxis().GetTitle() == genNameX
        assert self.gen.samples['tt'].datas[0].GetYaxis().GetTitle() == genNameY

        for n, d in self.gen.samples.items():
            hists = {'XL':d.datasX[0], 'XT':d.datasY[0]}
            ac, eac = lib.asymmetry(hists[self.observables[0]])
            roo.wimport_const(w, 'Ac_'+n, ac)
            roo.wimport_const(w, 'err_Ac_'+n, eac)

    def import_model(self, w, whichs={}, name=""):
        for part in ['','qcd']:
            if part not in whichs:
                whichs[part] = dict((i, '_both') for i in ['wj', 'st'])
                whichs[part].update({'tt': ''})
        whichs[''].update({'dy': '_symm'})

        [roo.factory(w, "SUM::%s_mj( expect_%sqcd_data * %sqcd_data_both, %s )" %
                     (name+lepton, lepton, lepton,
                      ','.join(['expect_%s_%s * %s_%s%s' %
                                (lepton + 'qcd', key, lepton + 'qcd', key, value)
                                for key, value in whichs['qcd'].items()])))
         for lepton in self.channels]

        [roo.factory(w, "SUM::%smodel_%s( expect_%sqcd_data * %sqcd_data_both, %s )" %
                     (name,lepton, lepton, lepton,
                      ','.join(['expect_%s_%s * %s_%s%s' %
                                (lepton + part, key, lepton + part, key, value)
                                for part, which in whichs.items()
                                for key, value in which.items()
                                ])))
         for lepton in self.channels]

        roo.factory(w, "SIMUL::%smodel(channel, %s)" %
                    (name, ', '.join("%s=%smodel_%s" %
                                     (lepton, name, lepton) for lepton in self.channels)))

    def import_alt_model(self, channelDict):
        w = self.w

        sample = "ttalt"
        leptons = ['el', 'mu']
        channels = dict((L, channelDict[(L,'top')]) for L in leptons)

        #xs
        xs = channels['el'].samples[sample].xs
        roo.wimport_const(w, 'xs_%s_hat' % sample, xs)
        roo.factory(w, "d_xs_%s[0,-0.5,1.5]" % sample)
        roo.factory(w, "expr::xs_%s('(1+@0)*@1',{d_xs_%s, xs_%s_hat})" % (3 * (sample,)))

        # eff
        [roo.wimport_const(w, 'eff_%s_%s' % (lepton, sample), channel.samples[sample].eff)
         for lepton, channel in channels.items() ]

        [self.import_shape(w, lepton, sample, data)
         for lepton, channel in channels.items()
         for sample, data in channel.samples.items()]

        self.import_model(w, whichs={'':dict((i, '_both')
                                             for i in ['dy', 'wj', 'st', 'ttalt'])},
                          name="alt")

    def import_expressions(self, w):
        [roo.factory(w, "sum::expect_%s_notqcd(%s)" %
                     (L, ','.join(['expect_%s_%s' % (L, s)
                                   for s in ['wj', 'st', 'tt']])))
         for L in self.channels_qcd]

        [roo.factory(w, "sum::expect_%(n)s_mj(expect_%(n)sqcd_data,expect_%(n)sqcd_notqcd)" %
                     {'n': L}) for L in self.channels]

    @roo.quiet
    def import_data(self, altData=None):
        w = self.w
        if altData:
            altData.SetName('data')
            roo.wimport(w, altData)
            self.altData = altData
            for d in altData.split(w.arg('channel')):
                d.SetName('data_'+d.GetName())
                roo.wimport(w, d)
            return
        obs_ = w.argSet(','.join(self.observables))
        obs = r.RooArgList(obs_)

        dataHists = dict([(L,chan.samples['data'].datas[0]) for L,chan in self.channels.items()])
        datas = [(L, r.RooDataHist('data_' + L, 'N_{obs}^{%s}' % L, obs, data))  for L, data in dataHists.items()]
        if not self.quiet:
            print self.observables[0]
            print '\n'.join('%s: %.6f (%.6f)'%((L,)+lib.asymmetry(hist.ProjectionX())) for L,hist in dataHists.items())
            print

        [roo.wimport(w, dat) for _, dat in datas]
        args = [r.RooFit.Import(*dat) for dat in datas]
        roo.wimport(w, r.RooDataHist('data', 'N_{obs}', obs,
                                     r.RooFit.Index(w.arg('channel')), *args))

    def print_n(self,logfile):
        scale = 1000.
        w = self.w
        length = 24

        chans = ['el','mu']
        cross = ['tt', 'wj', 'mj', 'st', 'dy']
        print>>logfile, '\t','&','\t&\t'.join(r'\multicolumn{2}{c}{N_{%s}}'%xs for xs in cross), 'Total', 'Observed'
        for chan in chans:
            tot = 0
            tote2 = 0
            print>>logfile, chan,'&',
            for xs in cross:
                val = w.arg('expect_%s_%s' % (chan, xs)).getVal()
                delta = w.arg('d_xs_%s'%xs)
                factor= w.arg('factor_%sqcd'%chan)
                relerr = max(0.0000001,
                             delta.getError() / (1+delta.getVal()) if delta else
                             factor.getError() / factor.getVal() if xs=='mj' else 0)
                tot += val
                tote2 += (val*relerr)**2
                print>>logfile, lib.roundString(val/scale,relerr*val/scale).rjust(length / 3), '&',
            print>>logfile, lib.roundString(tot/scale,math.sqrt(tote2)/scale).rjust(length / 3), '&',
            print>>logfile, self.channels[chan].samples['data'].datas[0].Integral()/scale, r'\\'
        print>>logfile


    @roo.quiet
    def chi2(self,lep):
        w = self.w
        mod = w.pdf('model_%s'%lep)
        dhist = unqueue(self.channels[lep].samples['data'].datas[0], True)
        exp = mod.generateBinned(w.argSet(','.join(self.observables)), 0, True, False)
        hist = exp.createHistogram(','.join(self.observables), 25, 5)
        chi2 = 5*[0]
        for iX in range(25):
            for iY in range(5):
                    ibin = hist.GetBin(iX+1,iY+1)
                    P = hist.GetBinContent(ibin)
                    Q = dhist.GetBinContent(ibin)
                    chi2[iY] += (P-Q)**2 / P
        del dhist
        del hist
        return sum(chi2)

    def expected_histogram(self, pdfname, expect=0):
        w = self.w
        mod = w.pdf(pdfname)
        args = ','.join(self.observables)
        exp = mod.generateBinned(w.argSet(args), expect, True, False)
        hist = exp.createHistogram(args, 25, 5)
        hist.SetName(pdfname+'genHist')
        return hist

    def data_hist(self, lep):
        w = self.w
        data = w.data('data') if not hasattr(self,'altData') else self.altData
        tmp = next(d for d in data.split(w.arg('channel')) if d.GetName()==lep)
        dhist = tmp.createHistogram('altaData'+'_hist_'+lep, w.arg(self.observables[0]), r.RooFit.Binning(25,-1,1), r.RooFit.YVar(w.arg(self.observables[1]), r.RooFit.Binning(5,-1,1)))
        return dhist


    @staticmethod
    def proj(h):
        return [h.ProjectionX( h.GetName()+'x%d'%i,i+1,i+1) for i in range(5)]

    @roo.quiet
    def visualize(self, printName=''):
        w = self.w
        titles = ['X_{%s}'%self.observables[0],'#Delta']
        for v,t in zip(self.observables,titles) :
            w.arg(v).SetTitle(t)

        r.gROOT.ProcessLine(".L lib/tdrstyle.C")
        r.setTDRStyle()
        r.TGaxis.SetMaxDigits(4)
        canvas = r.TCanvas()
        canvas.Divide(5,2,0,0)
        canvas.Print(printName+'[')

        for j,lep in enumerate(['el','mu']):
            dhist = self.data_hist(lep)
            hist = self.expected_histogram('model_'+lep)

            print self.channels[lep].samples.keys()
            hists = dict([(s, self.expected_histogram(lep+'_'+s+('_both' if s in ['wj','st'] else '_symm' if s=='dy' else ''), 
                                                 w.arg('expect_%s_%s'%(lep,s)).getVal()))
                     for s in self.channels[lep].samples if s not in ['data']])
            hists['mj'] = self.expected_histogram(lep+'_mj', w.arg('expect_%s_mj'%lep).getVal())

            stack = [('tt',),('wj',),('mj',), ('dy','st')][::-1]
            colors = [r.kViolet,r.kGreen+1,r.kRed,r.kGray][::-1]
            stackers = []
            for s,c in zip(stack,colors):
                h = hists[s[0]]
                for n in s[1:]: h.Add(hists[n])
                h.SetFillColor(c)
                h.SetLineColor(c)
                stackers.append(h)


            def proj(h): return [lib.symmAnti(hist) for hist in self.proj(h)]
            model = proj(hist)
            data = proj(dhist)
            comps = [proj(h) for h in stackers]
            hstack = [r.THStack('stack%d'%i,'') for i in range(5)]
            maximum = 1.1 * max(d[0].GetMaximum() for d in data)
            amaximum = 1.7 * max([max(abs(h[1].GetMaximum()),abs(h[1].GetMinimum())) for hist in [data,model] for h in hist])
            antis = []
            for i in range(5):
                canvas.cd(i+1)
                for c in comps: 
                    hstack[i].Add(c[i][0])
                model[i][0].SetMinimum(0)
                model[i][0].SetLineColor(r.kBlue)
                data[i][0].SetMinimum(0)
                data[i][0].SetMaximum(maximum)
                data[i][0].Draw()
                model[i][0].Draw('hist same')
                hstack[i].Draw('hist same')
                data[i][0].Draw('same')

                canvas.cd(i+6)
                data[i][1].SetMinimum(-amaximum)
                data[i][1].SetMaximum(amaximum)
                data[i][1].Draw()
                model[i][1].SetLineWidth(2)
                model[i][1].SetLineColor(r.kBlue)
                model[i][1].Draw('hist same')
                for h,c in zip(comps,colors)[-1:]:
                    h[i][1].ResetAttFill()
                    h[i][1].Draw('hist same')
            sys.stdout.write(' ')
            canvas.Print(printName)

        print
        canvas.Print(printName+']')
        print printName, 'written'

        return


    def Ac_raw(self, channel, model=None):
        hist = (self.data_hist(channel) if not model 
                else self.expected_histogram(model + '_'+channel))

        return [lib.asymmetry(h)[0] for h in self.proj(hist)]

