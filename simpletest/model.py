import lib
from lib import roo
import ROOT as r


class simpleModel(object):
    @roo.quiet
    def __init__(self, template):

        symm,anti = lib.symmAnti(template)
        for item in ['template','symm','anti']:
            setattr(self, item, eval(item))
        
        self.nBinsX = symm.GetNbinsX()

        w = r.RooWorkspace('Workspace')
        roo.factory(w, "x[-1,1]")
        argset = w.argSet("x")
        arglist = r.RooArgList(argset)
        roo.wimport(w, r.RooDataHist('dh_template_both','', arglist, template))
        roo.wimport(w, r.RooDataHist('dh_template_symm','', arglist, symm))
        roo.wimport(w, r.RooHistPdf('pdf_template_both','',argset, w.data('dh_template_both')))
        roo.wimport(w, r.RooHistPdf('pdf_template_symm','',argset, w.data('dh_template_symm')))
        self.import_asymmetry(w)
        self.w = w

    def import_asymmetry(self, w):
        alpha_max = 0.99 * lib.alphaMax(self.symm, self.anti)
        roo.factory(w, "alpha[0, -%f, %f]"%(alpha_max,alpha_max))
        roo.factory(w, 'SUM::model( alpha * pdf_template_both, pdf_template_symm)')
        ac, eac = lib.asymmetry(self.template)
        roo.wimport_const(w, "Ac_template", ac)
        roo.factory(w, "prod::Ac(Ac_template, alpha)")
        

    @roo.quiet
    def import_data(self, data, name=""):
        roo.wimport_const(self.w,"Ac_data"+name, lib.asymmetry(data)[0])
        roo.wimport(self.w,
                    r.RooDataHist('data'+name,
                                  'N_{obs}',
                                  r.RooArgList(self.w.argSet('x')),
                                  data))


    @roo.quiet
    def minimize(self, name=""):
        w = self.w
        nll = w.pdf('model').createNLL(w.data('data'+name))
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        minu.setStrategy(2)
        minu.migrad()


    def val(self,name):
        return self.w.arg(name).getVal()
