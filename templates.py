from inputs import channel_data
import utils
import ROOT as r
import math
r.gROOT.SetBatch(1)
r.gROOT.ProcessLine(".L tdrstyle.C")
r.setTDRStyle()
r.tdrStyle.SetErrorX(r.TStyle().GetErrorX())
r.tdrStyle.SetPadTopMargin(0.065)
r.TGaxis.SetMaxDigits(3)
r.tdrStyle.SetEndErrorSize(6)
#r.tdrStyle.SetPadRightMargin(0.06)

def unqueue(h):
    return utils.unQueuedBins(h,5,[-1,1],[-1,1])


channels = dict([(lep, channel_data(lep, 'top', signal='fitTopQueuedBin5_TridiscriminantWTopQCD')) for lep in ['el', 'mu']])
linetypes = [1, 2]

comps = ['ttgg','ttag','ttqg','ttqq']
colors = [r.kBlack, r.kGreen, r.kBlue, r.kRed]

projections = {}

for lt,(lep,ch) in zip(linetypes,channels.items()):
    for comp,color in zip(comps,colors): 
        symm = unqueue(ch.samples[comp].datas[1])
        anti = unqueue(ch.samples[comp].datas[2])
        [h.SetLineColor(color) for h in [symm,anti]]
        [h.SetLineStyle(lt) for h in [symm,anti]]
        [h.SetMarkerColor(color) for h in [symm,anti]]
        if lep=='mu': [h.SetMarkerStyle(4) for h in [symm,anti]]
        
        sf = 1./symm.Integral()
        symm.Scale(sf)
        anti.Scale(sf)

        projections[(lep,comp)] = [(symm.ProjectionX('symmx'+lep+comp), anti.ProjectionX('antix'+lep+comp)),
                                   (symm.ProjectionY('symmy'+lep+comp), anti.ProjectionY('antiy'+lep+comp)),
                                   (symm.ProjectionZ('symmz'+lep+comp), anti.ProjectionZ('antiz'+lep+comp))
                               ]

def extrema(A,B, func):
    if type(A) not in [list,tuple]:
        return func(A,B)
    return type(A)([extrema(a,b, func) for a,b in zip(A,B)])

def MAX(A,B) : return extrema(A,B, max)
def MIN(A,B) : return extrema(A,B, min)

maxs = reduce( MAX, [[tuple([i.GetMaximum()+i.GetBinError(4) for i in pair]) for pair in L] for L in projections.values()], 3*[(0,0)])
mins = reduce( MIN, [[tuple([i.GetMinimum()-i.GetBinError(4) for i in pair]) for pair in L] for L in projections.values()], 3*[(0,0)])

c = r.TCanvas()
fn = 'graphics/template.pdf'
c.Print(fn+'[')

for i,label in enumerate(['X_{L}','X_{T}']):
    for j,sublabel in enumerate(['symmetrized ','antisymmetrized ']):
        init = False
        for k,comp in enumerate(reversed(comps)):
            for lep in channels:
                h = projections[(lep,comp)][i][j]
                if j: 
                    h.SetBinError(3,0)
                    h.SetMinimum(-0.0065)
                    h.SetMaximum(0.0065)
                else:
                    h.SetMinimum(0)
                    h.SetMaximum(0.265)
                h.GetXaxis().SetNdivisions(5,4,0,False)
                #h.SetMaximum(1.1*maxs[i][j])
                #h.SetMinimum(1.1*mins[i][j])
                h.SetMarkerSize(1.5 - 0.3*math.sqrt(k) - (0 if lep=='mu' else 0.3))
                h.SetLineWidth(4-k)
                h.GetXaxis().SetTitleOffset(0.9)
                h.GetYaxis().SetTitle(sublabel + ' probability')
                h.GetXaxis().SetTitle(label)
                h.Draw('same e1' if init else 'e1')
                init = True
        c.Print(fn)


c.Print(fn+']')

