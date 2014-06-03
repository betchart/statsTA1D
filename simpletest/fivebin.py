import ROOT as r

def fivebin(N, aout, central=0, ain=0, middle=0):
    assert abs(aout) < 1
    assert abs(ain) < 1
    assert middle <= central < 1, "c: %f; m: %f" % (central, middle)
    assert 0 <= middle < 1

    bins = 5 if middle else 4 if ain else 3 if central else 2

    if not (middle or ain):
        middle = central

    args = (N,aout,central,ain,middle)
    name = "fivebin_N%d_Aout%.2f_C%.2f_Ain%.2f_M%.2f" % args
    title = "N=%d, A_{out}=%.2f, C=%.2f, A_{in}=%.2f, M=%.2f" % args
    h = r.TH1D(name, title, bins, -1, 1)

    Nin = N * (central-middle)
    Nout = N * (1 - central)
    
    h.SetBinContent(1, 0.5*Nout*(1-aout))
    h.SetBinContent(bins, 0.5*Nout*(1+aout))
    if bins>2:
        h.SetBinContent(bins/2 + 1, N * middle)
    if ain:
        h.SetBinContent(bins/2, 0.5*Nin*(1-ain))
        h.SetBinContent(1 + bins - bins/2, 0.5*Nin*(1+ain))

    return h
