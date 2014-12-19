partitions = ['full', 'hiM', 'loM', 'hiY', 'loY']
from asymmNames import genNames

def measurement_pars(partition='full', var='XL'):
    fields = ('R0_', 'signal')
    base = 3
    hemicycle = 5
    asymmetry = (base, 
                 genNames[var].replace('gen','fit')+'_TridiscriminantWTopQCD')

    pars = dict(zip(fields, asymmetry))
    N = pars['R0_']
    cycle = 2*hemicycle
    pDirs = [N,
             N + 1*cycle, (N, N + 1*cycle),
             N + 2*cycle, (N, N + 2*cycle), (N, N + 2*cycle), (N, N + 1*cycle)]
    pars.update({'R0_': dict(zip(partitions,pDirs))[partition]})
    pars.update({'label':'_'.join([var,partition])})
    return pars


def central():
    return {'d_lumi': 0,
            'd_xs_dy': 0,
            'd_xs_st': 0,
            'tag': 'ph_sn_jn_20',
            'genDirPre': 'R01',
            'genPre': '',
            'sigPre': '',
            'dirIncrement': 0,
            'd_wbb':0,
            'label': 'central',
            'twossigma': {}
            }


def systematics():
    sys =  ([
            {'label': 'BTAG', 'partsuffix': 'SF'},

            {'label': 'JER_up', 'tag': 'ph_su_jn_20'},
            {'label': 'JER_dn', 'tag': 'ph_sd_jn_20'},

            {'label': 'JES_up', 'tag': 'ph_sn_ju_20'},
            {'label': 'JES_dn', 'tag': 'ph_sn_jd_20'},

            {'label': 'PU_up', 'dirIncrement': 1, 'sigPre': '001_'},
            {'label': 'PU_dn', 'dirIncrement': 1, 'sigPre': '000_'},

            {'label': 'PT', 'dirIncrement': 4, 'sigPre': '001_', 'genPre': '001_', 'genDirPre':'R02'},

            {'label': 'WBB_up', 'd_wbb':+0.2},
            {'label': 'WBB_dn', 'd_wbb':-0.2},

            {'label': 'lumi_up', "d_lumi": +0.044},
            {'label': 'lumi_dn', "d_lumi": -0.044},
            
            {'label': 'DY_up', "d_xs_dy": +0.20},
            {'label': 'DY_dn', "d_xs_dy": -0.20},
            
            {'label': 'ST_up', "d_xs_st": +0.20},
            {'label': 'ST_dn', "d_xs_st": -0.20},

            {'label': 'RT_up', "Rst":2.0},  # from https://cds.cern.ch/record/1528574
            {'label': 'RT_dn', "Rst":1.5},  # Rt-ch. = sigma(t-ch., top) / sigma(t-ch., anti-top) = 1.76 +/- 0.27

            {'label': 'WJ_up', 'twossigma': {'d_xs_wj': +1.0}},
            {'label': 'WJ_dn', 'twossigma': {'d_xs_wj': -1.0}},

            {'label': 'TT_up', 'twossigma': {'d_xs_tt': +1.0}},
            {'label': 'TT_dn', 'twossigma': {'d_xs_tt': -1.0}},

            {'label': 'QCDe_up', 'twossigma': {'factor_elqcd': +1.0}},
            {'label': 'QCDe_dn', 'twossigma': {'factor_elqcd': -1.0}},

            {'label': 'QCDm_up', 'twossigma': {'factor_muqcd': +1.0}},
            {'label': 'QCDm_dn', 'twossigma': {'factor_muqcd': -1.0}},

            {'label': 'Q_dn', 'genPre': '053_', 'sigPre': '053_'},
            {'label': 'Q_up', 'genPre': '054_', 'sigPre': '054_'},

            {'label': 'as_dn', 'genPre': '055_', 'sigPre': '055_'},
            {'label': 'as_up', 'genPre': '056_', 'sigPre': '056_'},

            #{'label': 'thr30', 'tag': 'ph_sn_jn_30'}
            ] +
            [{'label': 'el%d' % i,
              'dirIncrement': 2,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': 'mu%d' % i,
              'dirIncrement': 3,
              'sigPre': '%03d_' % i,} for i in range(4)
            ]+
            [{'label': 'PD_%02d' % i,
              'genPre': '%03d_' % i,
              'sigPre': '%03d_' % i} for i in range(1, 53)
            ]
        )
    return sys


if __name__ == '__main__':
    sys = systematics()
    for s in sys:
        cen = central()
        cen.update(s)
        print cen
    print
    print central()
