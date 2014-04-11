

def ensemble_specs():
    base = {'alpha':1.0, 'lumiFactor':1.0, 'Nens':100}

    alphas = [-3,-2,-1.5,-1.25,-1.0,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1.0,1.25,1.5,2,3]
    lumis = [0.25,0.5,2.0,4.0,8.0]

    ups = ([{'label':'A%+.2f'%a, 'alpha':a} for a in alphas] +
           [{'label':'L%.2f'%l, 'lumiFactor':l} for l in lumis])

    specs = [dict(base) for u in ups]
    for s,u in zip(specs,ups): s.update(u)
    return specs


def calibration_specs():
    return [{'which':'mg', 'sample':'calib_mg.pu.sf'},
            {'which':'mn', 'sample':'calib_mn.pu.sf'},
            {'which':'ZP', 'sample':'calib_ZP.pu.sf'},
            {'which':'A2K', 'sample':'calib_A2K.pu.sf'},
            {'which':'R2K', 'sample':'calib_R2K.pu.sf'},
            {'which':'A.2K', 'sample':'calib_A200.pu.sf'},
            {'which':'R.2K', 'sample':'calib_R200.pu.sf'},
            {'which':'L.2K', 'sample':'calib_L200.pu.sf'}]
