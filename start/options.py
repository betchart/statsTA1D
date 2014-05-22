from optparse import OptionParser
parser = OptionParser("usage: %prog [options]")
def argOrTrue(option, opt, value, parser) :
    peek = next(iter(parser.rargs),None)
    if peek and peek[0]!='-' : del parser.rargs[0]
    setattr(parser.values, option.dest, peek if peek and peek[0]!='-' else True)

parser.add_option("--XL", dest="XL", default=False, action='store_true', help='measurement of asymmetry in XL')
parser.add_option("--XT", dest="XT", default=False, action='store_true', help='measurement of asymmetry in XT')
parser.add_option("--onlyel", dest='onlyel', default=False, action='store_true', help='measurement only with electron+jets channel')
parser.add_option("--onlymu", dest='onlymu', default=False, action='store_true', help='measurement only with muon+jets channel')
parser.add_option("--partitions", dest="partitions", default=None, metavar='p1,p2,...', action="callback", callback=argOrTrue, help='specify list of partitions, or list partitions')
parser.add_option("--systematics", dest="systematics", default=None, metavar='sys1,sys2,...', action="callback", callback=argOrTrue, help='specify list of systematics, or list systematics')
parser.add_option("--ensembles", dest="ensembles", default=None, metavar='ens1,ens2,...', action="callback", callback=argOrTrue, help='specify list of ensembles or list ensembles')
parser.add_option("--calibrations", dest="calibrations", default=None, metavar='cal1,cal2,...', action="callback", callback=argOrTrue, help='specify list of calibrations or list calibrations')
parser.add_option("--ensSlice", dest="ensSlice", default=None, metavar='lo:hi', help='ensembles by slice notation')
parser.add_option("--calSlice", dest="calSlice", default=None, metavar='lo:hi', help='calibrations by slice notation')
parser.add_option("--templates", dest="templates", default=None, metavar='lo:hi', help='templates by slice notation')
parser.add_option("--visualize", dest='visualize', default=False, action='store_true', help='project the fits')
parser.add_option("--batch", dest='batch', default=False, action='store_true', help='run on the batch queue')
parser.add_option("--site", dest='site', default=None, metavar='ic', help='batch site')
parser.add_option("--chunk", dest='chunk', default=10, metavar='N', type="int", help='number of jobs in a batch chunk')
parser.add_option("--nobg", dest='nobg', default=False, action='store_true',  help='do not include background samples in model')
parser.add_option("--rebin", dest='rebin', default=False, action='store_true',  help='rebin the asymmetry observable')
parser.add_option("--no3D", dest='no3D', default=False, action='store_true', help='rebin the tridiscriminant down to 1 bin')

def opts() :
    options,args = parser.parse_args()
    if options.partitions==None:
        parser.print_help()
        exit()
    if not options.XL^options.XT:
        parser.print_help()
        exit()
    if options.onlyel and options.onlymu:
        parser.print_help()
        exit()
    return options

def default(options = []) :
    options,args = parser.parse_args(options)
    return options
