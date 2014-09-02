#!/usr/bin/env python

import sys
import lib

pat = sys.argv[1]
end = '' if len(sys.argv)<3 else sys.argv[2]

nums = [int(f.replace(pat,'').replace(end+'.root','')) for f in lib.getCommandOutput("ls %s*.root"%pat)["stdout"].split()]
print max(nums), ': ', ','.join([str(s) for s in sorted(set(range(max(nums))) - set(nums))])
