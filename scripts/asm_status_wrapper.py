#!/usr/bin/env python
import os
import sys
import shutil
file_in = []
with open(sys.argv[1], 'r') as handle:
	for i in handle:
		file_in.append(i.strip())
fx_in = ",".join(file_in)
statswrapper = shutil.which("statswrapper.sh")
print(statswrapper + " in=" + fx_in)
