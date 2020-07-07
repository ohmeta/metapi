#!/usr/bin/env python3

import os
import sys
from site import addsitedir

sys.path.append(os.path.dirname(__file__))

from metapi.corer import main
print(sys.path)
main()
