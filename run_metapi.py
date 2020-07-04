#!/usr/bin/env python3

import os
import sys
from site import addsitedir

from metapi import corer

addsitedir(os.path.dirname(os.path.dirname(corer.__file__)))
print(sys.path)

if __name__ == '__main__':
    print(corer.__file__)
    corer.main()
