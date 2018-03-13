#!/ust/bin/env python

import sys

with open(sys.argv[1], 'r') as handle:
    num = 0
    for line in handle:
        num += 1
        if num == int(sys.argv[2]):
            print(str(num) + ":\t" + line)
            break

# awk 'NR==YOUR_LINE{print}' file_name
# sed -n -e YOUR_LINEp file_name
# perl -wnl -e '$. == YOUR_LINE and print and exit;'
# less -SN +TOUR_LINEg file_name (nice!)