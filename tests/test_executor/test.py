#!/usr/bin/env python


from executor import execute

output = str(snakemake.output)
execute(f'''touch {output}''')