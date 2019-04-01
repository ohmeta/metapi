#!/usr/bin/env python

import sys
import shutil
import argparse

STATSWRAPPER_TEMPLATE = '''{stats} \
in={input_list} \
minscaf={minscaf} > {output}'''


class statswrapper:
	def __init__(self, input_list, minscaf, output):
		self.stats = shutil.which("statswrapper.sh")
		self.input_list = ",".join(input_list)
		self.minscaf = minscaf
		self.output = output


def gen_shell(ilist, mlen, split, prefix, output):
	files = open(ilist, 'r').readlines()
	total = len(files)
	assert total >= split, "can't split"
	step = total // split
	m = total % split
	count = 0
	sub_files = []
	cmds = []
	for i in range(0, total, step):
		count += 1
		if count <= split:
			sub_files = [f.strip() for f in files[i:(i + step)]]
			output_ = "%s.%d.tsv" % (prefix, count)
			cmd = STATSWRAPPER_TEMPLATE.format_map(
				vars(statswrapper(sub_files, mlen, output_)))
			cmds.append(cmd)

		if (count > split) and (m > 0):
			sub_files += [f.strip() for f in files[(total - m):total]]
			output_ = "%s.%d.tsv" % (prefix, split)
			cmd = STATSWRAPPER_TEMPLATE.format_map(
				vars(statswrapper(sub_files, mlen, output_)))
			cmds[split - 1] = cmd

	with open(output, 'w') as oh:
		for i in cmds:
			oh.write(i + "\n")


def main():
	parser = argparse.ArgumentParser("assembler status wrapper")
	parser.add_argument('-l', '--list', type=str, help='input assembly file list')
	parser.add_argument('-m', '--min_len', type=int, default=0, help='minimal contig/scaffold length')
	parser.add_argument('-s', '--split', type=int, default=1, help='split input file')
	parser.add_argument('-p', '--prefix', type=str, default="asm_stats", help="assembly status output prefix")
	parser.add_argument('-o', '--output', type=str, default=sys.stdout, help='write cmd to file, default: stdout')
	args = parser.parse_args()

	gen_shell(args.list, args.min_len, args.split, args.prefix, args.output)


if __name__ == '__main__':
	main()
