#!/usr/bin/python

import sys

if len(sys.argv) != 3:
	sys.stderr.write("usage: %s <infile> <each_nth>\n" % sys.argv[0])

f = open(sys.argv[1], 'r')
each = int(sys.argv[2])

i = 0

for l in f:
	i += 1
	if i == each:
		sys.stdout.write(l)
		i = 0
