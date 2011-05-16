#!/usr/bin/python
## @file
## Do a test run and ask for parameters

import re
import sys
import subprocess
import time

## parameters
params = {}

# Read parameter file
f = open('params.txt', 'r')

for line in f:
	if re.match(r'^#', line):
		continue # read over comment
	
	m = re.match(r'(\w+)\s+(\S+)\s*.*', line)

	if m:
		name, default_value = m.group(1,2)
		value = raw_input("%s [%s]: " % (name, default_value))
		if len(value.strip()) == 0:
			value = default_value
		params[name] = value
	elif len(line.strip()) == 0:
		pass # ignore empty line
	else:
		sys.stderr.write("Warning: Line not recognized: \"%s\"" % line)

print ''
print "======================================"
print " Starting simulation..."
print "======================================\n"

p = subprocess.Popen(('./cylindric'), stdin=subprocess.PIPE)

for n, v in params.items():
	p.stdin.write("%s %s\n" % (n, v))

try:
	p.communicate()
except KeyboardInterrupt:
	print ''
	print "======================================"
	print " Aborted! "
	print "======================================\n"
	print "Trying to terminate..."

	try:
		p.terminate()
	except:
		print "  -> Failed..."
	
	time.sleep(3)

	print "Trying to kill..."

	try:
		p.terminate()
	except:
		print "  -> Failed..."
