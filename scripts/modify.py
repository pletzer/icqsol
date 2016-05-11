#!/usr/bin/env python

import argparse
import re
import importlib
import sys
import fileinput

parser = argparse.ArgumentParser(description='Modify file.')
parser.add_argument('--action_file', dest='action_file', default='',
                    help='file containing a list of actions to perform')

args, unk = parser.parse_known_args()
module_without_suffix = args.action_file.split('.')[0]
do = importlib.import_module(module_without_suffix)


def modify():

	insertBeforeFlg = True
	if not do.actions.has_key('insert before'):
		insertBeforeFlg = False

	for line in fileinput.input(unk):
		for oldpat, newpat in do.actions['replace']:
			m = re.search(oldpat, line)
			if m:
				line = re.sub(oldpat, newpat, line)
				break

		if insertBeforeFlg:
			for p in do.actions['insert before']:
				m = re.match(p[0], line)
				if insertBeforeFlg and m:
					line = p[1] + '\n' + line
		   			insertBeforeFlg = False

		print line,

if __name__ == '__main__': 
	modify()
