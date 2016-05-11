#!/usr/bin/env python

actions = {
	'insert before': [
		(
		'^import',
		'from __future__ import print_function',
		),
	],

	'replace': [
		(
		 'print\s+([^\(].*)\s*\,\s*\n\s*$',
	     'print(\\1, end="")\n',
	    ),
		(
		 'print\s+([^\(].*[^\,])\s*\n\s*$',
	     'print(\\1, end="\\\\n")\n',
	    ),
	],
}
