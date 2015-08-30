from icqsol.color.icqColorMap import ColorMap

fmin = 0.0
fmax = 1.0
cm = ColorMap(fmin=fmin, fmax=fmax)

n = 10
for i in range(n):
	f = fmin + (fmax - fmin)*i/float(n)
	print 'hot       f = {}: {}'.format(f, cm.hot(f))
	print 'cold      f = {}: {}'.format(f, cm.cold(f))
	print 'blackbody f = {}: {}'.format(f, cm.blackbody(f))
