from pymol import cmd

def color_destab(posfile,binding=False):
	#re_pos = re.compile('([A-Z])\s+([0-9]+)\s+[a-zA-Z]+\s+([0-9]+)\s+([0-9]+)')
	re_pos = re.compile('([A-Z])\s+([0-9]+)\s+([0-9]+)')
	cols = [ ('destab%i' %i, [float(i)/20,0,0]) for i in range(20) ]
	print cols
	for name,col in cols: cmd.set_color(name,col)
	cmd.color('gray','all')
	for match in re_pos.finditer( open(posfile).read() ):
		#chain,pos,ndestab,nbind = match.groups()
		chain,pos,ndestab = match.groups()
		which_n = ndestab
		if binding: which_n = nbind
		print chain, pos, which_n
		print cols[int(which_n)]
		cmd.color( cols[int(which_n)][0], 'c. %s and i. %s' %(chain,pos) )

cmd.extend('color_destab',color_destab)
