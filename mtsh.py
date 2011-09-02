#!/usr/bin/env python
import os,sys,time,multiprocessing,subprocess

'''
DESCRIPTION
	basic multithreaded queue processing script
	pipe list of commands and number of processes (-j) into this script
	(not sophisticated enough to be gracefully interrupted)
'''

def ShellWorker(cmd):
	print cmd
	proc = subprocess.Popen(["bash","-c",cmd])
	while proc.poll() == None:
		time.sleep(0.1)
	return proc.poll()

class MultiShell:
	def __init__(self,commands,nproc=1):
		self.pool = multiprocessing.Pool(nproc)
		results = self.pool.map(ShellWorker,commands)

if __name__ == "__main__":
	from optparse import OptionParser
	p=OptionParser()
	p.add_option('-j','--nprocs',default=4,type='int')
	opt,args=p.parse_args()

	commands=[]
	# commands from pipe
	if len(args) < 1:
		for l in sys.stdin:
			commands.append(l.strip())
		if len(commands) > 0:
			print 'commands read from stdin'
		else:
			print 'no commands received'; sys.exit()
	# or commands from a file
	elif len(args) == 1 and os.path.exists(args[0]):
		print 'reading commands from file %s' %args[0]
		for l in open(args[0]):
			commands.append(l.strip())
	# or commands from argument list
	else:
		print 'reading commands from argument list'
		commands=args

	m=MultiShell(commands,opt.nprocs)
