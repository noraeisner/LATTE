
#imports
import os
import sys
from argparse import ArgumentParser
from glob import glob

if __name__ == '__main__':
	ap = ArgumentParser(description='Lightcurve Analysis Tool for Transiting Exoplanets')
	ap.add_argument('--nickname', type=str, help='give the target a memorable name', default='no')

	args = ap.parse_args()

	indir = "./LATTE_output"
	tic = 9966358

	existing_files = glob("{}/*{}*".format(indir, tic))

	if (len(existing_files) > 0): 
		print ("This file already exists therefore SKIP. To overwrite files run this code with --o in the command line.")

	#if not args.nickname == 'no':
	#	os.system("mv {}/{} {}/{}_{}".format(indir, tic, indir, tic, args.nickname))

# End.