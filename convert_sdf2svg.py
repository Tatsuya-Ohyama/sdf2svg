#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Program to convert .sdf file to 2D structure image file
"""

import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem, Draw



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Program to convert .sdf file to 2D structure image file", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest="INPUT_SDF", metavar="INPUT.sdf", help="source .sdf file")
	parser.add_argument("-o", dest="OUTPUT_PREFIX", metavar="OUTPUT_PREFIX", default="", help="output image file")
	parser.add_argument("-f", dest="OUTPUT_FORMAT", metavar="FORMAT", default=".png", help="output format")
	parser.add_argument("--width", dest="OUTPUT_WIDTH", metavar="WIDTH", type=int, default=300, help="image width (Default: 300)")
	parser.add_argument("--height", dest="OUTPUT_HEIGHT", metavar="HEIGHT", type=int, default=300, help="image height (Default: 300)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	args = parser.parse_args()

	# check input file
	if not os.path.isfile(args.INPUT_SDF):
		sys.stderr.write("ERROR: no such file `{}`.\n".format(args.INPUT_SDF))
		sys.exit(1)

	# check output position
	output_dir = os.path.dirname(args.OUTPUT_PREFIX)
	if len(output_dir) != 0 and not os.path.isdir(output_dir):
		os.makedirs(output_dir)

	# output image
	for obj_mol in Chem.SDMolSupplier(args.INPUT_SDF):
		if obj_mol is None:
			continue

		mol_name = obj_mol.GetProp('_Name')
		output_file = "{}{}{}".format(args.OUTPUT_PREFIX, mol_name, args.OUTPUT_FORMAT)

		if args.FLAG_OVERWRITE == False:
			if os.path.isfile(output_file):
				sys.stderr.write("WARNING: file `{}` exists. Do you want to overwrite it? (y/N) ".format(output_file))
				sys.stderr.flush()

				user_choice = sys.stdin.readline().strip()
				if user_choice not in ["y", "Y"]:
					sys.stderr.write("INFO: skipped `{}`\n".format(mol_name))
					continue

		# convert 3D to 2D
		AllChem.Compute2DCoords(obj_mol)

		# output image
		Draw.MolToFile(obj_mol, output_file, size=(args.OUTPUT_WIDTH, args.OUTPUT_HEIGHT))
		print("output:", output_file)
