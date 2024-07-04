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
	parser.add_argument("-f", dest="OUTPUT_FORMAT", metavar="FORMAT", default=".png", help="output format (starts with period)")
	parser.add_argument("-k", dest="PROP_NAME", metavar="PROP_NAME", required="-l" not in sys.argv, help="property name of unique name")
	parser.add_argument("--width", dest="OUTPUT_WIDTH", metavar="WIDTH", type=int, default=300, help="image width (Default: 300)")
	parser.add_argument("--height", dest="OUTPUT_HEIGHT", metavar="HEIGHT", type=int, default=300, help="image height (Default: 300)")
	parser.add_argument("--keep-3D", dest="FLAG_KEEP_3D", action="store_true", default=False, help="output 3D structure")
	parser.add_argument("--label", dest="FLAG_ADD_LABEL", action="store_true", default=False, help="add label (Default: False)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("-l", dest="FLAG_SHOW_PROP_NAMES", action="store_true", default=False, help="show PROP_NAME list and exit")
	args = parser.parse_args()

	# check input file
	if not os.path.isfile(args.INPUT_SDF):
		sys.stderr.write("ERROR: no such file `{}`.\n".format(args.INPUT_SDF))
		sys.exit(1)

	if args.FLAG_SHOW_PROP_NAMES:
		list_prop_names = set([])
		for obj_mol in Chem.SDMolSupplier(args.INPUT_SDF):
			if obj_mol is None:
				continue

			list_prop_names |= set(obj_mol.GetPropNames())

		for prop_name in list_prop_names:
			print(prop_name)
		sys.exit(0)

	# check output position
	output_dir = os.path.dirname(args.OUTPUT_PREFIX)
	if len(output_dir) != 0 and not os.path.isdir(output_dir):
		os.makedirs(output_dir)

	# output image
	for mol_i, obj_mol in enumerate(Chem.SDMolSupplier(args.INPUT_SDF), 1):
		if obj_mol is None:
			continue

		if args.PROP_NAME not in obj_mol.GetPropNames():
			sys.stderr.write("WARNING: Not found PROP_NAME `{}` at {} th molecule. Skipped...\n".format(args.PROP_NAME, mol_i))
			continue

		mol_name = obj_mol.GetProp(args.PROP_NAME)
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
		if not args.FLAG_KEEP_3D:
			AllChem.Compute2DCoords(obj_mol)

		# output image
		Draw.MolToFile(obj_mol, output_file, size=(args.OUTPUT_WIDTH, args.OUTPUT_HEIGHT), legend=obj_mol.GetProp(args.PROP_NAME))
		print("output:", output_file)
