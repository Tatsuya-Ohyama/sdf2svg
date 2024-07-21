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
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import re

RE_GRID = re.compile(r"^(\d+)x(\d+)$")



# =============== function =============== #
def output_grid(output_file, list_obj_mol, list_label=None):
	useSVG = False
	ftype = "wb"
	if os.path.splitext(output_file)[1] == ".svg":
		useSVG = True
		ftype = "w"

	img = Draw.MolsToGridImage(
		list_obj_mol,
		molsPerRow=n_col,
		subImgSize=(args.SIZE, args.SIZE),
		legends=list_label,
		useSVG=useSVG
	)

	with open(output_file, ftype) as obj_output:
		obj_output.write(img)



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Program to convert .sdf file to 2D structure image file", formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument("-i", dest="INPUT_SDF", metavar="INPUT.sdf", help="source .sdf file")
	parser.add_argument("-o", dest="OUTPUT_PREFIX", metavar="OUTPUT_PREFIX", default="", help="output image file")
	parser.add_argument("-f", "--format", dest="OUTPUT_FORMAT", metavar="FORMAT", default=".png", choices=[".png", ".svg"], help="output format (starts with period)")
	parser.add_argument("-p", "--property", dest="PROP_NAME", metavar="PROP_NAME", required="-l" not in sys.argv, help="property name of unique name (for using output filepath and label)")
	parser.add_argument("-s", "--size", dest="SIZE", metavar="SIZE", type=int, default=300, help="image size (Default: 300)")
	parser.add_argument("--keep-3D", dest="FLAG_KEEP_3D", action="store_true", default=False, help="output 3D structure")
	parser.add_argument("--label", dest="FLAG_ADD_LABEL", action="store_true", default=False, help="add label (Default: False)")
	parser.add_argument("--grid", dest="GRID_OPTION", metavar="MxN", default=None, help="draw molecules in grid by M x N (Default: None)")
	parser.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser.add_argument("-l", "--list", dest="FLAG_SHOW_PROP_NAMES", action="store_true", default=False, help="show PROP_NAME list and exit")
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
	if args.GRID_OPTION is None:
		# for each structure
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
			legend = None
			if args.FLAG_ADD_LABEL:
				legend = obj_mol.GetProp(args.PROP_NAME)

			Draw.MolToFile(obj_mol, output_file, size=(args.SIZE, args.SIZE), legend=legend)

			print("output:", output_file)

	else:
		# for grid
		obj_match = RE_GRID.search(args.GRID_OPTION)
		if not obj_match:
			sys.stderr.write("ERROR: invalid grid option.\n")
			sys.exit(1)

		n_row, n_col = int(obj_match.group(1)), int(obj_match.group(2))
		grid_max = n_row * n_col
		list_obj_mol_grid = []
		n_output = 0
		for mol_i, obj_mol in enumerate(Chem.SDMolSupplier(args.INPUT_SDF), 1):
			if obj_mol is None:
				continue

			if args.PROP_NAME not in obj_mol.GetPropNames():
				sys.stderr.write("WARNING: Not found PROP_NAME `{}` at {} th molecule. Skipped...\n".format(args.PROP_NAME, mol_i))
				continue

			# convert 3D to 2D
			if not args.FLAG_KEEP_3D:
				AllChem.Compute2DCoords(obj_mol)

			if len(list_obj_mol_grid) < grid_max:
				list_obj_mol_grid.append(obj_mol)

			elif len(list_obj_mol_grid) == grid_max:
				n_output += 1
				output_file = "{}{}-{}{}".format(args.OUTPUT_PREFIX, args.GRID_OPTION, n_output, args.OUTPUT_FORMAT)

				if args.FLAG_OVERWRITE == False:
					if os.path.isfile(output_file):
						sys.stderr.write("WARNING: file `{}` exists. Do you want to overwrite it? (y/N) ".format(output_file))
						sys.stderr.flush()

						user_choice = sys.stdin.readline().strip()
						if user_choice not in ["y", "Y"]:
							sys.stderr.write("INFO: skipped `{}`\n".format(output_file))
							list_obj_mol_grid = []
							continue

				if args.FLAG_ADD_LABEL:
					legends = [v.GetProp(args.PROP_NAME) for v in list_obj_mol_grid]

				output_grid(output_file, list_obj_mol_grid, list_label=legends)
				print("output:", output_file)
				list_obj_mol_grid = []

		if len(list_obj_mol_grid) != 0:
			n_output += 1
			output_file = "{}{}-{}{}".format(args.OUTPUT_PREFIX, args.GRID_OPTION, n_output, args.OUTPUT_FORMAT)

			if args.FLAG_ADD_LABEL:
				legends = [v.GetProp(args.PROP_NAME) for v in list_obj_mol_grid]

			output_grid(output_file, list_obj_mol_grid, list_label=legends)
			print("output:", output_file)
