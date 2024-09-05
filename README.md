# sdf2svg.py

## 概要
Program to convert .sdf file to 2D structure image file


## 使用方法
```sh
$ sdf2svg.py [-h] [-i INPUT.sdf] [-o OUTPUT_PREFIX] [-f FORMAT] [-p PROP_NAME] [-s SIZE] [--keep-3D] [--label] [--grid MxN] [-O] [-l]
```

Option

* `-h`, `--help`
	: show this help message and exit
* `-i INPUT.sdf`
	: source .sdf file
* `-o OUTPUT_PREFIX`
	: output image file
* `-f FORMAT`, `--format FORMAT`
	: output format (starts with period)
* `-p PROP_NAME`, `--property PROP_NAME`
	: property name of unique name (for using output filepath and label)
* `-s SIZE`, `--size SIZE`
	: image size (Default: 300)
* `--keep-3D`
	: output 3D structure
* `--label`
	: add label (Default: False)
* `--grid MxN`
	: draw molecules in grid by M x N (Default: None)
* `-O`
	: overwrite forcibly
* `-l`, `--list`
	: show PROP_NAME list and exit



## 動作要件
* Python3
	* rdkit


## License
The MIT License (MIT)

Copyright (c) 2024 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama

## ChangeLog
### Ver. 1.5 (2024/09/05)
* Change `-p` from a required option to an optional one.
* Suppress output file name duplication.
* Support large ring structure.

### Ver. 1.4 (2024/07/22)
* Change `-k` option to `-p` option.
* Remove `--width` and `--height` options and add `-s` option.
* Add `--keep-3D` option.
* Add `--label` option.
* Add `--grid` option.

### Ver. 1.3 (2024/06/21)
* Rename program name (convert_sdf2svg.py -> sdf2svg.py)

### Ver. 1.2 (2024/06/20)
* Add -k and -l options.

### Ver. 1.1 (2024/06/05)
* Change `-o` to an optoinal attribute.

### Ver. 1.0 (2024/06/04)
* Released.
