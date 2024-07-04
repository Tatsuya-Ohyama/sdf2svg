# sdf2svg.py

## 概要
Program to convert .sdf file to 2D structure image file


## 使用方法
```sh
$ sdf2svg.py [-h] [-i INPUT.sdf] [-o OUTPUT_PREFIX] [-f FORMAT] -k PROP_NAME [--width WIDTH] [--height HEIGHT] [--keep-3D] [-O] [-l]
```

Option

* `-h`, `--help`
	: show this help message and exit
* `-i INPUT.sdf`
	: source .sdf file
* `-o OUTPUT_PREFIX`
	: output prefix
* `-f FORMAT`
	: output format (starts with period) (".png" and ".svg")
* `-k PROP_NAME`
	: property name of unique name
* `--width WIDTH`
	: image width (Default: 300)
* `--height HEIGHT`
	: image height (Default: 300)
* `--keep-3D`
	: output 3D structure
* `-O`
	: overwrite forcibly
* `-l`
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
### Ver. 1.4 (2024/07/03)
* Add `--keep-3D` option.

### Ver. 1.3 (2024/06/21)
* Rename program name (convert_sdf2svg.py -> sdf2svg.py)

### Ver. 1.2 (2024/06/20)
* Add -k and -l options.

### Ver. 1.1 (2024/06/05)
* Change `-o` to an optoinal attribute.

### Ver. 1.0 (2024/06/04)
* Released.
