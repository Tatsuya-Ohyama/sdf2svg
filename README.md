# convert_sdf2svg.py

## 概要
Program to convert .sdf file to 2D structure image file


## 使用方法
```sh
$ convert_sdf2svg.py [-h] [-i INPUT.sdf] [-o OUTPUT_PREFIX] [-f FORMAT] [--width WIDTH] [--height HEIGHT] [-O]
```

* Option
* `-h`, `--help`
	: show this help message and exit
* `-i INPUT.sdf`
	: source .sdf file
* `-o OUTPUT_PREFIX`
	: output .svg file
* `-f FORMAT`
	: output format (".png" and ".svg")
* `--width WIDTH`
	: image width (Default: 300)
* `--height HEIGHT`
	: image height (Default: 300)
* `-O`
	: overwrite forcibly


## 動作要件
* Python3
	* rdkit


## License
The MIT License (MIT)

Copyright (c) 2024 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 1.1 (2024/06/05)
* Change `-o` to an optoinal attribute.

### Ver. 1.0 (2024/06/04)
* Released.
