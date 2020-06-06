# CIMP - C Image Manipulation Program

## Overview

**CIMP** (originally iplib) is a simple CLI image manipulation program, written in standard ANSI C.
Developed as the final assignment for 2020 programming course. This repo includes also the project documentation (italian!).

## Makefile

The provided makefile contains several useful recipes for building the library:

* build for production: `make build`
* build for testing: `make test`
* perform memory tests (requires valgrind installed on your system): `make test-run`
* clean generated files: `make clean`

## Usage

After compiling, run `main_iplib` with no arguments. You'll receive a descriptive prompt on how to use the CLI utility.

## Licensing

The underlining bmp library provided with this project is an adaptation of the [Bitmap API](https://github.com/wernsey/bitmap) developed by [Werner Stoop](https://github.com/wernsey). The Bitmap API is released under the [MIT LICENSE](./LICENSE_BMP).

Moreover, several parts of the code were originally provided as part of the assignment. Head to specific files for the respective copyright notice.

This project is released under the [MIT LICENSE](./LICENSE).

## Authors

* [Ina Popescu](https://github.com/Ina-pps)
* [Matteo Agnoletto](https://github.com/EPMatt)
* [Lorenzo Armando Donatelli](https://github.com/Donnyz)
