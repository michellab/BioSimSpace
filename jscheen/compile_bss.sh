#!/bin/bash

current_dir=$(pwd)
cd /home/jscheen/BioSimSpace/python

~/sire.app/bin/ipython setup.py install
cd $current_dir