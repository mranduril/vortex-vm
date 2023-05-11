#!/usr/bin/env python3
# coding=utf-8
from __future__ import print_function

import sys
import os
import os.path as path
import re
import argparse
from datetime import datetime

script_dir = path.dirname(path.realpath(__file__))

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--input", default='none', help='Verilog header')
parser.add_argument('-o', "--output", default='none', help='C header')

args = parser.parse_args()

if args.input == 'none' or args.output == 'none':
    print('Error: invalid arguments')
    sys.exit()

translation_rules = [
    # preprocessor directives
    (re.compile(r'`include\s+.*$'), r''),
    (re.compile(r'`ifdef'), r'#ifdef'),
    (re.compile(r'`ifndef'), r'#ifndef'),
    (re.compile(r'`else'), r'#else'),
    (re.compile(r'`define'), r'#define'),    
    (re.compile(r'`endif'), r'#endif'),

    # macro expansion
    (re.compile(r"`([A-Za-z_][$_0-9A-Za-z]*)"), r'\1'),

    # literals
    (re.compile(r"\d+'d(\d+)"), r'\1'),
    (re.compile(r"\d+'b([01]+)"), r'0b\1'),
    (re.compile(r"128'h([\da-fA-F_]+)"), r'"\1"'),
    (re.compile(r"\d+'h([\da-fA-F]+)"), r'0x\1')    
]

with open(args.output, 'w') as f:
    print('''
// auto-generated by gen_config.py. DO NOT EDIT
// Generated at {date}

// Translated from {input}:
'''[1:].format(date=datetime.now(), input=args.input), file=f)
    with open(args.input, 'r') as r:
        lineno = 0
        for line in r:
            for pat, repl in translation_rules:
                match = pat.search(line)
                if match:                        
                    line = re.sub(pat, repl, line)
                    #print("*** match @" + str(lineno) + ": " + match.group() + " => " + line)
            f.write(line)
            lineno = lineno + 1
    print('''
'''[1:], file=f)



