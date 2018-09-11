#!/usr/bin/env python
import argparse
parses=argparse.ArgumentParser(description="trim scaffold or random from fa.")
parses.add_argument('-i','--input',dest='input', default=None, help='read in ref fa.')
parses.add_argument('-l','--label',dest='label',default='scaffold', 
                    help='the label which is contained by a chr and thr chr will be trimmed')
parses.add_argument('-o','--out',dest='output', default=None, help='wirteout')

args=parses.parse_args()

readin = args.input
writeto = args.output
label = args.label.upper()

with open(readin) as rf, open(writeto, 'w') as wf:
    writeOut = False
    for line in rf:
        if '>' in line:
            if label not in line.upper():
                writeOut = True
            else:
                writeOut = False
        if writeOut:
            wf.write(line)
            