# -*- coding: utf-8 -*-
'''
Python 3.10.6
12 September 2022
'''

from datetime import datetime
import os

def update_log(message):
    logfile = 'pfastr.log'
    if not os.path.exists(logfile):
        with open(logfile, 'w') as outf: pass
    output = '{0} {1} \n'.format(message, datetime.now())
    with open(logfile, 'a') as outf:
        outf.write(output)

def print_prediction(prediction):
    print('serotype prediction for ' + prediction[0])
    print(prediction[1])
    print('probability {}'.format(prediction[2]))
    print(prediction[3])

def write_results(outdir, prediction):
    if not os.path.exists(outdir): os.system('mkdir ' + outdir)
    outpath = outdir + '/prediction.txt'
    with open(outpath, 'w+') as pred_out:
        for x in prediction:
            pred_out.write(str(x) + '\n')
