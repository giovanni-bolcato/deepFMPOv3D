# -*- coding: utf-8 -*-
#

import argparse
import datetime

now = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

parser = argparse.ArgumentParser(description='fragment molecules and builds improved analogs')
parser.add_argument('-fragmentFile', dest='fragment_file', metavar='-f', type=str, required=True, help='molecule file (.smi) to be fragmented')
parser.add_argument('-leadSeriesFile', dest='lead_file', metavar='-l', type=str, required=True, help='lead molecules to be optimized, supplied as .smi')
parser.add_argument('-outputFile', dest='out_file', metavar='-o', type=str, required=True, default="myMPOcompounds_{}.csv".format(str(now)), help='results end up here')


def read_args():
    args = parser.parse_args()

    fragment_file = args.fragment_file
    lead_file = args.lead_file
    out_file = args.out_file

    print (args)

    return args
