# -*- coding: utf-8 -*-
'''
Python 3.10.6
4 September 2020
'''

import csv

class ThresholdCall:

    def __init__(self, serotype, prob):
        self.serotype = serotype
        self.prob = float(prob)
        self.threshold_file = 'ref/threshold/prob_thresholds.csv'
        self.thresholds = None
        self.valid = True
        self.run_check()

    def import_thresholds(self):
        threshold_tbl = {}
        with open(self.threshold_file, 'r') as inf:
            reader = csv.reader(inf)
            next(reader) #skip header
            for row in reader:
                threshold_tbl[row[0]] = float(row[1])
        self.thresholds = threshold_tbl

    def run_check(self):
        self.import_thresholds()
        # defaults to 0.5 prob threshold
        prob_min = self.thresholds.get(self.serotype, 0.5)
        if self.prob < prob_min:
            self.valid = False
