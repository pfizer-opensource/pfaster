# -*- coding: utf-8 -*-
"""
Python 3.10.6
30 April 2020
"""

import pickle

import numpy as np
from sklearn.ensemble import RandomForestClassifier


class RFClassify:

    def __init__(self, mash_counts, model):
        self.mash_counts = mash_counts
        self.model = model
        self.serotype = None
        self.predict()

    def load_model(self):
        with open(self.model, 'rb') as rf_model:
            classifier = pickle.load(rf_model)
        return classifier

    def predict(self):
        classifier = self.load_model()
        features = np.array( list(self.mash_counts) ).reshape(1, -1)
        st = classifier.predict(features)[0]
        prob = np.max( classifier.predict_proba(features) )
        self.serotype = (st, prob)
