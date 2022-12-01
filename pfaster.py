# -*- coding: utf-8 -*-
'''
Python 3.10.6
12 September 2022
'''

from predict import Mash, call_features, threshold, model_predict
from tools import cmd_parser
from tools import exporters as exp


def screen_fasta(fasta, db = 'sketch_k70.pkl'):
    screener = Mash.MashScreen('ref/sketch/{}'.format(db))
    screener.screen(fasta)
    matches = screener.ref_counts
    return list(matches.values())

def run_model(mash_results, model = 'model.rfm'):
    model = 'model/{}'.format(model)
    st = model_predict.RFClassify(mash_results, model).serotype
    return st # tuple - serotype, probability

def call_serotype(fasta, outdir):
    exp.update_log('screening ' + fasta)
    #format: [isolate, serotype, probability, flags]
    prediction = [fasta, 'NT', 0, ''] #initialize as NT
    try:
        hash_matches = screen_fasta(fasta) # run Mash screen
        exp.update_log('MashScreen completed')
        pred = run_model(hash_matches) # pass results to RF model
        exp.update_log('Model probabilities computed')
        pred_sero = pred[0][2:] # removes 'Pn' prefix
        prob = round(pred[1], 2)
        if pred_sero in {'24B', '24F'}:
            prediction[3] = 'AMBIGUOUS - 24B OR 24F'
        # run ORF check step
        ambig_types = ('35B', '35D', '18B', '18C', '15B', '15C')
        if pred_sero in ambig_types:
            msg = 'ambiguous serotype for {} - checking ORF'.format(fasta)
            exp.update_log(msg)
            feature_check = call_features.ORFchecker(fasta, pred_sero[:2])
            pred_sero = feature_check.serotype
            prediction[3] = prediction[3] + feature_check.flag
        prediction[1] = pred_sero
        prediction[2] = prob
    except:
        prediction = [fasta, 'not typed', 'N/A', 'failed to predict serotype']
        return prediction
    #probability thresholding
    high_confidence = threshold.ThresholdCall(prediction[1], prediction[2]).valid
    if not high_confidence:
        low_conf_flag = '{} called at low confidence;'.format(prediction[1])
        prediction[3] = prediction[3] + low_conf_flag
        prediction[1] = 'NT'
        update_log('confidence threshold not met for {}'.format(prediction[0]))
    exp.update_log('{0} called as {1}'.format(prediction[0], prediction[1]))
    if outdir: exp.write_results(outdir, prediction)
    return prediction

def main():
    args = cmd_parser.CommandLine().args
    st_prediction = call_serotype(**args)
    exp.print_prediction(st_prediction)

if __name__ == '__main__':
    main()
