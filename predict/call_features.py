# -*- coding: utf-8 -*-
'''
Python 3.10.6
24 August 2020
'''

import os
import re
import xml.etree.ElementTree as ET

import screed

class ORFchecker:

    def __init__(self, fasta, serogroup):
        self.fasta = fasta
        self.serogroup = serogroup #15,18,35
        self.blst = 'dummy.xml'
        self.serotype = serogroup
        self.flag = ''
        self.run()

    # run blast search for causal gene
    def blast(self):
        ref = 'ref/blast/causal_genes.fasta'
        self.blst = self.fasta.split('.f')[0] + '_blast.xml'
        cmd = 'blastn -db {0} -query {1} -out {2} -outfmt 5'.format(ref, self.fasta, self.blst)
        os.system(cmd)

    # extract aligned sequence from blast output
    def pull_results(self):
        tree = ET.parse(self.blst)
        root = tree.getroot()
        # path to blast results
        hp = './BlastOutput_iterations/Iteration/Iteration_hits/Hit/Hit_hsps/Hsp'
        attributes = ('Hsp_hit-from', 'Hsp_hit-to', 'Hsp_qseq')
        outputs = []
        for attribute in attributes:
            path = '{0}/{1}'.format(hp, attribute)
            element = ET.tostring(root.findall(path)[0])
            output = re.search('>(.+?)<', str(element)).group(1)
            outputs.append(output)
        blast_result = {'start':int(outputs[0]),
                        'end':int(outputs[1]),
                        'seq':outputs[2]}
        return blast_result

    # reverse complement sequence if needed and remove dashes
    def curate_sequence(self, blast_result):
        seq = blast_result['seq'].replace('-', '')
        if blast_result['start'] > blast_result['end']:
            seq = screed.rc(seq)
        return seq

    # checks for presence of a stop codon
    def integrity_check(self, sequence):
        stops = {'TAA', 'TGA', 'TAG'}
        complete = 1
        if len(sequence)%3 != 0:
            complete = 0
        i = 0
        while i <= len(sequence)-3 and complete:
            codon = sequence[i:i+3]
            if codon in stops:
                complete = 0
            i += 3
        return complete

    # remove blast output files
    def cleanup(self):
        if os.path.exists(self.blst):
            os.system('rm ' + self.blst)

    def run(self):
        try:
            self.blast()
            blast_result = self.pull_results()
            corrected = self.curate_sequence(blast_result)
            prot_functional = self.integrity_check(corrected)
            st_tbl = {'15':{1:'15B', 0:'15C'}, '18':{1:'18C', 0:'18B'}, '35':{1:'35B', 0:'35D'}}
            soi = st_tbl[self.serogroup] #serotypes of interest
            self.serotype = soi[prot_functional]
        except:
            self.flag = 'ORF PREDICTION FAILED: SEROTYPE AMBIGUOUS;'
        finally:
            self.cleanup()
