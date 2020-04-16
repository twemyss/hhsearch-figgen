#!/usr/bin/env python

import math
import requests
import numpy as np

class ParsedHMM:
    """ Represents a parsed HMM file

    This stores representations of HHSuite HMMs.
    """

    def __init__(self, config):
        """ Parse HMM file given in configuration, and load into memory

        :param config: parsed JSON configuration object
        """
        # Define the important parameters as state variable
        self.hmm_string = ''
        self.alphabet = []
        self.states = []
        self.probs = []
        self.state_probs = []
        self.clustal_colours = []
        self.height_array = []
        self.ss = ''
        self.ss_probs = ''
        self.nulls = []

        # Process HMM
        print(config['master']['hmm_file'])
        with open(config['master']['hmm_file']) as hmm:
            in_hmm_section = False
            in_psipred_ss = False
            in_psipred_ss_probs = False
            for line in hmm:
                if line.startswith("LENG"):
                    self.length = int(line.split()[1])
                if line.startswith("NAME"):
                    if ".fa" in line or "_" in line:
                        self.name = config['master']['name']
                    else:
                        self.name = line.split()[1]
                if line.startswith("FILT"):
                    self.num_seqs = float(line.split()[1])
                if line.startswith("HMM"):
                    in_hmm_section = True
                if in_hmm_section:
                    self.hmm_string += line

                # Process null
                if line.startswith("NULL"):
                    splitLine = line.split()
                    for i in range(1, len(splitLine)):
                        self.nulls.append(2 ** (int(splitLine[i]) / -1000))

                # Clear current section if we're at a new FASTA header
                if line.startswith(">"):
                    in_psipred_ss = False
                    in_psipred_ss_probs = False

                # NB, unlike for HMM, we don't want the ss to include the FASTA header, so we check section before
                if in_psipred_ss:
                    self.ss += line
                if "ss_pred PSIPRED predicted secondary structure" in line:
                    in_psipred_ss = True
                if in_psipred_ss_probs:
                    self.ss_probs += line
                if "ss_conf PSIPRED confidence values" in line:
                    in_psipred_ss_probs = True

        # Remove linebreaks/whitespace from ss and ss_probs
        self.ss = ''.join(self.ss.split())
        self.ss_probs = ''.join(self.ss_probs.split())

        # Now process the HMM string we have
        hmm_lines = self.hmm_string.split("\n")
        for idx, line in enumerate(hmm_lines):
            # Extract the fields within the HMM file
            if line.startswith("HMM"):
                for hmm_line in range(len(line.split())):
                    if hmm_line != 0:
                        self.alphabet.append(line.split()[hmm_line])
                for hmm_line in range(len(hmm_lines[idx + 1].split())):
                    self.states.append(hmm_lines[idx + 1].split()[hmm_line])

        # Now move on to the match states themselves
        for l in range(self.length):
            alphabet_line = hmm_lines[l * 3 + 3].split()
            stats_line = hmm_lines[l * 3 + 4].split()
            current_alphabet_probs = []
            for a in range(len(alphabet_line)):
                if a > 1 and a < len(alphabet_line) - 1:
                    if alphabet_line[a] is '*':
                        f = 0
                    else:
                        f = 2 ** (int(alphabet_line[a]) / -1000)
                    current_alphabet_probs.append(f)
            self.probs.append(current_alphabet_probs)
            current_states = []
            for s in range(len(stats_line)):
                if stats_line[s] is '*':
                    f = 0
                else:
                    f = 2 ** (int(stats_line[s]) / -1000)
                current_states.append(f)
            self.state_probs.append(current_states)


        if 'output' not in config or config['output']['conservation_plot']['type'] == 'traditional':
            # Go through each position and calculate the height
            for height_idx in range(len(self.probs)):
                entropy = self.get_shannon_entropy(height_idx)

                height = entropy
                self.height_array.append(height)
            for aa_idx in range(len(self.probs)):
                height = self.height_array[aa_idx]
                heights = [0, 0, 0, 0, 0, 0, 0, 0]
                for j in range(len(self.probs[aa_idx])):
                    aa_prob = self.probs[aa_idx][j]
                    aa_name = self.alphabet[j]
                    for k in range(len(config['colours'])):
                        if aa_name in config['colours'][k]['aa']:
                            if aa_prob > self.nulls[j]:
                                heights[k] += aa_prob * height
                self.clustal_colours.append(heights)
            print("Heights:")
            print(self.height_array)
            print("Max height:")
            print(np.amax(self.height_array))
            print("Background:")
            print(self.nulls)
            print(np.sum(self.nulls))
        elif config['output']['conservation_plot']['type'] == 'skylign':
            # Get the HMM file and send it to skylign via post, using arguments
            print("Making request")
            r = requests.post('http://skylign.org/',
                              params={'file': open(config['master']['alignment_a3m'], 'r')},
                              headers={"Accept": "Application/json"}, timeout=5)
            print(r.text)

        return

    def getKullbackLeiblerDistance(self, position):
        return

    def get_shannon_entropy(self, height_idx):
        position = self.probs[height_idx]
        shannon_entropy = 0
        for i in range(len(position)):
            if position[i] > 0:
                if position[i] > self.nulls[i]:
                    shannon_entropy += position[i] * math.log2(position[i] / self.nulls[i])
        return shannon_entropy

    def get_small_sample_correction(self):
        return 1 / math.log(2) * (20 - 1) / (2 * self.num_seqs)
