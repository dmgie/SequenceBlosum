# Given a file with multiple FASTA sequences, calculate the log odds ratio for each pair of amno acids

import argparse
import math
from itertools import combinations

parser = argparse.ArgumentParser(description = 'Calculate the log odds ratio for each pair of amino acids')
parser.add_argument('filename', help = 'FASTA file')
parser.add_argument('-t','--temp-scale', help = 'Whether to enable the temperature scale', action = 'store_true')
args = parser.parse_args()

# Simple sequence class
class Sequence:
    """
    A simple sequence class that can be used to store a sequence and its header, as well as the optimum temperature for that organism
    """
    def __init__(self, header, sequence, temp_scale):
        self.header = header
        self.sequence = sequence
        self.temp = temp_scale

    def getHeader(self):
        return self.header

    def getSequence(self):
        return self.sequence


# FASTA reader that tracks sequence as Sequence objects
def readFASTA(filename):
    # Create a dictionary to store the sequences
    sequences = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Loop through the lines
        for line in lines:
            # If the line is a header, create a new sequence object
            if line[0] == '>':
                header = line[1:].strip()
                # number at the end of the header is the temperature
                if args.temp_scale:
                    temp = int(header.split()[-1])
                    header = header.split()[0]
                else:
                    temp = 0
                sequences[header] = Sequence(header, '', temp)
            # If the line is a sequence, add it to the current sequence object
            else:
                sequences[header].sequence += line.strip()
    # print('Read', len(sequences), 'sequences')
    return sequences


def calculateTempScaling(sequences):
    """
    Given a dictionary of sequences, calculate the temperature scaling factor for the set of sequences 
    Input: dictionary of Sequence objects
    Output: temperature scaling factor
    """
    # Get all the temperatures from the sequences 
    seqTemps = []
    for seq in sequences.values():
        seqTemps.append(seq.temp)

    # Calculate the variance of the temperatures
    meanTemp = sum(seqTemps) / len(seqTemps)
    variance = 0
    for temp in seqTemps:
        variance += (temp - meanTemp) ** 2 
    variance /= len(seqTemps)

    return variance

# Calculate log odds ratio
def calcLogOdds(sequences):
    """
    Calculate the log odds ratio according to the BLOSUM matrix formula from (Heinkoff and Henikoff (1992))
    Input: dictionary of Sequence objects
    Output: dictionary of log odds ratios
    """
    # Get the length of the aligned sequences by looking at the first sequence 
    seqLength = len(list(sequences.values())[0].sequence)


    # Sequences as list
    sequenceList = sequences.values()

    ### Joint probability of amino acids pairs and counts of each amino acids
    aaCount = {}
    pairCount = {}
    for i in range(seqLength):
        # Count number of occurences of each amino acids
        for seq in sequenceList:
            aa = seq.sequence[i]
            if aa not in aaCount:
                aaCount[aa] = 0
            aaCount[aa] += 1
        # Go through each pair of sequences and count pairs occurence
        for seq1, seq2 in combinations(sequenceList, 2):
            # Count the pairs
            aa1 = seq1.sequence[i]
            aa2 = seq2.sequence[i]
            pair = (aa1, aa2)
            pair = tuple(sorted(pair))
            if pair not in pairCount:
                pairCount[pair] = 0
            pairCount[pair] += 1

    
    print(aaCount)
    totalAAs = sum(aaCount.values())
    aaFreq = {}
    for aa in aaCount:
        aaFreq[aa] = aaCount[aa] / totalAAs
    
    # Calculate total number of posible pairs => len * n(n-1)/2 (or rather n choose 2)
    # Normalise frequency of each pair of amino acids i.e (num. seq choose 2) = observed frequency
    totalPairs = seqLength * len(sequenceList) * (len(sequenceList) - 1) / 2
    pairFreq = {}
    for pair in pairCount:
        pairFreq[pair] = pairCount[pair] / totalPairs

    logOddsRatios = {}
    # Calculate log odds ratio for each pair of amino acids
    for pair in pairFreq:
        aa1 = pair[0]
        aa2 = pair[1]
        expectedFreq = 0
        if aa1 == aa2:
            expectedFreq = aaFreq[aa1] * aaFreq[aa2]
        else:
            expectedFreq = 2 * aaFreq[aa1] * aaFreq[aa2]
        # logOddsRatio = 2 * math.log(pairFreq[pair] / (expectedFreq),2) 
        # TODO: Replace 2 with a scaling factor
        logOddsRatio = 2 * math.log(pairFreq[pair] / (expectedFreq)) 
        logOddsRatios[pair] = logOddsRatio

    # print(logOddsRatio)
    return logOddsRatios


# Print the log odds ratio between every amino acid in a nice matrix
def printLogOddsRatio(logOddsRatio):
    """
    Pretty print the substitution matrix in a nice format to the terminal, for easy piping to a file 
    Input: dictionary of log odds logOddsRatio 
    Output: None (prints to terminal)
    """
    # Get the list of amino acids contained from our fasta files
    aaList = []
    for pair in logOddsRatio:
        if pair[0] not in aaList:
            aaList.append(pair[0])
        if pair[1] not in aaList:
            aaList.append(pair[1])
    aaList.sort()
    # Print the header
    print(' ', end = '')
    for aa in aaList:
        print('\t', aa, end = '')
    print()
    # Print the matrix in a tabular format
    for aa1 in aaList:
        for aa2 in aaList:
            # If the pair is in the dictionary, print the log odds ratio, rounded to an int
            if (aa1, aa2) in logOddsRatio:
                print('\t', int(round(logOddsRatio[(aa1, aa2)], 0)), end = '')
            else:
                print('\t', 0, end = '')
        print()



# Main function
def main():
    sequences = readFASTA(args.filename)
    logOddsRatio = calcLogOdds(sequences)
    # print(calculateTempScaling(sequences))
    printLogOddsRatio(logOddsRatio)

if __name__ == '__main__':
    main()
