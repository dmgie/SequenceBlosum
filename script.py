# Given a file with multiple FASTA sequences, calculate the log odds ratio for each pair of amno acids

import argparse
import math

parser = argparse.ArgumentParser(description = 'Calculate the log odds ratio for each pair of amino acids')
parser.add_argument('filename', help = 'FASTA file with multiple sequences')
args = parser.parse_args()

# Simple sequence class
class Sequence:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence

    def getHeader(self):
        return self.header

    def getSequence(self):
        return self.sequence


# FASTA reader that tracks sequence as Sequence objects
def readFASTA(filename):
    # Read the file 
    sequences = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        # Create a dictionary to store the sequences
        # Loop through the lines
        for line in lines:
            # If the line is a header, create a new sequence object
            if line[0] == '>':
                header = line[1:].strip()
                sequences[header] = Sequence(header, '')
            # If the line is a sequence, add it to the current sequence object
            else:
                sequences[header].sequence += line.strip()
    print('Read', len(sequences), 'sequences')
    return sequences


# Create a function that calculates the log odds ratio for each pair of amino acids
# It is the join probability divided by the probability of each letter
def calculateLogOddsRatio(sequences):
    # Create a dictionary to store the log odds ratio
    logOddsRatio = {}
    # Loop through the sequences
    for sequence in sequences.values():
        # Loop through the sequence
        for i in range(len(sequence.sequence) - 1):
            # Get the current and next amino acids
            currentAA = sequence.sequence[i]
            nextAA = sequence.sequence[i + 1]
            # If the current and next amino acids are not in the dictionary, add them
            if (currentAA, nextAA) not in logOddsRatio:
                logOddsRatio[(currentAA, nextAA)] = 0
            # Increment the log odds ratio for the current and next amino acids
            logOddsRatio[(currentAA, nextAA)] += 1
    # Calculate the log odds ratio for each pair of amino acids
    for pair in logOddsRatio:
        # Get the number of times the current and next amino acids appear together
        currentNextAA = logOddsRatio[pair]
        # Get the number of times the current amino acid appears
        currentAA = 0
        for pair2 in logOddsRatio:
            if pair2[0] == pair[0]:
                currentAA += logOddsRatio[pair2]
        # Get the number of times the next amino acid appears
        nextAA = 0
        for pair2 in logOddsRatio:
            if pair2[1] == pair[1]:
                nextAA += logOddsRatio[pair2]
        # Calculate log odds ratio
        logOddsRatio[pair] = math.log((currentNextAA / (currentAA * nextAA)) + 1)
    return logOddsRatio


# Another way to do it
def logOdds(sequences):
    # Count occurence of each amino acid 
    aaCount = {}
    for sequence in sequences.values():
        for aa in sequence.sequence:
            if aa not in aaCount:
                aaCount[aa] = 0
            aaCount[aa] += 1 
    # Count occurence of each pair of amino acids occuring together in different sequences
    aaPairCount = {}
    for sequence in sequences.values():
        for i in range(len(sequence.sequence) - 1):
            aaPair = sequence.sequence[i:i + 2]
            if aaPair not in aaPairCount:
                aaPairCount[aaPair] = 0
            aaPairCount[aaPair] += 1
    # Calculate log odds ratio 
    logOddsRatio = {}
    for aaPair in aaPairCount:
        logOddsRatio[aaPair] = math.log((aaPairCount[aaPair] / (aaCount[aaPair[0]] * aaCount[aaPair[1]])))


    print(logOddsRatio)
    return logOddsRatio
    # For each pair of possible amino acids
    # - Calculate the joint probability of these occuring between different sequences



# Print log odds ratio between each pair of amino acids as a matrix
def printLogOddsRatio(logOddsRatio):
    # Get the list of amino acids
    aminoAcids = []
    for pair in logOddsRatio:
        if pair[0] not in aminoAcids:
            aminoAcids.append(pair[0])
        if pair[1] not in aminoAcids:
            aminoAcids.append(pair[1])
    # Sort the list of amino amino 
    aminoAcids.sort()

    # Print the header 
    print(' ', end = '')
    for aa in aminoAcids:
        print(aa, end = ' ')
    print()


    # Print the log odds ratio for each pair of amino amino acids 
    for aa1 in aminoAcids:
        for aa2 in aminoAcids:
            if aa1+aa2 in logOddsRatio:
                print(logOddsRatio[aa1+aa2], end = ' ')
            else:
                print(0, end = ' ')



    # Old function
    # for aa1 in 'ACDEFGHIKLMNPQRSTVWY':
    #     for aa2 in 'ACDEFGHIKLMNPQRSTVWY':
    #         print(logOddsRatio[(aa1, aa2)], end = '\t')



# Main function
def main():
    # Read the FASTA file 
    sequences = readFASTA(args.filename)
    # Calculate the log odds ratio for each pair of amino acids
    logOddsRatio = logOdds(sequences)
    # Print log odds ratio between each pair of amino acids as a matrix
    printLogOddsRatio(logOddsRatio)

if __name__ == '__main__':
    main()
