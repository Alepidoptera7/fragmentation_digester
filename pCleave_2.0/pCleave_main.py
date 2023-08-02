#!/usr/bin/env python3
# Name: Quin Lamothe (alamothe)
# Group Members: None

import sys
import re
import matplotlib.pyplot as pls
import numpy as np
import os
import textwrap

class cleaver:
    """This class contains all methods pertaining to developing a restriction digest. Methods include means to process
    brand name restriction enzymes, degenerate nucleotides and graphing/print capabilities for top and bottom DNA
    strands.

    Input: Sequences and headers from files.
    Output: Restriction digests as graphic and console display
    """
    def __init__(self, restriction, mark):
        """This member will hold the initial values of all objects contained.

        Input: Console commands
        Output: none
        """

        self.enzymes = {'EcoR1': 'GAATTC', 'BamHI': 'GGATTC', 'HIND': 'AAGCTT', 'Nspl': 'RCATG'}

        self.degenerate_nucleotides = {'R': ['G', 'A'], 'Y': ['T', 'C'], 'M': ['A', 'C'], 'K': ['G', 'T'],
                                       'S': ['G', 'C'], 'W': ['A', 'T'], 'H': ['A', 'C', 'T'], 'B': ['G', 'T', 'C'],
                                       'V': ['G', 'C', 'A'], 'D': ['G', 'A', 'T'], 'N': ['G', 'A', 'T', 'C']}

        self.cleaveList = []
        self.topStrand = []
        self.bottomStrand = []
        self.seq = ' '
        self.head = ' '
        self.topSeq = ' '
        self.bottomSeq = ' '
        self.degenerateList = []

        #command line parameters
        if restriction:
            self.restriction = restriction
        else:
            print('You must enter a restriction specificity with -r=, or a path with -p=')

        self.mark = mark

    def the_Cleave(self, seq):
        """This method will generate the locations of restriction sites for use in other functions. The restriction
        specificity provided by the user may include a number of degenerate nucleotides.

        Input: Headers and Sequences
        Output: Location of restriction sites along top and bottom strand
        """

        #searching for the user input in the dictionary of brand name restriction enzymes
        if self.restriction in self.enzymes.keys():
            self.restriction = self.enzymes[self.restriction]

        #searching the dictionary of degenerate nucleotides to discover if the user input contains a degenerate
        for degenerate in self.degenerate_nucleotides.keys():
            if degenerate in self.restriction:
                for item in self.degenerate_nucleotides[degenerate]:
                    for nuc in item:
                        self.degenerateList.append(self.restriction.replace(degenerate, nuc))

        #if there are infact degenerate nucleotides in the specificity, this will return true
        if self.degenerateList:
            for item in self.degenerateList:
                self.restriction = item
                #each degenerate outcome will be treated as specificity
                if self.restriction in seq:
                    self.cleaveList += [i.start() for i in re.finditer(self.restriction, seq)]

        #base case, searching the sequence for specificity
        if self.restriction in seq:
            self.cleaveList = [i.start() for i in re.finditer(self.restriction, seq)]
        else:
            self.cleaveList = []

        return self.cleaveList

    def gcContent(self):
        """This member is designed only to develop the GC content of a given sequence.

        Input: the current sequence
        Output: GC content rounded to two decimal places
        """
        g = 0
        c = 0
        totalNucs = 0
        #sequence is searched for G and C elements, with instances counted
        for nuc in self.seq:
            if nuc is 'G':
                g += 1
            if nuc is 'C':
                c += 1
            totalNucs += 1

        return round((g+c)/totalNucs * 100, 2)

    def coordinator(self, head, seq):
        """This objects is designed to moderate sequences delivered as parameters to the
        function which processes file data, in that it will pass the top strand through
        the function, and then the reverse complement which is manufactured within this object.

        Input: The file data.
        Output: the list of genes as obtained from processing both top and bottom strands.
        """

        self.head = head
        self.seq = seq.replace('-', '')

        self.topSeq = self.seq
        self.topStrand = self.the_Cleave(self.topSeq)

        #developing the bottom strand and passing it to the processing member
        self.bottomSeq = self.seq.lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G').replace(' ', '')
        self.bottomStrand = self.the_Cleave(self.bottomSeq)

    def print_and_pop(self):
        """This member is designed to display the restriction digest in graphic form as a line plot.
        Top and bottom strand will be displayed parallel. All informational output is displayed by this member.

        Input:Restriction sites on top and bottom strands.
        Output:Each sequence represented by its own line plot graphic.
        """

        topMark = ""
        bottomMark = ""

        #command line input can suppress console output
        if self.mark:
            for i in range(len(self.topSeq)):
                if i in self.topStrand:
                    topMark += '*'
                else:
                    topMark += '.'

            for i in range(len(self.bottomSeq)):
                if i in self.bottomStrand:
                    bottomMark += "*"
                else:
                    bottomMark += "."

        #lists are generated, each element containing 100 characters
        topMarkWrap = textwrap.wrap(topMark, 100)
        topWrap = textwrap.wrap(self.topSeq, 100)
        bottomWrap = textwrap.wrap(self.bottomSeq, 100)
        bottomWrapMark = textwrap.wrap(bottomMark, 100)

        #All items are printed together, such that sequences match corresponding restriciton site markers
        print(self.head, ' ', 'GC content: ', self.gcContent(), "%")
        for i, j, k, l in zip(topMarkWrap, topWrap, bottomWrap, bottomWrapMark):
            print(i)
            print(j)
            print(len(i)*'|')
            print(k)
            print(l)
            print(' ')

        #sequence header as window name
        pls.figure(num=self.head)

        #top strand graph
        pls.subplot(211)
        pls.vlines(self.topStrand, 0, len(self.seq), colors='r')
        pls.xlabel('\'5- Top Strand Restriction Sites -3\'')
        pls.subplots_adjust(hspace=1.0, wspace=1.0)

        if len(self.seq) < 100:
            pls.tick_params(which='minor', length=2.5)
            pls.xticks(np.arange(0, len(self.seq), step=10.0), rotation=-90)
        else:
            pass
        cur_axes = pls.gca()
        cur_axes.axes.get_xaxis().set_visible(True)
        cur_axes.axes.get_yaxis().set_visible(False)

        #bottom strand graph
        pls.subplot(212)
        pls.vlines(self.bottomStrand, 0, len(self.seq), colors='r')
        pls.xlabel('3\'- Bottom Strand Restriction Sites -5\'')
        pls.subplots_adjust(hspace=1.0, wspace=1.0)

        if len(self.seq) < 100:
            pls.tick_params(which='minor', length=2.5)
            pls.xticks(np.arange(0, len(self.seq), step=10.0), rotation=-90)
        else:
            pass
        cur_axes = pls.gca()
        cur_axes.axes.get_xaxis().set_visible(True)
        cur_axes.axes.get_yaxis().set_visible(False)

        pls.show()

class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output')

        self.parser.add_argument('-r', '--restriction', action='store', nargs='?', const=True, default=False,
                                 help='the specified restriction site')
        self.parser.add_argument('-m', '--mark', action='store', nargs='?', const=True, default=True,
                                 help='mark the restriction locations on print')
        self.parser.add_argument('-p', '--path', action='store', nargs='?', const=True, default=False,
                                 help='a path full of FASTA files')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
    print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line =\
                fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
                header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence

def main(inCL=None):
    """This main is designed to handle file and path reading, and command line inputs.

    Input: File and path names, command line options
    Output: File information, command line options
    """

    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    site = cleaver(myCommandLine.args.restriction, myCommandLine.args.mark)

    #this allows the use of a path, entered by command line
    if myCommandLine.args.path:
        path_get = myCommandLine.args.path
        path_x = path_get
        directory = os.fsencode(str(path_x))
        files = os.listdir(directory)
        os.chdir(directory)

        #each file in the path is placed into the reading class
        for file in files:
            myReader = FastAreader(file)
            for head, seq in myReader.readFasta():
                site.coordinator(head, seq)
                site.print_and_pop()
    else:
        #this is the base case of reading one file
        myReader = FastAreader()  # make sure to change this to use stdin
        for head, seq in myReader.readFasta():
            site.coordinator(head, seq)
            site.print_and_pop()

if __name__ == "__main__":
    main()  # delete the list when you want to run with STDIN
