#!/usr/bin/env python

# (c) 2015, A. Riva, DiBiG, ICBR Bioinformatics
# University of Florida

import sys
import os.path
import ConfigParser

# Utility

def splitCommas(s):
    """Split a string s containing items separated by commas into a list
of strings. Spaces before or after each string are removed."""
    items = s.split(",")
    return [ i.strip(" ") for i in items ]

def cleanEOL(s):
    """Remove \r from the end of string `s'."""
    return s.rstrip("\r")

# This class implements a collection of fastq files from multiple samples,
# each one possibly having replicates. Samples can also be grouped into
# conditions, and conditions can be compared to each other in contrasts.
# The contents of a SampleCollection are initialized using a configuration file.

# The basic unit is a readset. Each readset contains a single fastq file, or
# two fastq files that should be handled together (in paired-end mode).
# A sample contains one or more readsets (interpreted as replicates).

class SampleCollection():
    conf = None                 # Configuration object
    readsets = []               # List of all readsets
    nreadsets = 0               # Number of readsets
    samples = []                # List of all samples
    nsamples = 0                # Number of samples
    conditions = []             # List of all conditions
    nconditions = 0             # Number of conditions
    contrasts = []              # List of all contrasts
    ncontrasts = 0              # Number of contrasts
    # For use as an iterator
    _current = 0

    def __init__(self, conf):
        """Initialize this SampleCollection from the configuration object `conf'."""
        self.conf = conf
        self.readsets = []
        self.samples = []
        self.conditions = []
        self.contrasts = []
        self.initializeConditions()
        if self.nconditions == 0:
            self.initializeSamples()
        self.parseContrasts()

    def describe(self):
        print "{} conditions: {}".format(self.nconditions, self.conditions)
        print "{} samples:    {}".format(self.nsamples, self.samples)
        print "{} readsets:   {}".format(self.nreadsets, self.readsets)
        print "{} contrasts:  {}".format(self.ncontrasts, self.contrasts)

    def getConf(self, entry, section="General", required=False):
        try:
            return self.conf.get(section, entry)
        except ConfigParser.NoOptionError:
            if required:
                print "Configuration error: entry {} is required in section [{}]".format(entry, section)
                sys.exit()
            else:
                return None

    def fixPath(self, path):
        """Add ../ in front of `path' unless it is absolute."""
        if path[0] == "/":
            return path
        else:
            return "../" + path

    def checkPath(self, p):
        """If path `p' exists, return the result of fixPath() on it. Otherwise, signal an error and exit."""
        if p == None:
            return p
        elif not os.path.isfile(p):
            print "Error: file {} does not exist or is not readable.".format(p)
            sys.exit()
        else:
            return self.fixPath(p)

### Conditions

    def initializeConditions(self):
        """Initialize conditions from the configuration file."""
        condnames = self.getConf("conditions")
        if condnames != None:
            condnames = splitCommas(condnames)
            for c in condnames:
                # print "condition: " + c
                condsamples = self.getConf("samples", section=c)
                if condsamples == None:
                    break
                condsamples = splitCommas(condsamples)
                condition = {'name': c,
                             'samples': condsamples}
                for s in condsamples:
                    # print "  sample: " + s
                    if not self.findSample(s):
                        self.addSample(s, role='default')
                condinputs = self.getConf("inputs", section=c)
                if condinputs != None:
                    condinputs = splitCommas(condinputs)
                    condition['inputs'] = condinputs
                    for s in condinputs:
                        # print "  input: " + s
                        if not self.findSample(s):
                            self.addSample(s, role='input')
                    condition['samples'] = condsamples + condinputs

                self.conditions.append(condition)
                self.nconditions += 1

    def findCondition(self, name):
        """Return the condition called `name'."""
        for c in self.conditions:
            if c['name'] == name:
                return c
        return None

    def conditionSamples(self, name, role='default'):
        """Return all samples for condition `name'. By default, only samples
with 'default' role are returned. The `role' argument can be used to return
samples with a different role, or all if `role' is None."""
        if type(name).__name__ == 'str':
            c = self.findCondition(name)
        else:
            c = name
        if c == None:
            return None
        else:
            samples = []
            for s in c['samples']:
                smp = self.findSample(s)
                if role == None or role == smp['role']:
                    samples.append(smp)
            return samples

    def conditionBAMs(self, name, role='default', key='bam'):
        """Returns the list of BAM files for all the samples in this condition. By default, 
only BAM files for samples with 'default' role are returned. The `role' argument can be used 
to return samples with a different role, or all if `role' is None."""
        return [ s[key] for s in self.conditionSamples(name, role=role) ]

### Samples

    def initializeSamples(self):
        """If conditions are not being used, look for samples in the General section."""
        samplenames = self.getConf("samples")
        if samplenames != None:
            samplenames = splitCommas(samplenames)
            for s in samplenames:
                self.addSample(s)

    def addSample(self, name, role='default'):
        # print "adding sample {} with role {}".format(name, role)
        sample = {'name': name, 'role': role, 'readsets': self.parseReadsets(name)}
        self.samples.append(sample)
        self.nsamples += 1
        return sample

    def parseReadsets(self, samplename):
        rs = []

        # First try single-end
        f1 = self.getConf("fastq", samplename)
        if f1:
            r = {'name': "{}_r1".format(samplename),
                 'left': cleanEOL(f1),
                 'paired': False,
                 'bad': False}
            self.readsets.append(r)
            self.nreadsets += 1
            rs.append(r)
            return rs

        # Single-end with replicates
        i = 1
        while True:
            e1 = "r{}_fastq".format(i)
            f1 = self.getConf(e1, samplename)
            if f1 == None:
                break
            r = {'name': "{}_r{}".format(samplename, i),
                 'left': cleanEOL(f1),
                 'paired': False,
                 'bad': False}
            self.readsets.append(r)
            self.nreadsets += 1
            rs.append(r)
            i += 1
        if len(rs) > 0:
            return rs

        # Paired-end without replicates
        f1 = self.getConf("left", samplename)
        f2 = self.getConf("right", samplename)
        if f1 and f2:
            r = {'name': "{}_r1".format(samplename),
                 'left': cleanEOL(f1),
                 'right': cleanEOL(f2),
                 'paired': True,
                 'bad': False}
            self.readsets.append(r)
            self.nreadsets += 1
            rs.append(r)
            return rs

        # Paired-end with replicates
        i = 1
        while True:
            e1 = "r{}_left".format(i)
            e2 = "r{}_right".format(i)
            f1 = self.getConf(e1, samplename)
            f2 = self.getConf(e2, samplename)
            if f1 == None or f2 == None:
                break
            r = {'name': "{}_r{}".format(samplename, i),
                 'left': cleanEOL(f1),
                 'right': cleanEOL(f2),
                 'paired': True,
                 'bad': False}
            self.readsets.append(r)
            self.nreadsets += 1
            rs.append(r)
            i +=1
        return rs

    def findSample(self, name):
        """Return the sample called `name'."""
        for smp in self.samples:
            if smp['name'] == name:
                return smp
        return None

    def sampleReadsets(self, name):
        """Return the readsets of sample `name'."""
        smp = self.findSample(name)
        if smp == None:
            return None
        else:
            return smp['readsets']

### Readsets

    def findReadset(self, name):
        for rs in self.readsets:
            if rs['name'] == name:
                return rs
        return None

### Contrasts

    def splitContrast(self, c):
        p = c.find("^")
        if p == -1:
            print "Configuration error: contrast `{}' should have the form testsample^controlsample.".format(c)
            sys.exit()
        test = c[:p]
        ctrl = c[p+1:]
        name = "{}.vs.{}".format(test, ctrl)
        return {'test': test, 'control': ctrl, 'name': name, 'label': c}

    def parseContrasts(self):
        base = self.getConf("contrasts")
        if base != None:
            base = splitCommas(base)
            self.contrasts = [ self.splitContrast(c) for c in base ]
            self.ncontrasts = len(self.contrasts)
            return self.contrasts
        else:
            return None

    def showSamples(self):
        print "{} total samples".format(self.nsamples)
        for s in self.samples:
            print "\nSample '{}' ({} readsets)".format(s['name'], len(s['readsets']))
            for rs in s['readsets']:
                if rs['paired']:
                    print "  {}".format(rs['name'])
                    print "    left = {}\n    right = {}".format(rs['left'], rs['right'])
                else:
                    print "  {}".format(rs['name'])
                    print "    fastq = {}\n".format(rs['left'])

# Verify that the SampleCollection is consistent

    def verify(self, verbose=True):
        """Verify that this SampleCollection is consistent: all fastq files should exist, all samples
referenced in conditions should exist, all samples or conditions referenced in contrasts should exist."""
        good = True
        
        ## Check samples in conditions
        for c in self.conditions:
            for s in c['samples']:
                if not self.findSample(s):
                    good = False
                    if verbose:
                        print "Condition {} contains non-existent sample {}.".format(c['name'], s)

        ## Check samples or conditions in contrasts
        for c in self.contrasts:
            ctrl = c['control']
            test = c['test']
            if not (self.findSample(ctrl) or self.findCondition(ctrl)):
                good = False
                if verbose:
                    print "Contrast {}^{} contains non-existent sample or condition {}.".format(ctrl, test, ctrl)
            if not (self.findSample(test) or self.findCondition(test)):
                good = False
                if verbose:
                    print "Contrast {}^{} contains non-existent sample or condition {}.".format(ctrl, test, test)
        
        ## Check fastq files in readsets
        for rs in self.readsets:
            if not os.path.isfile(rs['left']):
                good = False
                if verbose:
                    print "Readset {} references missing file {}".format(rs['name'], rs['left'])
            if rs['paired']:
                if not os.path.isfile(rs['right']):
                    good = False
                    if verbose:
                        print "Readset {} references missing file {}".format(rs['name'], rs['right'])
        return good

# Using the collection as an iterator over (good) readsets

    def __iter__(self):
        self.__idx = 0
        return self

    def next(self):
        while True:
            if self.__idx >= self.nreadsets:
                raise StopIteration
            x = self.readsets[self.__idx]
            self.__idx += 1
            if not ('bad' in x and x['bad']):
                return x

if __name__ == "__main__":
    conf = ConfigParser.ConfigParser()
    conf.read(sys.argv[1])
    sc = SampleCollection(conf)
    sc.describe()
    print sc.verify()
    sc.showSamples()
