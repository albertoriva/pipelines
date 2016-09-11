# An actor for ChIPseq experiments

import Logger
from MultiSampleActor import MultiSampleActor
from SampleCollection import SampleCollection

class ChipSeqActor(MultiSampleActor):
    
    def mergeBAMs(self, run=True, condition=True):
        """If a sample contains multiple replicates, merge its BAM files into a single one. Sets the 'bamfile' entry in the sample to the merged BAM file, or the existing one in case of no replicates. If `condition' is True, also merges the BAM files for all samples in each condition into a single BAM file."""
        nmerge = 0
        for smp in self.sc.samples:
            if len(smp['readsets']) == 1:
                rs = smp['readsets'][0]
                smp['bamfile'] = rs['bamfile']
            else:
                smp['bamfile'] = smp['name'] + ".bam"
                # Only do the merge if the merged file is missing or old
                if self.missingOrStale(smp['bamfile'], other=smp['readsets'][0]['bamfile']):
                    cmdline = "generic.qsub bamtools merge -out " + smp['bamfile']
                    for rs in smp['readsets']:
                        cmdline += " -in " + rs['bamfile']
                    cmdline += " module:bamtools"
                    if run:
                        self.submit(cmdline, done="merge.@.done")
                        nmerge += 1
        self.wait(("merge.@.done", nmerge))
        
        ## Merge BAMs for all samples in each condition
        nmerge = 0

        for cond in self.sc.conditions:

            ## First samples with default role
            condbams = self.sc.conditionBAMs(cond['name'])
            if len(condbams) == 1:
                cond['bamfile'] = condbams[0]
            else:
                cond['bamfile'] = cond['name'] + ".bam"
                if self.missingOrStale(cond['bamfile'], other=condbams[0]):
                    cmdline = "generic.qsub bamtools merge -out " + cond['bamfile']
                    for sb in condbams:
                        cmdline += " -in " + sb
                    cmdline += " module:bamtools"
                    if run:
                        self.submit(cmdline, done="merge.@.done")
                        nmerge += 1

            ## Then samples with 'input' role
            condibams = self.sc.conditionBAMs(cond['name'], role='input')
            if len(condbams) == 1:
                cond['ibamfile'] = condibams[0]
            else:
                cond['ibamfile'] = cond['name'] + "_input.bam"
                if self.missingOrStale(cond['ibamfile'], other=condbams[0]):
                    cmdline = "generic.qsub bamtools merge -out " + cond['ibamfile']
                    for sb in condibams:
                        cmdline += " -in " + sb
                    cmdline += " module:bamtools"
                    if run:
                        self.submit(cmdline, done="merge.@.done")
                        nmerge += 1

        self.wait(("merge.@.done", nmerge))

    def bamToWig(self):
        """Convert the BAM file(s) in each condition to WIG files."""
        nwig = 0
        for cond in self.sc.conditions:
            bams = self.sc.conditionBAMs(cond['name'], role=None)
            for bam in bams:
                wig = self.setFileExt(bam, ".wig", remove=[".bam"])
                self.submit("bamtowig.qsub {} {}".format(bam, wig), done="btw.@.done")
                nwig += 1
        return nwig
