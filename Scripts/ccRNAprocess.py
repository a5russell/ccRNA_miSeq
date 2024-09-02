
import pandas as pd


#hardcoded admittedly, start with samfiles, unzipped fastq don't worry about collapsing or calling distances, honestly want to peek at "raw-ish" data first
#make a list of pandas to concatenate at end. Improves runtime.
#requires all files are sorted
def assignReads(*, matchChart, TSOseq, R1fastq, R2fastq, R1samfile, R2samfile, outFile):
    seqs = pd.read_csv(matchChart)
    prevsam1name = ''
    prevsam2name = ''
    with open(R1samfile) as R1sf, open(R2samfile) as R2sf, open(R1fastq) as R1fq, open(R2fastq) as R2fq:
        #read in name, strip newline character
        R1name = R1fq.readline().split(' ')[0][1:]
        output = []
        readCount = 0
        prevsam1name = ''
        prevsam2name = ''
        while R1name != '':
            sam1 = R1sf.readline().split('\t')
            sam2 = R2sf.readline().split('\t')
            R2name = R2fq.readline().split(' ')[0][1:]
            sam1name = sam1[0]
            sam2name = sam2[0]
            #iterate through slight multimappers
            while(sam1name == prevsam1name):
                sam1 = R1sf.readline().split('\t')
                sam1name = sam1[0]
            while(sam2name == prevsam2name):
                sam2 = R2sf.readline().split('\t')
                sam2name = sam2[0]
            names = set([R1name, R2name, sam1name, sam2name])
            try:
                if R1name != R2name:
                    raise ValueError
            except:
                    print(' '.join(['names between files do not match or are out of in fastq, set of unique names are '] + [name for name in names]))
            #at times, reads may have been removed during trimming, so if all names not congruent, first check two fastq match as they should.
            #then, check both bamfiles if they match R1 or R2 names, and advance either as needed. Advance all until match
            if len(names) != 1:
                #advance until at least 3 things match
                while R1name != sam1name or R2name != sam2name:
                    R1name = R1fq.readline().split(' ')[0][1:]
                    R2name = R2fq.readline().split(' ')[0][1:]
                if sam1name != sam2name:
                    #see if sam2 is ahead
                    if R1name == sam1name:
                        while sam1name != sam2name and R1name != '':
                            #advance fastq files to next record, as well as sam1, until sam2 is matched
                            R2name = R2fq.readline()
                            R2name = R2fq.readline()
                            R2name = R2fq.readline().split(' ')[0][1:]
                            R1name = R1fq.readline()
                            R1name = R1fq.readline()
                            R1name = R1fq.readline().split(' ')[0][1:]
                            sam1 = R1sf.readline().split('\t')
                            sam1name = sam1[0]
                    else:
                        while sam1name != sam2name and R1name != '':
                            #advance fastq files to next record, as well as sam1, until sam1 is matched
                            R2name = R2fq.readline()
                            R2name = R2fq.readline()
                            R2name = R2fq.readline().split(' ')[0][1:]
                            R1name = R1fq.readline()
                            R1name = R1fq.readline()
                            R1name = R1fq.readline().split(' ')[0][1:]
                            sam2 = R2sf.readline().split('\t')
                            sam2name = sam2[0]
            if R1name != '':
                seg1 = sam1[2]
                seg2 = sam2[2]
                r1 = R1fq.readline()
                r2 = R2fq.readline()
                #exclude unmapped

                if seg1 != '*' and seg2 != '*':
                    r1end = seqs[seqs.Segment == seg1].R1_seq.iloc[0]
                    r2end = seqs[seqs.Segment == seg2].R2_seq.iloc[0]
                    position1 = r1.find(r1end)
                    position2 = r2.find(r2end)
                    TSO = r1.find(TSOseq)
                    #require that there are ends and this wasnt an internal thing
                    if position1 != -1 and position2 != -1:
                        position1 = position1 + len(r1end)
                        position2 = position2 + len(r2end)
                        if TSO != -1:
                            TSO = TSO - position1
                        output += [pd.DataFrame({'Segment_1':[seg1], 'Segment_2':[seg2], 'R1_dist':[position1], 'R2_dist':[position2],
                                    'TSO_space':[TSO]})]
                prevsam1name = sam1name
                prevsam2name = sam2name
                R2name = R2fq.readline()
                R2name = R2fq.readline()
                R1name = R1fq.readline()
                R1name = R1fq.readline()
                R1name = R1fq.readline().split(' ')[0][1:]

                readCount += 1
        output = pd.concat(output)
        output.to_csv(outFile)
