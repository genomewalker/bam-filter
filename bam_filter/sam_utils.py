import pysam
import numpy as np
import re
import pandas as pd
from multiprocessing import Pool
from scipy import stats
import pysam
import tqdm


# Function to calculate evenness of coverage
def coverage_evenness(coverage):
    """
    Calculate the evenness of coverage
    """
    # get coverage evenness
    # covEvenness = (
    #     (breadth * 100) / (100.0 * expBreadth * expBreadth) if expBreadth > 0 else 0
    # )
    # C = mean(X)
    # D2 = X[X<=C]
    # N = len(X)
    # n = len(D2)
    # E = 1 - (n - sum(D2) / C) / N
    C = float(round(np.mean(coverage)))
    D2 = [x for x in coverage if x <= C]
    if len(D2) == 0:  # pragma: no cover
        covEvenness = 1.0
    else:
        covEvenness = 1.0 - (len(D2) - sum(D2) / C) / len(coverage)

    return covEvenness


def extractFromBam(params):
    """
    Worker function per chromosome
    loop over a bam file and create tuple with lists containing metrics:
    two definitions of the edit distances to the reference genome scaled by aligned read length
    """
    bam, chromosome = params
    samfile = pysam.AlignmentFile(bam, "rb")
    editDistancesNM = []
    editDistancesMD = []
    aniNM = []
    aniMD = []
    readLength = []
    readAlignedLength = []
    nReads = 0
    refLen = int(samfile.get_reference_length(chromosome))

    for read in samfile.fetch(reference=chromosome, multiple_iterators=False):
        nReads += 1
        aniNM.append((1 - ((read.get_tag("NM") / read.query_alignment_length))) * 100)
        aniMD.append(
            (
                (
                    1
                    - (
                        sum(
                            [
                                len(item)
                                for item in re.split("[0-9^]", read.get_tag("MD"))
                            ]
                        )
                        + sum(  # Parse MD string to get mismatches/deletions
                            [item[1] for item in read.cigartuples if item[0] == 1]
                        )
                    )  # Parse cigar to get insertions
                    / read.query_alignment_length
                )
            )
            * 100
        )

        editDistancesNM.append(read.get_tag("NM"))
        editDistancesMD.append(
            (
                sum([len(item) for item in re.split("[0-9^]", read.get_tag("MD"))])
                + sum(  # Parse MD string to get mismatches/deletions
                    [item[1] for item in read.cigartuples if item[0] == 1]
                )
            )
        )
        readLength.append(read.query_length)
        readAlignedLength.append(read.query_alignment_length)
    if nReads > 1:
        # get bases covered by reads pileup
        covPos = [
            pileupcolumn.n
            for pileupcolumn in samfile.pileup(
                chromosome, start=None, stop=None, region=None, stepper="nofilter"
            )
        ]

        basesCovered = int(len(covPos))
        # get SD from covered bases
        covSd = np.std(covPos)
        # get average coverage
        meanCoverage = sum(covPos) / refLen
        meanCoverageCovered = sum(covPos) / basesCovered

        breadth = basesCovered / refLen
        expBreadth = 1 - np.exp(-meanCoverage)
        breadthExpRatio = breadth / expBreadth

        covEvenness = coverage_evenness(covPos)

        cV = covSd / meanCoverage

        return (
            chromosome,
            int(nReads),
            np.mean(readLength),
            np.mean(readAlignedLength),
            np.mean(editDistancesNM),
            np.mean(editDistancesMD),
            np.mean(aniNM),
            np.mean(aniMD),
            basesCovered,
            meanCoverage,
            meanCoverageCovered,
            refLen,
            breadth,
            expBreadth,
            breadthExpRatio,
            cV,
            covEvenness,
        )
    else:
        return (
            chromosome,
            nReads,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            refLen,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
            np.nan,
        )


def processBam(bam, threads=1):
    """
    Processing function: calls pool of worker functions
    to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    """
    samfile = pysam.AlignmentFile(bam, "rb")
    if not samfile.has_index():
        pysam.index(bam)
        samfile = pysam.AlignmentFile(
            bam, "rb"
        )  # Need to reload the samfile after creating index
    chromosomes = samfile.references[0:100]
    datadf = pd.DataFrame()
    params = zip([bam] * len(chromosomes), chromosomes)
    try:
        print("Hello")
        p = Pool(threads)

        for result in tqdm.tqdm(
            p.imap_unordered(extractFromBam, params),
            total=len(chromosomes),
            leave=False,
            ncols=80,
            desc=f"Genomes processed",
        ):
            (
                chromosome,
                nReads,
                readLength,
                readAlignedLength,
                editDistancesNM,
                editDistancesMD,
                aniNM,
                aniMD,
                basesCovered,
                meanCoverage,
                meanCoverageCovered,
                refLen,
                breadth,
                expBreadth,
                breadthExpRatio,
                cV,
                covEvenness,
            ) = result
            datadf = datadf.append(
                {
                    "chromosome": chromosome,
                    "nReads": nReads,
                    "readLength": readLength,
                    "readAlignedLength": readAlignedLength,
                    "editDistancesNM": editDistancesNM,
                    "editDistancesMD": editDistancesMD,
                    "aniNM": aniNM,
                    "aniMD": aniMD,
                    "basesCovered": basesCovered,
                    "meanCoverage": meanCoverage,
                    "meanCoverageCovered": meanCoverageCovered,
                    "referenceLength": refLen,
                    "breadth": breadth,
                    "expBreadth": expBreadth,
                    "breadthExpRatio": breadthExpRatio,
                    "cV": cV,
                    "covEvenness": covEvenness,
                },
                ignore_index=True,
                sort=False,
            )

        p.close()
        p.join()

    except KeyboardInterrupt:
        print("Terminating worker threads")
        p.terminate()
        p.join()
        sys.exit()
    return datadf.convert_dtypes()


def filterReferenceBAMfile(bam, refs_dict, outBAMfile=None):
    (refNames, refLengths) = zip(*refs_dict.items())
    outBAMfile = pysam.Samfile(
        outBAMfile,
        "wb",
        referencenames=list(refNames),
        referencelengths=list(refLengths),
    )
    samfile = pysam.Samfile(bam, "rb")
    for read in samfile.fetch(multiple_iterators=True):
        if samfile.getrname(read.rname) in refNames:
            read.rname = refNames.index(samfile.getrname(read.rname))
            read.mrnm = read.rname
            outBAMfile.write(read)
    samfile.close()
    outBAMfile.close()

