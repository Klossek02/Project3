import os
import csv
import pysam
import vcfpy
from collections import defaultdict
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

ROOT = os.path.expanduser("~/project3")
BAM_IN = os.path.join(ROOT, "SRR_final_sorted.bam")
VCF_OUT = os.path.join(ROOT, "vcf_output.vcf")
RESULT_DIR = os.path.join(ROOT, "sv_results")

OTHER_VCFS = {
    "Delly": os.path.join(ROOT, "delly.vcf"),
    "BreakDancer": os.path.join(ROOT, "breakdancer.vcf"),
    "Pindel": os.path.join(ROOT, "pindel.vcf")
}

TEST_DIR = os.path.join(ROOT, "tests")
TEST_VCFS = {
    "test_basic.vcf":            (3, [100, 200, 300]),
    "test_min_length.vcf":       (2, [120, 300]),
    "test_missing_svtype.vcf":   (3, [50, 100, 150]),
    "test_multiple_chrom.vcf":   (4, [70, 250, 50, 100])}

RESULT_DIR  = "project3/sv_results"
MIN_MAPQ       = 30
MIN_SUPPORT    = 3
CLUSTER_WINDOW = 500
MIN_SV_LENGTH  = 50
TOLERANCE      = 1000

os.makedirs(RESULT_DIR, exist_ok=True)

# ——— CREATING/ CALLING SVs FROM BAM ———
def call_sv_from_bam(bam_path):
    bam  = pysam.AlignmentFile(bam_path, "rb")
    disc = defaultdict(list)
    for read in bam.fetch():
        if read.mapping_quality < MIN_MAPQ:
            continue
        if (not read.is_proper_pair) and (not read.is_unmapped) and (not read.mate_is_unmapped):
            key = (
                read.reference_name,
                read.next_reference_name,
                read.is_reverse,
                read.mate_is_reverse
            )
            disc[key].append(read)
    bam.close()

    sv_calls = []
    for reads in disc.values():
        reads.sort(key=lambda r: r.reference_start)
        cluster = [reads[0]]
        for r in reads[1:]:
            if r.reference_start - cluster[-1].reference_start <= CLUSTER_WINDOW:
                cluster.append(r)
            else:
                if len(cluster) >= MIN_SUPPORT:
                    sv_calls.append(_cluster_to_sv(cluster))
                cluster = [r]
        if len(cluster) >= MIN_SUPPORT:
            sv_calls.append(_cluster_to_sv(cluster))
    return sv_calls

"""
def _cluster_to_sv(reads):
    chrom = reads[0].reference_name
    starts = [r.reference_start for r in reads]
    ends = [r.reference_end   for r in reads]
    start = min(starts) + 1
    end = max(ends)
    svtype = "DEL"   # simple deletion
    return (chrom, svtype, start, end)
"""

def _cluster_to_sv(reads):
    """
    Turning a cluster of discordant reads into one SV call, classifying
    by orientation:
      - different chromosomes --> BND (translocation)
      - same orientation (both F-->F or both R-->R) --> INV
      - F-->R orientation --> DEL
      - R-->F orientation --> DUP
    """

    chrom1 = reads[0].reference_name
    chrom2 = reads[0].next_reference_name

    starts = [r.reference_start for r in reads]
    ends   = [r.reference_end   for r in reads]
    start  = min(starts) + 1
    end    = max(ends)

    # classifying by chrom and orientation
    if chrom1 != chrom2:
        svtype = "BND"
    else:
        # True if the read itself is mapped forward
        fwd      = not reads[0].is_reverse
        # True if its mate is mapped forward
        mate_fwd = not reads[0].mate_is_reverse

        if fwd == mate_fwd:
            # both forward or both reverse --> inversion
            svtype = "INV"
        elif fwd and not mate_fwd:
            # F-->R --> deletion signature
            svtype = "DEL"
        elif not fwd and mate_fwd:
            # R-->F --> tandem duplication signature
            svtype = "DUP"
        else:
            # fallback
            svtype = "BND"

    return (chrom1, svtype, start, end)


def write_sv_to_vcf(sv_calls, vcf_path):
    header_lines = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]
    with open(vcf_path, 'w') as f:
        for h in header_lines:
            f.write(h + '\n')
        for chrom, svtype, start, end in sorted(sv_calls, key=lambda x: (x[0], x[2])):
            alt  = f"<{svtype}>"
            info = f"SVTYPE={svtype};END={end}"
            f.write(f"{chrom}\t{start}\t.\tN\t{alt}\t.\tPASS\t{info}\n")

# ——— READ BACK VCF INTO INTERVALS ———
def get_sv_intervals(vcf_path):
    sv_list, lengths = [], []
    with vcfpy.Reader.from_path(vcf_path) as r:
        for rec in r:
            chrom  = rec.CHROM
            svtype = rec.INFO.get('SVTYPE', 'NA')
            start  = int(rec.POS)
            end    = int(rec.INFO.get('END', start))
            length = abs(end - start)
            if length >= MIN_SV_LENGTH:
                sv_list.append((chrom, svtype, start, end))
                lengths.append(length)
    return sv_list, lengths

# ——— UNIT TESTS ———
def test_get_sv_intervals():
    print("\n== TEST get_sv_intervals() ==")
    for fname, (exp_n, exp_lens) in TEST_VCFS.items():
        path = os.path.join(TEST_DIR, fname)
        svs, lens = get_sv_intervals(path)
        assert len(svs) == exp_n, f"{fname}: expected {exp_n} SVs, got {len(svs)}"
        assert lens == exp_lens, f"{fname}: expected lengths {exp_lens}, got {lens}"
        print(f"  {fname} OK")
    print("test_get_sv_intervals finished successfully!\n")

# ——— COMPARE & PLOT ———
def compare_and_plot(output_vcf, other_vcfs):
    write_sv_to_vcf(call_sv_from_bam(BAM_IN), output_vcf)

    # loading every set
    all_sv, all_len = {}, {}
    for name, path in {**other_vcfs, 'Output': output_vcf}.items():
        sl, ln = get_sv_intervals(path)
        all_sv[name]   = sl
        all_len[name]  = ln

    # fuzzy match helper
    def match(sv, lst):
        c1,t1,s1,e1 = sv
        return any(c1==c2 and t1==t2 and abs(s1-s2)<=TOLERANCE and abs(e1-e2)<=TOLERANCE
                   for c2,t2,s2,e2 in lst)

    # uniques
    only_output = [sv for sv in all_sv['Output'] if not match(sv, all_sv['Delly']) and not match(sv, all_sv['BreakDancer']) and not match(sv, all_sv['Pindel'])]
    only_delly = [sv for sv in all_sv['Delly'] if not match(sv, all_sv['Output']) and not match(sv, all_sv['BreakDancer']) and not match(sv, all_sv['Pindel'])]
    only_breakdancer= [sv for sv in all_sv['BreakDancer'] if not match(sv, all_sv['Output']) and not match(sv, all_sv['Delly']) and not match(sv, all_sv['Pindel'])]
    only_pindel = [sv for sv in all_sv['Pindel'] if not match(sv, all_sv['Output']) and not match(sv, all_sv['Delly']) and not match(sv, all_sv['BreakDancer'])]

    # pairwise overlaps (only A ∩ B for now)
    OD = [sv for sv in all_sv['Output'] if match(sv, all_sv['Delly']) and not match(sv, all_sv['BreakDancer']) and not match(sv, all_sv['Pindel'])]
    OB = [sv for sv in all_sv['Output'] if match(sv, all_sv['BreakDancer']) and not match(sv, all_sv['Delly']) and not match(sv, all_sv['Pindel'])]
    OP = [sv for sv in all_sv['Output'] if match(sv, all_sv['Pindel']) and not match(sv, all_sv['Delly']) and not match(sv, all_sv['BreakDancer'])]
    DB = [sv for sv in all_sv['Delly'] if match(sv, all_sv['BreakDancer']) and not match(sv, all_sv['Output']) and not match(sv, all_sv['Pindel'])]
    DP = [sv for sv in all_sv['Delly'] if match(sv, all_sv['Pindel']) and not match(sv, all_sv['Output']) and not match(sv, all_sv['BreakDancer'])]
    BP = [sv for sv in all_sv['BreakDancer'] if match(sv, all_sv['Pindel']) and not match(sv, all_sv['Output']) and not match(sv, all_sv['Delly'])]

    # all-three
    ALL = [sv for sv in all_sv['Output'] if match(sv, all_sv['Delly']) and match(sv, all_sv['BreakDancer']) and match(sv, all_sv['Pindel'])]

    # histogram
    plt.figure(figsize=(8,5))
    for tool, lens in all_len.items():
        plt.hist(lens, bins=30, alpha=0.5, label=tool)
    plt.legend(); plt.title("SV length comparison")
    plt.savefig(f"{RESULT_DIR}/histogram.png"); plt.close()

    # Venn (Output / Delly / Pindel)
    venn3([
        set(map(str, all_sv['Output'])),
        set(map(str, all_sv['Delly'])),
        set(map(str, all_sv['BreakDancer']))
    ], set_labels=('Output','Delly','BreakDancer'))
    plt.title("SV overlap")
    plt.savefig(f"{RESULT_DIR}/venn.png"); plt.close()

    def to_csv(lst, fn):
        with open(fn,'w',newline='') as f:
            w = csv.writer(f)
            w.writerow(['Chr','Type','Start','End'])
            w.writerows(lst)

    mapping = {
        "only_output.csv": only_output,
        "only_delly.csv": only_delly,
        "only_breakdancer.csv": only_breakdancer,
        "only_pindel.csv": only_pindel,
        "output_delly.csv": OD,
        "output_breakdancer.csv": OB,
        "output_pindel.csv": OP,
        "delly_breakdancer.csv":DB,
        "delly_pindel.csv": DP,
        "breakdancer_pindel.csv":BP,
        "common_all.csv": ALL
    }
    for fn, lst in mapping.items():
        to_csv(lst, os.path.join(RESULT_DIR, fn))

    print("--- SV SUMMARY ---")
    print(f"Only output: {len(only_output)}")
    print(f"Only Delly: {len(only_delly)}")
    print(f"Only BreakDancer: {len(only_breakdancer)}")
    print(f"Only Pindel: {len(only_pindel)}")
    print(f"Output ∩ Delly: {len(OD)}")
    print(f"Output ∩ BD: {len(OB)}")
    print(f"Output ∩ Pindel: {len(OP)}")
    print(f"Delly ∩ BD: {len(DB)}")
    print(f"Delly ∩ Pindel: {len(DP)}")
    print(f"BD ∩ Pindel: {len(BP)}")
    print(f"All three: {len(ALL)}")
    print(f"Results in {RESULT_DIR}/")

if __name__ == "__main__":
    test_get_sv_intervals()
    compare_and_plot(VCF_OUT, OTHER_VCFS)