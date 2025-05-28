import os
import csv
import pysam
import vcfpy
from collections import defaultdict
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import warnings
import pandas as pd

from plots import plot_len_by_type, plot_jaccard_heat, plot_svtype_bar, plot_upset

# ——— GLOBAL SETTINGS ——————————————————————————————

warnings.filterwarnings("ignore")

ROOT = os.path.expanduser(".")
BAM_IN = os.path.join(ROOT, "data/SRR_final_sorted.bam")
VCF_OUT = os.path.join(ROOT, "vcf_output.vcf")
RESULT_DIR = os.path.join(ROOT, "results")

OTHER_VCFS = {
    "Delly": os.path.join(ROOT, "delly.vcf"),
    "BreakDancer": os.path.join(ROOT, "breakdancer.vcf"),
    "Pindel": os.path.join(ROOT, "pindel.vcf")
}

TEST_DIR = os.path.join(ROOT, "data")
TEST_VCFS = {
    "test_basic.vcf": (3, [100, 200, 300]),
    "test_min_length.vcf": (2, [120, 300]),
    "test_missing_svtype.vcf": (3, [50, 100, 150]),
    "test_multiple_chrom.vcf": (4, [70, 250, 50, 100])}

MIN_MAPQ = 30
MIN_SUPPORT = 3
CLUSTER_WINDOW = 500
MIN_SV_LENGTH = 50
TOLERANCE = 1000

os.makedirs(RESULT_DIR, exist_ok=True)

# ——— CREATING/ CALLING SVs FROM BAM ———
def call_sv_from_bam(bam_path):
    bam  = pysam.AlignmentFile(bam_path, "rb")
    disc = defaultdict(list)
    for read in bam.fetch():
        if read.mapping_quality < MIN_MAPQ:
            continue
        if (not read.is_proper_pair) and (not read.is_unmapped) and (not read.mate_is_unmapped):
            key = (read.reference_name, read.next_reference_name, read.is_reverse, read.mate_is_reverse)
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
    ends = [r.reference_end   for r in reads]
    start = min(starts) + 1
    end = max(ends)

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
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]

    with open(vcf_path, 'w') as f:
        for h in header_lines:
            f.write(h + '\n')
        for idx, (chrom, svtype, start, end) in enumerate(sorted(sv_calls, key=lambda x: (x[0], x[2]))):
            alt = f"<{svtype}>"
            svlen = end - start + 1
            info = f"SVTYPE={svtype};END={end};SVLEN={svlen}"
            id_ = f"SV{idx:05d}"  # e.g. SV00001 and so on
            f.write(f"{chrom}\t{start}\t{id_}\tN\t{alt}\t.\tPASS\t{info}\n")


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
def compare_and_plot(output_vcf, other_vcfs, result_dir=RESULT_DIR):
    write_sv_to_vcf(call_sv_from_bam(BAM_IN), output_vcf)

    # loading every set
    all_sv, all_len = {}, {}
    for name, path in {**other_vcfs, 'Output': output_vcf}.items():
        sl, ln = get_sv_intervals(path)
        all_sv[name] = sl
        all_len[name] = ln

    overlap_cnt, out_cnt, overlap = overlap_score(
        all_sv['Output'],
        all_sv['Delly'],
        all_sv['BreakDancer'],
        all_sv['Pindel'],
        TOLERANCE)

    print("\n=== OVERLAP METRIC ===")
    print(f"Output variants: {out_cnt}")
    print(f"Overlapping with others: {overlap_cnt}")
    print(f"OVERLAP score: {overlap:.3f}\n")

    overlap_df = pd.DataFrame({
        'overlap_count': [overlap_cnt],
        'output_count': [out_cnt],
        'overlap_score': [overlap]
    })
    
    # Save to CSV
    overlap_df.to_csv(f"{result_dir}/overlap_scores.csv", index=False)



    plot_len_by_type(all_sv, result_dir)
    plot_jaccard_heat(all_sv, TOLERANCE, result_dir)
    plot_svtype_bar(all_sv, result_dir)

    try:
        plot_upset(all_sv, TOLERANCE, result_dir)
    except ImportError:
        print("matplotlib-upset not installed – skipping UpSet plot")

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
    plt.legend()
    plt.title("SV length comparison")
    plt.savefig(f"{result_dir}/histogram.png")
    plt.close()

    # Venn (Output / Delly / Pindel)
    venn3([
        set(map(str, all_sv['Output'])),
        set(map(str, all_sv['Delly'])),
        set(map(str, all_sv['BreakDancer']))
    ], set_labels=('Output','Delly','BreakDancer'))
    plt.title("SV overlap")
    plt.savefig(f"{result_dir}/venn.png")
    plt.close()

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
        "common_all.csv": ALL}

    for fn, lst in mapping.items():
        to_csv(lst, os.path.join(result_dir, fn))

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
    print(f"Results in {result_dir}/")


def overlap_score(output_set, delly_set, bd_set, pindel_set, tol):
    union_ref = delly_set + bd_set + pindel_set

    inter = [
        sv for sv in output_set
        if any(
            # we use the same 'match'
            (sv[0] == r[0] and # chromosome
             sv[1] == r[1] and # type
             abs(sv[2]-r[2]) <= tol and abs(sv[3]-r[3]) <= tol)
            for r in union_ref)]

    overlap_cnt = len(inter)
    output_cnt  = len(output_set) or 1  
    score       = overlap_cnt / output_cnt
    return overlap_cnt, output_cnt, score


def fake_read(rname, pos, is_rev, mrname, mpos, mate_rev):
    header = pysam.AlignmentHeader.from_dict({
        "SQ": [{"SN": "chr1", "LN": 1_000_000},
               {"SN": "chr2", "LN": 1_000_000}]})

    a = pysam.AlignedSegment(header)
    a.reference_id = 0 if rname == "chr1" else 1
    a.reference_start = pos
    a.cigarstring = "100M"
    a.is_reverse = is_rev
    a.next_reference_id = 0 if mrname == "chr1" else 1
    a.next_reference_start= mpos
    a.mate_is_reverse = mate_rev
    a.mapping_quality = 60
    return a


def run_with_params(tol, min_sup, min_sv_len):
    global TOLERANCE, MIN_SUPPORT, MIN_SV_LENGTH
    TOLERANCE = tol
    MIN_SUPPORT = min_sup
    MIN_SV_LENGTH = min_sv_len
    
    compare_and_plot_with_parametrized_dir(VCF_OUT, OTHER_VCFS)

def compare_and_plot_with_parametrized_dir(output_vcf, other_vcfs):
    dir = f"{RESULT_DIR}/tol={TOLERANCE},min_sup={MIN_SUPPORT},min_sv_len={MIN_SV_LENGTH}"
    os.makedirs(dir, exist_ok=True)
    compare_and_plot(f"{dir}/{output_vcf}", other_vcfs, result_dir=dir)



if __name__ == "__main__":
    test_get_sv_intervals()
    for tol in [200, 500, 1000]:
        for sup in [2, 3, 4]:
            for svl in [50, 100]:
                print(f"\n### Testing T = {tol}, SUP = {sup}, MIN_LEN = {svl}")
                run_with_params(tol, sup, svl)

    # Best for SUP=3, MINLEN=50, TOLERANCE=1000 - OVERLAP = 0.044
    # But SUP=4, MINLEN=100, TOLERANCE=500–1000 with 20% overlap is also good
