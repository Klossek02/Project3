import vcfpy
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import csv
import os
import warnings
warnings.filterwarnings("ignore")

# input SV files 
vcf_files = {
    'Delly': 'project3/delly.vcf',
    'BreakDancer': 'project3/breakdancer.vcf',
    'Pindel': 'project3/pindel.vcf'
}

TOLERANCE = 200
MIN_SV_LENGTH = 50

os.makedirs("project3/sv_results", exist_ok=True)

def get_sv_intervals(vcf_path):
    sv_list, lengths = [], []
    with vcfpy.Reader.from_path(vcf_path) as r:
        for rec in r:
            chrom = rec.CHROM
            svtype = rec.INFO.get('SVTYPE', 'NA')
            start = int(rec.POS)
            end   = int(rec.INFO.get('END', start))
            length = abs(end - start)
            if length >= MIN_SV_LENGTH:
                sv_list.append((chrom, svtype, start, end))
                lengths.append(length)
    return sv_list, lengths

# === UNIT TESTS ===
def test_get_sv_intervals():
    print("\n== TEST get_sv_intervals() ==")
    test_cases = {
        "project3/test_basic.vcf":        (3, [100,200,300]),
        "project3/test_min_length.vcf":   (2, [120,300]),      
        "project3/test_missing_svtype.vcf":(3, [50,100,150]),   
        "project3/test_multiple_chrom.vcf":(4, [70,250,50,100])
    }
    for path, (exp_n, exp_lens) in test_cases.items():
        sv_list, lengths = get_sv_intervals(path)
        assert len(sv_list) == exp_n, f"{path}: expected {exp_n} SV, is {len(sv_list)}"
        assert lengths == exp_lens,\
            f"{path}: expected length {exp_lens}, is {lengths}"
        print(f"  {os.path.basename(path)} OK")
    print("test_get_sv_intervals finished successfully!\n")

# === MAIN CODE ===
if __name__ == "__main__":
    test_get_sv_intervals()

    # collect SV from all tools
    sv_dict, length_dict = {}, {}
    for tool, fn in vcf_files.items():
        sl, ln = get_sv_intervals(fn)
        sv_dict[tool]    = sl
        length_dict[tool]= ln

    # fuzzy‐match
    def match_sv(sv, lst, tol=TOLERANCE):
        c1,t1,s1,e1 = sv
        for c2,t2,s2,e2 in lst:
            if c1==c2 and t1==t2 and abs(s1-s2)<=tol and abs(e1-e2)<=tol:
                return True
        return False

    only = lambda A,B,C: [sv for sv in sv_dict[A]
                         if not match_sv(sv,sv_dict[B]) and not match_sv(sv,sv_dict[C])]

    only_delly = only('Delly','BreakDancer','Pindel')
    only_breakdancer = only('BreakDancer','Delly','Pindel')
    only_pindel = only('Pindel','Delly','BreakDancer')

    DB = [sv for sv in sv_dict['Delly'] if match_sv(sv,sv_dict['BreakDancer']) and not match_sv(sv,sv_dict['Pindel'])]
    DP = [sv for sv in sv_dict['Delly'] if match_sv(sv,sv_dict['Pindel']) and not match_sv(sv,sv_dict['BreakDancer'])]
    BP = [sv for sv in sv_dict['BreakDancer'] if match_sv(sv,sv_dict['Pindel']) and not match_sv(sv,sv_dict['Delly'])]
    ALL= [sv for sv in sv_dict['Delly'] if match_sv(sv,sv_dict['BreakDancer']) and match_sv(sv,sv_dict['Pindel'])]

    
    plt.figure(figsize=(10,6))
    for t in vcf_files:
        plt.hist(length_dict[t], bins=40, alpha=0.5, label=t)
    plt.legend(); plt.xlabel('Length [bp]'); plt.ylabel('Number'); 
    plt.title('Histogram of SV length'); plt.tight_layout()
    plt.savefig('project3/sv_results/histogram_sv_lengths.png'); plt.close()

    venn3([
        set(map(str, sv_dict['Delly'])),
        set(map(str, sv_dict['BreakDancer'])),
        set(map(str, sv_dict['Pindel']))
    ], set_labels=('Delly','BreakDancer','Pindel'))
    plt.title('SV comparison'); plt.tight_layout()
    plt.savefig('project3/sv_results/sv_venn_fuzzy.png'); plt.close()


    def to_csv(lst, fn):
        with open(fn,'w',newline='') as f:
            w=csv.writer(f); w.writerow(['Chr','Typ','Start','End'])
            w.writerows(lst)
    to_csv(only_delly, 'project3/sv_results/only_delly.csv')
    to_csv(only_breakdancer, 'project3/sv_results/only_breakdancer.csv')
    to_csv(only_pindel, 'project3/sv_results/only_pindel.csv')
    to_csv(DB,'project3/sv_results/delly_breakdancer.csv')
    to_csv(DP, 'project3/sv_results/delly_pindel.csv')
    to_csv(BP, 'project3/sv_results/breakdancer_pindel.csv')
    to_csv(ALL, 'project3/sv_results/common_all.csv')

    print("\n==== SV SUMMARY ====")
    print(f"only Delly: {len(only_delly)}")
    print(f"only BreakDancer: {len(only_breakdancer)}")
    print(f"only Pindel: {len(only_pindel)}")
    print(f"Delly∩BD: {len(DB)}")
    print(f"Delly∩Pindel: {len(DP)}")
    print(f"BD∩Pindel: {len(BP)}")
    print(f"all three: {len(ALL)}")
    print("Wyniki w project3/sv_results/")
