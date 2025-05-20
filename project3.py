import vcfpy
from matplotlib_venn import venn3
import matplotlib.pyplot as plt
import csv
import os

vcf_files = {'Delly': 'project3/delly.vcf', 'BreakDancer': 'project3/breakdancer.vcf', 'Pindel': 'project3/pindel.vcf'}

TOLERANCE = 200    # position tolerance (bp) -- can be adjusted 
MIN_SV_LENGTH = 50 # minimal SV length for comparison 

os.makedirs("project3/sv_results", exist_ok = True)

def get_sv_intervals(vcf_path):

    sv_list = []
    lengths = []

    with vcfpy.Reader.from_path(vcf_path) as reader:
        for record in reader:
            chrom = record.CHROM
            svtype = record.INFO.get('SVTYPE', 'NA')
            start = int(record.POS)
            end = int(record.INFO.get('END', start))
            length = abs(end - start)
            if length >= MIN_SV_LENGTH:
                sv_list.append((chrom, svtype, start, end))
                lengths.append(length)
    return sv_list, lengths

# SV lists + their lengths
sv_dict = {}
length_dict = {}
for tool, path in vcf_files.items():
    sv_list, lengths = get_sv_intervals(path)
    sv_dict[tool] = sv_list
    length_dict[tool] = lengths

# fuzzy match (for better tolerance )
def match_sv(sv, sv_list, tolerance = TOLERANCE):
    chrom1, svtype1, start1, end1 = sv
    for chrom2, svtype2, start2, end2 in sv_list:
        if chrom1 == chrom2 and svtype1 == svtype2:
            if abs(start1 - start2) <= tolerance and abs(end1 - end2) <= tolerance:
                return True
    return False

# unique and common pairs 
only_delly = [sv for sv in sv_dict['Delly']
              if not match_sv(sv, sv_dict['BreakDancer']) and not match_sv(sv, sv_dict['Pindel'])]

only_breakdancer = [sv for sv in sv_dict['BreakDancer']
                    if not match_sv(sv, sv_dict['Delly']) and not match_sv(sv, sv_dict['Pindel'])]

only_pindel = [sv for sv in sv_dict['Pindel']
               if not match_sv(sv, sv_dict['Delly']) and not match_sv(sv, sv_dict['BreakDancer'])]

# pairs 
delly_breakdancer = [sv for sv in sv_dict['Delly'] if match_sv(sv, sv_dict['BreakDancer']) and not match_sv(sv, sv_dict['Pindel'])]
delly_pindel = [sv for sv in sv_dict['Delly'] if match_sv(sv, sv_dict['Pindel']) and not match_sv(sv, sv_dict['BreakDancer'])]
breakdancer_pindel = [sv for sv in sv_dict['BreakDancer'] if match_sv(sv, sv_dict['Pindel']) and not match_sv(sv, sv_dict['Delly'])]

# common for all 
common_all = [sv for sv in sv_dict['Delly'] if match_sv(sv, sv_dict['BreakDancer']) and match_sv(sv, sv_dict['Pindel'])]


plt.figure(figsize=(10,6))
for tool in vcf_files:
    plt.hist(length_dict[tool], bins=40, alpha=0.5, label=tool)
plt.xlabel('SV length [bp]')
plt.ylabel('No. of variants')
plt.title('Histogram of SV length for each tool')
plt.legend()
plt.tight_layout()
plt.savefig('project3/sv_results/histogram_sv_lengths.png')
plt.show()

# Venn diagram
venn3([
    set([str(sv) for sv in sv_dict['Delly']]),
    set([str(sv) for sv in sv_dict['BreakDancer']]),
    set([str(sv) for sv in sv_dict['Pindel']])
    ],
    set_labels=('Delly', 'BreakDancer', 'Pindel'))
plt.title('SV comparision')
plt.tight_layout()
plt.savefig('project3/sv_results/sv_venn_fuzzy.png')
plt.show()


def write_sv_to_csv(sv_list, filename):
    with open(filename, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Chr', 'Typ', 'Start', 'End'])
        for sv in sv_list:
            w.writerow(sv)

write_sv_to_csv(only_delly, 'project3/sv_results/only_delly.csv')
write_sv_to_csv(only_breakdancer, 'project3/sv_results/only_breakdancer.csv')
write_sv_to_csv(only_pindel, 'project3/sv_results/only_pindel.csv')
write_sv_to_csv(delly_breakdancer, 'project3/sv_results/delly_breakdancer.csv')
write_sv_to_csv(delly_pindel, 'project3/sv_results/delly_pindel.csv')
write_sv_to_csv(breakdancer_pindel, 'project3/sv_results/breakdancer_pindel.csv')
write_sv_to_csv(common_all, 'project3/sv_results/common_all.csv')

print("\ SV SUMMARY ")
print(f"SV no. only in Delly: {len(only_delly)}")
print(f"SV no. only in BreakDancer: {len(only_breakdancer)}")
print(f"SV no. only in Pindel: {len(only_pindel)}")
print(f"SV common no. for Delly and BreakDancer: {len(delly_breakdancer)}")
print(f"SV common no. for Delly and Pindel: {len(delly_pindel)}")
print(f"SV common no. for BreakDancer and Pindel: {len(breakdancer_pindel)}")
print(f"SV common no. for all 3 tools: {len(common_all)}")

