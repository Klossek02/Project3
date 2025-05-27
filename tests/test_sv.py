import sys
import os
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


from project3 import _cluster_to_sv, fake_read, compare_and_plot, call_sv_from_bam, BAM_IN
import pysam

def test_inv_called_for_same_orientation():
    cluster = [
        fake_read("chr1", 100, False, "chr1", 1000, False),
        fake_read("chr1", 150, False, "chr1", 1050, False),
        fake_read("chr1", 200, False, "chr1", 1100, False)
    ]
    chrom, svtype, start, end = _cluster_to_sv(cluster)
    assert svtype == "INV"

def test_compare_creates_outputs(tmp_path):
    delly = tmp_path / "delly_stub.vcf"
    bd    = tmp_path / "bd_stub.vcf"
    pindel = tmp_path / "pindel_stub.vcf"

    for f in [delly, bd, pindel]:
        f.write_text("""##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n""")

    output_vcf = tmp_path / "output.vcf"
    result_dir = tmp_path / "sv_results"
    os.makedirs(result_dir, exist_ok=True)

    compare_and_plot(output_vcf, {
        "Delly": delly,
        "BreakDancer": bd,
        "Pindel": pindel
    },  result_dir=result_dir)

    assert (result_dir / "venn.png").exists()
    assert sum(1 for _ in result_dir.glob("*.csv")) == 11

@pytest.mark.slow  # Mark this test as slow (for selective running)
@pytest.mark.skipif(
    not os.path.exists(BAM_IN),
    reason="SRR_final_sorted.bam not found – skipping integration test",
)
def test_human_bam():
    """
    Integration test using human sequencing data.

    Under the default parameters (e.g., MIN_SUPPORT=3, CLUSTER_WINDOW=500),
    the detector should find exactly 90 structural variants (SVs).
    If parameters are changed in the future, update the expected count here.
    """
    sv_calls = call_sv_from_bam(BAM_IN)
    assert len(sv_calls) == 90, f"Expected 90 SVs, but got {len(sv_calls)}"
