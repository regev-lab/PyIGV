import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyigv import Alignment, plot_alignments


def test_alignment_basic():
    """Test basic alignment creation"""
    target = "AAATAAA"
    query = "AAAGAAA"
    alignment = ["AAATAAA", "AAAGAAA"]

    aln = Alignment(alignment, target, query)

    # Check that the alignment was created
    assert aln.target == target
    assert aln.query == query
    assert aln.target_alignment == alignment[0]
    assert aln.query_alignment == alignment[1]


def test_alignment_with_gaps():
    """Test alignment with gaps (insertions/deletions)"""
    target = "AAATAAA"
    query = "AAAGAAA"
    alignment = ["AAA-TAAA", "AAAGAAA-"]

    aln = Alignment(alignment, target, query)

    assert aln.insertion_ct >= 0
    assert aln.deletion_ct >= 0
    assert aln.mutation_ct >= 0


def test_alignment_mutation_counting():
    """Test mutation counting"""
    target = "AAAA"
    query = "AAAT"
    alignment = ["AAAA", "AAAT"]

    aln = Alignment(alignment, target, query)

    # There should be 1 mismatch
    assert aln.mutation_ct == 1
    assert aln.insertion_ct == 0
    assert aln.deletion_ct == 0

def test_alignment_multi_mutation_counting():
    """Test mutation counting"""
    target = "AAACCCAAA"
    query = "AAATTTAAA"
    alignment = ["AAA---CCCAAA", "AAATTT---AAA"]

    aln = Alignment(alignment, target, query)

    # There should be 1 mismatch
    assert aln.mutation_ct == 3
    assert aln.insertion_ct == 0
    assert aln.deletion_ct == 0


def test_alignment_insertion():
    """Test insertion detection"""
    target = "AAAA"
    query = "AAAAA"
    alignment = ["AAAA-", "AAAAA"]

    aln = Alignment(alignment, target, query)

    assert aln.insertion_ct == 1


def test_alignment_deletion():
    """Test deletion detection"""
    target = "AAAAA"
    query = "AAAA"
    alignment = ["AAAAA", "AAAA-"]

    aln = Alignment(alignment, target, query)

    assert aln.deletion_ct == 1


def test_alignment_comparison():
    """Test alignment comparison (less than operator)"""
    target = "AAAA"
    query1 = "AAAA"
    query2 = "AAAT"

    aln1 = Alignment(["AAAA", "AAAA"], target, query1)
    aln2 = Alignment(["AAAA", "AAAT"], target, query2)

    # aln1 has fewer mutations, so it should be "less than" aln2
    assert aln1 < aln2


def test_alignment_str():
    """Test string representation"""
    target = "AAAA"
    query = "AAAT"
    alignment = ["AAAA", "AAAT"]

    aln = Alignment(alignment, target, query)
    str_repr = str(aln)

    assert "Target:" in str_repr
    assert "Query:" in str_repr
    assert "Edits:" in str_repr


def test_get_color_row():
    """Test color row generation"""
    target = "AAAA"
    query = "AAAT"
    alignment = ["AAAA", "AAAT"]

    aln = Alignment(alignment, target, query)
    color_row = aln.get_color_row()

    assert isinstance(color_row, list)
    assert len(color_row) > 0


def test_get_symbols():
    """Test symbol extraction"""
    target = "AAAA"
    query = "AAAT"
    alignment = ["AAAA", "AAAT"]

    aln = Alignment(alignment, target, query)
    symbols = aln.get_symbols()

    assert isinstance(symbols, list)
    assert len(symbols) > 0


def test_get_insertion_indices():
    """Test insertion indices extraction"""
    target = "AAAA"
    query = "AAAAA"
    alignment = ["AAAA-", "AAAAA"]

    aln = Alignment(alignment, target, query)
    indices = aln.get_insertion_indices()

    assert isinstance(indices, list)


def test_plot_alignments_basic():
    """Test basic plotting functionality"""
    target = "AAATAAA"
    query1 = "AAAGAAA"
    query2 = "AAACAAA"

    aln1 = Alignment(["AAATAAA", "AAAGAAA"], target, query1)
    aln2 = Alignment(["AAATAAA", "AAACAAA"], target, query2)

    alignments = [aln1, aln2]

    # Test plotting without saving
    fig = plot_alignments(alignments, barcode="TEST", overlap_name="test_overlap")

    assert fig is not None
    plt.close(fig)


def test_plot_alignments_with_pdf(tmp_path):
    """Test plotting with PDF output"""
    output_path = tmp_path / "test_plot.pdf"

    target = "AAATAAA"
    query1 = "AAAGAAA"

    aln1 = Alignment(["AAATAAA", "AAAGAAA"], target, query1)
    alignments = [aln1]

    # Create PDF and plot
    with PdfPages(str(output_path)) as pdf:
        fig = plot_alignments(alignments, pdf=pdf)
        plt.close(fig)

    # Check that the file exists
    assert output_path.exists()
    assert output_path.stat().st_size > 0


def test_plot_alignments_truncated():
    """Test truncated plotting mode"""
    target = "AAATAAA"
    query1 = "AAAGAAA"

    aln1 = Alignment(["AAATAAA", "AAAGAAA"], target, query1)
    alignments = [aln1]

    # Test truncated mode
    fig = plot_alignments(alignments, truncate=True)

    assert fig is not None
    plt.close(fig)
