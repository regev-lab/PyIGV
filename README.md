# PyIGV

Python alignment viewer library based on the Integrative Genomics Viewer (IGV) style for visualizing DNA/RNA sequence alignments.

## Overview

PyIGV provides a simple, intuitive way to visualize pairwise sequence alignments in Python. It displays alignments in an IGV-like format, with color-coded mismatches, insertions, and deletions.

## Installation

```bash
pip install pyigv
```

## Features

- **Color-coded visualization**: Mismatches are highlighted with base-specific colors (A=green, T=red, G=gold, C=blue)
- **Gap handling**: Automatically detects and visualizes insertions and deletions
- **Mutation counting**: Tracks the number of insertions, deletions, and substitutions
- **PDF export**: Save alignment visualizations to PDF files
- **Flexible display**: Option to show full alignments or truncated views

## Quick Start

```python
from pyigv import Alignment, plot_alignments

# Define your sequences
target = "AAATAAA"
query = "AAAGAAA"

# Create alignment (with gap notation using '-')
alignment = ["AAATAAA", "AAAGAAA"]

# Create Alignment object
aln = Alignment(alignment, target, query)

# Print alignment information
print(aln)
print(f"Mutations: {aln.mutation_ct}")
print(f"Insertions: {aln.insertion_ct}")
print(f"Deletions: {aln.deletion_ct}")

# Visualize multiple alignments
alignments = [aln]
plot_alignments(alignments, barcode="Sample1", overlap_name="Region1")
```

## Usage Examples

### Basic Alignment

```python
from pyigv import Alignment

# Perfect match
target = "AAAA"
query = "AAAA"
alignment = ["AAAA", "AAAA"]
aln = Alignment(alignment, target, query)
print(f"Mutations: {aln.mutation_ct}")  # Output: 0
```

### Alignment with Mismatch

```python
# Single mismatch at position 3
target = "AAAA"
query = "AAAT"
alignment = ["AAAA", "AAAT"]
aln = Alignment(alignment, target, query)
print(f"Mutations: {aln.mutation_ct}")  # Output: 1
```

### Alignment with Insertion

```python
# Insertion in query
target = "AAAA"
query = "AAAAA"
alignment = ["AAAA-", "AAAAA"]  # '-' indicates gap in target
aln = Alignment(alignment, target, query)
print(f"Insertions: {aln.insertion_ct}")  # Output: 1
```

### Alignment with Deletion

```python
# Deletion in query
target = "AAAAA"
query = "AAAA"
alignment = ["AAAAA", "AAAA-"]  # '-' indicates gap in query
aln = Alignment(alignment, target, query)
print(f"Deletions: {aln.deletion_ct}")  # Output: 1
```

### Plotting Multiple Alignments

```python
from pyigv import Alignment, plot_alignments
import matplotlib.pyplot as plt

target = "AAATAAA"

# Create multiple alignments
aln1 = Alignment(["AAATAAA", "AAAGAAA"], target, "AAAGAAA")
aln2 = Alignment(["AAATAAA", "AAACAAA"], target, "AAACAAA")
aln3 = Alignment(["AAATAAA", "AAAAAAA"], target, "AAAAAAA")

alignments = [aln1, aln2, aln3]

# Plot and display
plot_alignments(
    alignments,
    barcode="TEST001",
    overlap_name="Region_1"
)
plt.show()
```

### Saving to PDF

```python
from pyigv import plot_alignments
from matplotlib.backends.backend_pdf import PdfPages

# Create your alignments
alignments = [aln1, aln2, aln3]

# Save to PDF
with PdfPages("alignment_output.pdf") as pdf:
    plot_alignments(
        alignments,
        barcode="TEST001",
        overlap_name="Region_1",
        pdf=pdf
    )
```

### Truncated View

For alignments with many insertions, you can use truncated view to focus on the reference sequence:

```python
plot_alignments(
    alignments,
    barcode="TEST001",
    overlap_name="Region_1",
    truncate=True  # Removes insertions from display
)
```

## API Reference

### `Alignment` Class

#### Constructor

```python
Alignment(alignment: Sequence[str], target: str, query: str)
```

**Parameters:**
- `alignment`: A list/tuple of two strings representing the aligned sequences with gaps marked as '-'
- `target`: The target (reference) sequence
- `query`: The query sequence

#### Attributes

- `target`: Target sequence
- `query`: Query sequence
- `target_alignment`: Aligned target sequence with gaps
- `query_alignment`: Aligned query sequence with gaps
- `symbols`: Processed alignment symbols
- `edits`: Edit operations (I=insertion, D=deletion, M=mismatch, space=match)
- `insertion_ct`: Number of insertions
- `deletion_ct`: Number of deletions
- `mutation_ct`: Number of mismatches/substitutions

#### Methods

- `get_color_row(truncate=False)`: Get color codes for visualization
- `get_symbols(truncate=False)`: Get alignment symbols
- `get_insertion_indices()`: Get positions and lengths of insertions

### `plot_alignments` Function

```python
plot_alignments(alignments, barcode="", overlap_name="", pdf=None, truncate=False)
```

**Parameters:**
- `alignments`: List of Alignment objects to visualize
- `barcode`: Label for the barcode/sample
- `overlap_name`: Name of the overlapping region
- `pdf`: PdfPages object for saving to PDF (optional)
- `truncate`: If True, removes insertions from display

**Returns:**
- matplotlib Figure object

## Color Scheme

- **Green (A)**: Adenine mismatches
- **Red (T)**: Thymine mismatches
- **Gold (G)**: Guanine mismatches
- **Blue (C)**: Cytosine mismatches
- **Gray**: Matches
- **White**: Deletions

## Development

### Running Tests

```bash
# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest tests/
```

### Code Quality

```bash
# Format code
black src/ tests/

# Lint code
flake8 src/ tests/
```

## License

MIT License - see LICENSE file for details

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Citation

If you use PyIGV in your research, please cite:

```
PyIGV: Python alignment viewer library
https://github.com/regev-lab/PyIGV
```

## Support

For issues, questions, or contributions, please visit:
- Issue Tracker: https://github.com/regev-lab/PyIGV/issues
- Source Code: https://github.com/regev-lab/PyIGV
