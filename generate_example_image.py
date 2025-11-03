#!/usr/bin/env python3
"""Generate example image for README"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from src.pyigv import Alignment, plot_alignments

# Create example alignments
target = "AAACCCGGGTTTATATATAT"
queries = [
    "AAACCCGGGTTTATATATAT",  # Perfect match
    "AAAGCCGGGTTTATATATAT",  # One mismatch (C->G)
    "AAACCCGGGTTTTATATAT",   # One deletion
    "AAACCCGGGTTTATATATATAT", # One insertion
    "AAATTTGGGAAACCCCCCCC",  # Multiple changes
]

# Create alignments
alignments = [Alignment(target, query) for query in queries]

# Generate plot
fig = plot_alignments(alignments, title="Example Alignment Visualization", return_fig=True)

# Save as PNG with high DPI for clarity
plt.savefig('docs/images/example_output.png', dpi=150, bbox_inches='tight')
print("Example image saved to docs/images/example_output.png")

# Also generate a truncated view example
fig2 = plot_alignments(alignments, title="Example Truncated View", truncate=True, return_fig=True)
plt.savefig('docs/images/example_output_truncated.png', dpi=150, bbox_inches='tight')
print("Truncated example image saved to docs/images/example_output_truncated.png")

plt.close('all')
