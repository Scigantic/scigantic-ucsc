"""Read UCSC Genome Browser data from its public S3 bucket without downloading it."""
from .assemblies import Assembly, assemblies, find
from .conservation import TRACKS, chromosomes, conservation, conservation_bins
from .toga import Gene, toga_chromosomes, toga_genes
from .maf447 import alignment as multiz447_alignment, blocks as multiz447_blocks
from .multiz import Block, Row, alignment as multiz_alignment, blocks as multiz_blocks

__version__ = "0.3.0"

__all__ = [
    "Assembly", "Block", "Gene", "Row", "TRACKS",
    "assemblies", "multiz_alignment", "multiz447_alignment", "multiz447_blocks", "toga_chromosomes", "toga_genes", "chromosomes", "conservation", "conservation_bins", "find", "multiz_blocks",
]
