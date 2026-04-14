#!/usr/bin/env python3
#Author: Dede Kurniawan
#Email: dedekurniawan@genomics.cn or dedearkun2710@gmail.com

"""
Compute and visualize GC content vs sequencing depth per genomic window. One of robust way to identify contamination in the genome.

Usage: gc-depth-plot.py <fasta> <pandepth_output> [options]

Positional arguments:
  fasta                Genome FASTA file (gzipped is also fine)
  pandepth             Pandepth windowed depth file (.win.stat.gz)

Options:
  -h, --help           Show this help message and exit
  -w, --window WINDOW  Window size, must match pandepth -w value (default: 1000)
  -o, --output OUTPUT  Output plot file (.png or .pdf, default: gc-depth.png)
  --log-depth          Use logarithmic scale for the depth axis
  --plot-only TSV      Skip processing, re-plot from an existing combined TSV (from --output-data)
  --output-data FILE   Save merged GC and depth data to this TSV file (can be reused with --plot-only)

Author: Dede Kurniawan
Email: dedekurniawan@genomics.cn or dedearkun2710@gmail.com
"""

import argparse
import gzip
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from scipy.ndimage import gaussian_filter

class DocstringHelpParser(argparse.ArgumentParser):
    """ArgumentParser subclass that prints the module docstring verbatim for --help."""

    def format_help(self):
        doc = (__doc__ or '').strip()
        return f"{doc}\n" if doc else super().format_help()

def open_file(path):
    """Open a plain or gzip-compressed file, detected by magic bytes."""
    with open(path, 'rb') as f:
        magic = f.read(2)
    if magic == b'\x1f\x8b':
        return gzip.open(path, 'rt')
    return open(path, 'r')


def read_pandepth(path):
    """
    Read pandepth windowed output into a dict keyed by (chromosome, start_position).
    Accepts both gzipped and plain text files regardless of file extension.
    Skips comment lines starting with '#'.
    Expected columns: Chr, Start, End, Length, CoveredSite, TotalDepth, Coverage(%), MeanDepth
    """
    depth_dict = {}
    with open_file(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 8:
                continue
            chrom = parts[0]
            start = int(parts[1])
            mean_depth = float(parts[7])
            depth_dict[(chrom, start)] = mean_depth
    return depth_dict


def compute_gc_windows(fasta_path, depth_dict, window_size):
    """
    Stream through FASTA chromosome by chromosome, compute GC% per window,
    and match against depth_dict using (chromosome, start_position) keys.
    Windows with no matching depth entry are skipped.
    Returns (gc_array, depth_array) as float32 numpy arrays.
    """
    gc_list = []
    depth_list = []

    def process_chrom(chrom, seq):
        seq = seq.upper()
        seq_len = len(seq)
        pos = 0
        win_start = 1  # pandepth uses 1-based start positions
        while pos < seq_len:
            win_seq = seq[pos:pos + window_size]
            actual_len = len(win_seq)
            if actual_len == 0:
                break
            key = (chrom, win_start)
            if key in depth_dict:
                gc_count = win_seq.count('G') + win_seq.count('C')
                gc_pct = gc_count / actual_len * 100.0
                gc_list.append(gc_pct)
                depth_list.append(depth_dict[key])
            pos += window_size
            win_start += window_size

    current_chrom = None
    seq_parts = []

    with open_file(fasta_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_chrom is not None and seq_parts:
                    process_chrom(current_chrom, ''.join(seq_parts))
                current_chrom = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if current_chrom is not None and seq_parts:
            process_chrom(current_chrom, ''.join(seq_parts))

    return np.array(gc_list, dtype=np.float32), np.array(depth_list, dtype=np.float32)


def save_combined(output_path, gc_arr, depth_arr):
    """Save GC and depth arrays to a two-column TSV file (for --plot-only reuse)."""
    with open(output_path, 'w') as f:
        f.write('gc\tdepth\n')
        for gc, depth in zip(gc_arr, depth_arr):
            f.write(f'{gc:.4f}\t{depth:.4f}\n')


def load_combined(path):
    """Load a previously saved combined TSV produced by --output-data."""
    gc_list = []
    depth_list = []
    with open_file(path) as f:
        for i, line in enumerate(f):
            if i == 0 and line.startswith('gc'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            gc_list.append(float(parts[0]))
            depth_list.append(float(parts[1]))
    return np.array(gc_list, dtype=np.float32), np.array(depth_list, dtype=np.float32)


def create_visualization(gc_arr, depth_arr, output_file, log_depth=False):
    """Generate GC-Depth scatter plot with marginal histograms."""
    print("\n[Step 2] Creating visualization...")

    mask = (gc_arr >= 0) & (gc_arr <= 100) & (depth_arr > 0)
    gc = gc_arr[mask]
    depth = depth_arr[mask]

    if len(gc) == 0:
        print("Error: No valid data points after filtering.")
        sys.exit(1)

    avg_gc = float(np.mean(gc))
    avg_depth = float(np.mean(depth))
    median_gc = float(np.median(gc))
    median_depth = float(np.median(depth))

    print(f"  GC average: {avg_gc:.2f}%, median: {median_gc:.2f}%")
    print(f"  Depth average: {avg_depth:.2f}X, median: {median_depth:.2f}X")
    print(f"  Total windows: {len(gc)}")

    plt.style.use('seaborn-v0_8-whitegrid')
    mpl.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'DejaVu Sans', 'Helvetica'],
        'font.size': 14,
        'axes.labelsize': 16,
        'axes.titlesize': 14,
        'xtick.labelsize': 12,
        'ytick.labelsize': 12,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'savefig.facecolor': 'white',
        'axes.linewidth': 0.8,
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        'grid.alpha': 0.5,
        'grid.linestyle': '--',
        'grid.linewidth': 0.5,
        'grid.color': '#b0b0b0'
    })

    colors = {
        'histogram': '#b8b8b8',
        'histogram_edge': '#909090',
        'average_line': '#DC143C',
        'grid': '#c0c0c0'
    }

    cmap = LinearSegmentedColormap.from_list('gc_depth', ['#440154', '#21918c', '#fde725'])

    gc_min = float(np.min(gc))
    gc_max = float(np.max(gc))
    gc_range_min = gc_min - 5
    gc_range_max = gc_max + 5
    depth_max = float(np.quantile(depth, 0.99)) * 1.1

    fig = plt.figure(figsize=(11, 10), facecolor='white')
    gs = GridSpec(4, 4, figure=fig,
                  height_ratios=[1.2, 1.2, 1.2, 1.2],
                  width_ratios=[1, 1, 1, 1.3],
                  hspace=0.45, wspace=0.55)

    # Main scatter plot
    ax_main = fig.add_subplot(gs[1:4, 0:3])

    print("  Calculating point density...")
    x_bins = np.linspace(gc_range_min, gc_range_max, 100)
    y_bins = np.linspace(0, depth_max, 100)
    hist2d, x_edges, y_edges = np.histogram2d(gc, depth, bins=[x_bins, y_bins])
    hist2d_smooth = gaussian_filter(hist2d.T, sigma=1.5)

    x_idx = np.clip(np.digitize(gc, x_edges) - 1, 0, len(x_edges) - 2)
    y_idx = np.clip(np.digitize(depth, y_edges) - 1, 0, len(y_edges) - 2)
    density = hist2d_smooth[y_idx, x_idx]
    density_norm = (density - density.min()) / (density.max() - density.min() + 1e-10)

    ax_main.scatter(
        gc, depth,
        c=density_norm,
        s=15,
        alpha=0.7,
        cmap=cmap,
        edgecolors='none',
        rasterized=True
    )
    ax_main.axhline(y=median_depth, color=colors['average_line'],
                    linestyle='--', linewidth=1.5, alpha=0.8)
    ax_main.axvline(x=median_gc, color=colors['average_line'],
                    linestyle='--', linewidth=1.5, alpha=0.8)
    ax_main.set_xlabel(f'GC % (Average : {avg_gc:.2f} %)', fontsize=13, fontweight='bold')
    ax_main.set_ylabel(f'Depth (Average : {avg_depth:.2f} X)', fontsize=13, fontweight='bold')
    ax_main.set_xlim(gc_range_min, gc_range_max)
    if log_depth:
        ax_main.set_yscale('log')
        ax_main.set_ylim(float(np.min(depth)) * 0.8, depth_max)
    else:
        ax_main.set_ylim(0, depth_max)
    ax_main.grid(True, linestyle='--', alpha=0.5, color=colors['grid'], linewidth=0.5)
    ax_main.set_axisbelow(True)
    for spine in ax_main.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)
    ax_main.tick_params(axis='both', which='both', direction='out',
                        top=False, bottom=True, left=True, right=False,
                        length=5, width=0.8)

    # Top GC histogram
    ax_gc = fig.add_subplot(gs[0, 0:3])
    ax_gc.hist(gc, bins=60,
               color=colors['histogram'],
               edgecolor=colors['histogram_edge'],
               linewidth=0.3, alpha=0.9,
               range=(gc_range_min, gc_range_max))
    ax_gc.axvline(x=median_gc, color=colors['average_line'],
                  linestyle='--', linewidth=1.5, alpha=0.8)
    ax_gc.set_ylabel('Numbers', fontsize=13, fontweight='bold')
    ax_gc.set_xlabel('')
    ax_gc.set_xlim(gc_range_min, gc_range_max)
    ax_gc.xaxis.set_tick_params(labelbottom=True)
    ax_gc.grid(True, linestyle='--', alpha=0.5, color=colors['grid'], linewidth=0.5)
    ax_gc.set_axisbelow(True)
    for spine in ax_gc.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)
    ax_gc.tick_params(axis='both', which='both', direction='out',
                      top=False, bottom=True, left=True, right=False,
                      length=5, width=0.8)

    # Right depth histogram
    ax_depth = fig.add_subplot(gs[1:4, 3])
    depth_range = (float(np.min(depth)) * 0.8, depth_max) if log_depth else (0, depth_max)
    ax_depth.hist(depth, bins=60,
                  orientation='horizontal',
                  color=colors['histogram'],
                  edgecolor=colors['histogram_edge'],
                  linewidth=0.3, alpha=0.9,
                  range=depth_range)
    ax_depth.axhline(y=median_depth, color=colors['average_line'],
                     linestyle='--', linewidth=1.5, alpha=0.8)
    ax_depth.set_xlabel('Numbers', fontsize=13, fontweight='bold')
    ax_depth.tick_params(axis='x', labelsize=11)
    ax_depth.set_ylabel('')
    ax_depth.set_ylim(depth_range)
    ax_depth.yaxis.set_tick_params(labelleft=True)
    ax_depth.grid(True, linestyle='--', alpha=0.5, color=colors['grid'], linewidth=0.5)
    ax_depth.set_axisbelow(True)
    for spine in ax_depth.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(1.0)
    ax_depth.tick_params(axis='both', which='both', direction='out',
                         top=False, bottom=True, left=True, right=False,
                         length=5, width=0.8)

    fig.subplots_adjust(left=0.10, right=0.95, bottom=0.08, top=0.95,
                        hspace=0.5, wspace=0.5)

    print(f"  Saving plot to: {output_file}")
    plt.savefig(output_file, dpi=300, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.1)
    plt.close()

    return median_gc, median_depth, len(gc)


def main():
    parser = DocstringHelpParser(
        description='Compute and visualize GC content vs sequencing depth per genomic window.',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('fasta', nargs='?',
                        help='Genome FASTA file, chromosome or contig level (.fa/.fasta, plain or gzipped)')
    parser.add_argument('pandepth', nargs='?',
                        help='Pandepth windowed depth file (.win.stat.gz or plain text with the same format)')
    parser.add_argument('-w', '--window', type=int, default=1000,
                        help='Window size in bp, must match pandepth -w value (default: 1000)')
    parser.add_argument('-o', '--output', default='gc-depth.png',
                        help='Output plot file (.png or .pdf, default: gc-depth.png)')
    parser.add_argument('--log-depth', action='store_true',
                        help='Use logarithmic scale for the depth axis')
    parser.add_argument('--plot-only', metavar='TSV',
                        help='Skip processing, re-plot from an existing combined TSV (from --output-data)')
    parser.add_argument('--output-data', metavar='FILE',
                        help='Save merged GC and depth data to this TSV file (can be reused with --plot-only)')

    args = parser.parse_args()

    ext = os.path.splitext(args.output)[1].lower()
    if ext not in ('.png', '.pdf'):
        print("Error: Output file must use .png or .pdf extension.")
        sys.exit(1)

    print("gc-depth-plot")

    if args.plot_only:
        if not os.path.exists(args.plot_only):
            print(f"Error: TSV file not found: {args.plot_only}")
            sys.exit(1)
        print(f"\nPlot-only mode: {args.plot_only}")
        gc_arr, depth_arr = load_combined(args.plot_only)

    else:
        if not args.fasta or not args.pandepth:
            parser.error("fasta and pandepth arguments are required unless using --plot-only")
        if not os.path.exists(args.fasta):
            print(f"Error: FASTA file not found: {args.fasta}")
            sys.exit(1)
        if not os.path.exists(args.pandepth):
            print(f"Error: Pandepth file not found: {args.pandepth}")
            sys.exit(1)

        print(f"\nParameters:")
        print(f"  FASTA:       {args.fasta}")
        print(f"  Depth file:  {args.pandepth}")
        print(f"  Window size: {args.window} bp")
        print(f"  Output:      {args.output}")

        print("\n[Step 1] Processing data...")
        print("  Reading pandepth output...")
        depth_dict = read_pandepth(args.pandepth)
        print(f"  Loaded {len(depth_dict)} depth windows")

        print("  Computing GC content from FASTA...")
        gc_arr, depth_arr = compute_gc_windows(args.fasta, depth_dict, args.window)
        print(f"  Matched {len(gc_arr)} windows with both GC and depth data")

        if len(gc_arr) == 0:
            print("Error: No windows matched between FASTA and pandepth output.")
            print("  Check that sequence names (chromosomes/contigs) and window size match between inputs.")
            sys.exit(1)

        if args.output_data:
            save_combined(args.output_data, gc_arr, depth_arr)
            print(f"  Combined data saved to: {args.output_data}")

    median_gc, median_depth, n_windows = create_visualization(
        gc_arr, depth_arr, args.output, args.log_depth
    )

    print(f"\nResults:")
    print(f"- Visualization output : {args.output}")
    print(f"- Median GC            : {median_gc:.2f}%")
    print(f"- Median Depth         : {median_depth:.2f}X")
    print("\nDone!")


if __name__ == '__main__':
    main()
