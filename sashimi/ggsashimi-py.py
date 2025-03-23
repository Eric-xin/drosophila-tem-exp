#!/usr/bin/env python

import sys, re, copy, os, codecs, gzip, argparse, platform
from collections import OrderedDict
import pysam
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

__version__ = "1.1.5"

def get_version():
    prog = 'ggsashimi'
    version = '{} v{}'.format(prog, __version__)
    return version

def define_options():
    parser = argparse.ArgumentParser(description='Create sashimi plot for a given genomic region')
    parser.add_argument("-b", "--bam", type=str, required=True,
            help="""Individual bam file or file with a list of bam files.
            In the case of a list of files the format is tsv:
            1col: id for bam file,
            2col: path of bam file,
            3+col: additional columns""")
    parser.add_argument("-c", "--coordinates", type=str, required=True,
            help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
    parser.add_argument("-o", "--out-prefix", type=str, dest="out_prefix", default="sashimi",
            help="Prefix for plot file name [default=%(default)s]")
    parser.add_argument("-S", "--out-strand", type=str, dest="out_strand", default="both",
            help="Only for --strand other than 'NONE'. Choose which signal strand to plot: <both> <plus> <minus> [default=%(default)s]")
    parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
            help="Minimum number of reads supporting a junction to be drawn [default=1]")
    parser.add_argument("-j", "--junctions-bed", type=str, dest="junctions_bed", default="",
            help="Junction BED file name [default=no junction file]")
    parser.add_argument("-g", "--gtf",
            help="Gtf file with annotation (only exons is enough)")
    parser.add_argument("-s", "--strand", default="NONE", type=str,
            help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
    parser.add_argument("--shrink", action="store_true",
            help="Shrink the junctions by a factor for nicer display [default=%(default)s]")
    parser.add_argument("-O", "--overlay", type=int,
            help="Index of column with overlay levels (1-based)")
    parser.add_argument("-A", "--aggr", type=str, default="",
            help="""Aggregate function for overlay: <mean> <median> <mean_j> <median_j>.
                    Use mean_j | median_j to keep density overlay but aggregate junction counts [default=no aggregation]""")
    parser.add_argument("-C", "--color-factor", type=int, dest="color_factor",
            help="Index of column with color levels (1-based)")
    parser.add_argument("--alpha", type=float, default=0.5,
            help="Transparency level for density histogram [default=%(default)s]")
    parser.add_argument("-P", "--palette", type=str,
            help="Color palette file. TSV file with ≥1 columns, where the color is the first column. Both R color names and hexadecimal values are valid")
    parser.add_argument("-L", "--labels", type=int, dest="labels", default=1,
            help="Index of column with labels (1-based) [default=%(default)s]")
    parser.add_argument("--fix-y-scale", default=False, action="store_true", dest="fix_y_scale",
            help="Fix y-scale across individual signal plots [default=%(default)s]")
    parser.add_argument("--height", type=float, default=2,
            help="Height of the individual signal plot in inches [default=%(default)s]")
    parser.add_argument("--ann-height", type=float, default=1.5, dest="ann_height",
            help="Height of annotation plot in inches [default=%(default)s]")
    parser.add_argument("--width", type=float, default=10,
            help="Width of the plot in inches [default=%(default)s]")
    parser.add_argument("--base-size", type=float, default=14, dest="base_size",
            help="Base font size of the plot in points [default=%(default)s]")
    parser.add_argument("-F", "--out-format", type=str, default="pdf", dest="out_format",
            help="Output file format: <pdf> <svg> <png> <jpeg> <tiff> [default=%(default)s]")
    parser.add_argument("-R", "--out-resolution", type=int, default=300, dest="out_resolution",
            help="Output file resolution in DPI. Applies only to raster output formats [default=%(default)s]")
    parser.add_argument("--debug-info", action="store_true",
            help="Show system information useful for debugging purposes")
    parser.add_argument('--version', action='version', version=get_version())
    return parser

def parse_coordinates(c):
    c = c.replace(",", "")
    chrom = c.split(":")[0]
    start, end = c.split(":")[1].split("-")
    start, end = int(start) - 1, int(end)
    return chrom, start, end

def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions):
    if CIGAR_op == "M":
        for i in range(pos, pos + CIGAR_len):
            if i < start or i >= end:
                continue
            ind = i - start
            a[ind] += 1
    if CIGAR_op in ["I", "S"]:
        return pos
    if CIGAR_op == "D":
        pass
    if CIGAR_op == "N":
        don = pos
        acc = pos + CIGAR_len
        if don > start and acc < end:
            junctions[(don, acc)] = junctions.get((don, acc), 0) + 1
    pos = pos + CIGAR_len
    return pos

def flip_read(s, samflag):
    if s in ["NONE", "SENSE"]:
        return 0
    if s == "ANTISENSE":
        return 1
    if s == "MATE1_SENSE":
        if int(samflag) & 64:
            return 0
        if int(samflag) & 128:
            return 1
    if s == "MATE2_SENSE":
        if int(samflag) & 64:
            return 1
        if int(samflag) & 128:
            return 0
    return 0

def read_bam(f, c, s):
    chrom, start, end = parse_coordinates(c)
    a = {"+" : [0] * (end - start)}
    junctions = {"+": OrderedDict()}
    if s != "NONE":
        a["-"] = [0] * (end - start)
        junctions["-"] = OrderedDict()
    samfile = pysam.AlignmentFile(f)
    for read in samfile.fetch(chrom, start, end):
        if read.is_unmapped:
            continue
        samflag = read.flag
        read_start = read.reference_start + 1
        CIGAR = read.cigarstring
        if any(op in CIGAR for op in ["H", "P", "X", "="]):
            continue
        read_strand = ["+", "-"][flip_read(s, samflag) ^ bool(samflag & 16)]
        if s == "NONE":
            read_strand = "+"
        CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
        CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]
        pos = read_start
        for n, op in enumerate(CIGAR_ops):
            length = int(CIGAR_lens[n])
            pos = count_operator(op, length, pos, start, end, a[read_strand], junctions[read_strand])
    samfile.close()
    return a, junctions

def get_bam_path(index, path):
    if os.path.isabs(path):
        return path
    base_dir = os.path.dirname(index)
    return os.path.join(base_dir, path)

def read_bam_input(f, overlay, color, label):
    if f.endswith(".bam"):
        bn = os.path.basename(f).replace(".bam", "")
        yield bn, f, None, None, bn
        return
    with codecs.open(f, encoding='utf-8') as openf:
        for line in openf:
            parts = line.strip().split("\t")
            bam = get_bam_path(f, parts[1])
            overlay_level = parts[overlay-1] if overlay else None
            color_level = parts[color-1] if color else None
            label_text = parts[label-1] if label else None
            yield parts[0], bam, overlay_level, color_level, label_text

def prepare_for_plot(a, junctions, c, m):
    _, start, _ = parse_coordinates(c)
    x = [i + start for i in range(len(a))]
    y = a
    dons, accs, yd, ya, counts = [], [], [], [], []
    for (don, acc), n in junctions.items():
        if n < m:
            continue
        dons.append(don)
        accs.append(acc)
        counts.append(n)
        yd.append(a[don - start - 1] if 0 <= don - start - 1 < len(a) else 0)
        ya.append(a[acc - start + 1] if 0 <= acc - start + 1 < len(a) else 0)
    return x, y, dons, accs, yd, ya, counts

def intersect_introns(data):
    data = sorted(data)
    it = iter(data)
    a, b = next(it)
    for c, d in it:
        if b > c:
            b = min(b, d)
            a = max(a, c)
        else:
            yield a, b
            a, b = c, d
    yield a, b

def shrink_density(x, y, introns):
    new_x, new_y = [], []
    shift = 0
    start = 0
    for a_val, b_val in introns:
        try:
            end_index = x.index(a_val) + 1
        except ValueError:
            continue
        new_x += [i - shift for i in x[start:end_index]]
        new_y += y[start:end_index]
        start = x.index(b_val)
        l = (b_val - a_val)
        shift += l - (l ** 0.7)
    new_x += [i - shift for i in x[start:]]
    new_y += y[start:]
    return new_x, new_y

def shrink_junctions(dons, accs, introns):
    new_dons = [0]*len(dons)
    new_accs = [0]*len(accs)
    real_introns = dict()
    shift_acc = 0
    shift_don = 0
    s = set()
    junctions = list(zip(dons, accs))
    for a_val, b_val in introns:
        l = b_val - a_val
        shift_acc += l - int(l ** 0.7)
        real_introns[a_val - shift_don] = a_val
        real_introns[b_val - shift_acc] = b_val
        for i, (don, acc) in enumerate(junctions):
            if a_val >= don and b_val <= acc:
                if (don, acc) not in s:
                    new_dons[i] = don - shift_don
                    new_accs[i] = acc - shift_acc
                else:
                    new_accs[i] = acc - shift_acc
                s.add((don, acc))
        shift_don = shift_acc
    return real_introns, new_dons, new_accs

def read_palette(f):
    palette = ["#ff0000", "#00ff00", "#0000ff", "#000000"]
    if f:
        with open(f) as openf:
            palette = [line.split("\t")[0].strip() for line in openf if line.strip()]
    return palette

def read_gtf(f, c):
    exons = OrderedDict()
    transcripts = OrderedDict()
    chrom, start, end = parse_coordinates(c)
    end = end - 1
    opener = gzip.open(f, 'rt') if f.endswith(".gz") else open(f)
    with opener as openf:
        for line in openf:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            el_chr, _, el, el_start, el_end, _, strand, _, tags = parts
            if el_chr != chrom:
                continue
            if el not in ("transcript", "exon"):
                continue
            match = re.search('transcript_id "([^"]+)"', tags)
            if not match:
                print("ERROR: 'transcript_id' attribute is missing in the GTF file.")
                sys.exit(1)
            transcript_id = match.group(1)
            el_start, el_end = int(el_start) - 1, int(el_end)
            if el == "transcript":
                if el_end > start and el_start < end:
                    transcripts[transcript_id] = (max(start, el_start), min(end, el_end), strand)
                continue
            if el == "exon":
                if (start < el_start < end or start < el_end < end):
                    exons.setdefault(transcript_id, []).append((max(el_start, start), min(end, el_end), strand))
    return transcripts, exons

def make_introns(transcripts, exons, intersected_introns=None):
    new_transcripts = copy.deepcopy(transcripts)
    new_exons = copy.deepcopy(exons)
    introns = OrderedDict()
    if intersected_introns:
        for tx, (tx_start, tx_end, strand) in new_transcripts.items():
            total_shift = 0
            for a_val, b_val in intersected_introns:
                l = b_val - a_val
                shift = l - int(l ** 0.7)
                total_shift += shift
                for i, (exon_start, exon_end, strand) in enumerate(new_exons.get(tx, [])):
                    new_exon_start, new_exon_end = new_exons[tx][i][:2]
                    if a_val < exon_start:
                        if b_val > exon_end:
                            if i == len(new_exons[tx]) - 1:
                                total_shift = total_shift - shift + (exon_start - a_val) * (1 - int(l ** -0.3))
                            shift = (exon_start - a_val) * (1 - int(l ** -0.3))
                            new_exon_end = new_exons[tx][i][1] - shift
                        new_exon_start = new_exons[tx][i][0] - shift
                    if b_val <= exon_end:
                        new_exon_end = new_exons[tx][i][1] - shift
                    new_exons[tx][i] = (new_exon_start, new_exon_end, strand)
            tx_start = min(tx_start, sorted(new_exons.get(tx, [[sys.maxsize]]))[0][0])
            new_transcripts[tx] = (tx_start, tx_end - total_shift, strand)
    for tx, (tx_start, tx_end, strand) in new_transcripts.items():
        intron_start = tx_start
        ex_end = 0
        for ex_start, ex_end, strand in sorted(new_exons.get(tx, [])):
            intron_end = ex_start
            if tx_start < ex_start:
                introns.setdefault(tx, []).append((intron_start, intron_end, strand))
            intron_start = ex_end
        if tx_end > ex_end:
            introns.setdefault(tx, []).append((intron_start, tx_end, strand))
    return {'transcripts': new_transcripts, 'exons': new_exons, 'introns': introns}

def median(lst):
    s_lst = sorted(lst)
    n = len(s_lst)
    if n % 2 == 1:
        return s_lst[n//2]
    else:
        return sum(s_lst[n//2 - 1:n//2 + 1]) / 2.

def mean(lst):
    return sum(lst) / len(lst)

def aggregate_tracks(bam_data, overlay_dict, aggr, intersected_introns):
    aggregated = {}
    if not overlay_dict:
        for id, data in bam_data.items():
            x, y, dons, accs, yd, ya, counts = data
            if intersected_introns:
                x, y = shrink_density(x, y, intersected_introns)
                _, new_dons, new_accs = shrink_junctions(dons, accs, intersected_introns)
                dons, accs = new_dons, new_accs
            aggregated[id] = (x, y, dons, accs, yd, ya, counts)
    else:
        for overlay, ids in overlay_dict.items():
            combined_x = None
            combined_y = []
            combined_dons = []
            combined_accs = []
            combined_yd = []
            combined_ya = []
            combined_counts = []
            for id in ids:
                if id not in bam_data:
                    continue
                x, y, dons, accs, yd, ya, counts = bam_data[id]
                if intersected_introns:
                    x, y = shrink_density(x, y, intersected_introns)
                    _, new_dons, new_accs = shrink_junctions(dons, accs, intersected_introns)
                    dons, accs = new_dons, new_accs
                if combined_x is None:
                    combined_x = x
                    combined_y.append(y)
                else:
                    combined_y.append(y)
                combined_dons.extend(dons)
                combined_accs.extend(accs)
                combined_yd.extend(yd)
                combined_ya.extend(ya)
                combined_counts.extend(counts)
            if aggr and not aggr.endswith("_j"):
                if aggr == "mean":
                    agg_y = [mean(vals) for vals in zip(*combined_y)]
                elif aggr == "median":
                    agg_y = [median(list(vals)) for vals in zip(*combined_y)]
                else:
                    agg_y = combined_y[0]
                aggregated[overlay] = (combined_x, agg_y, combined_dons, combined_accs, combined_yd, combined_ya, combined_counts)
            else:
                aggregated[overlay] = (combined_x, combined_y[0], combined_dons, combined_accs, combined_yd, combined_ya, combined_counts)
    return aggregated

def draw_arc(ax, x1, x2, y1, y2, count, total_count, color, above=True):
    mid_x = (x1 + x2) / 2.0
    height = (x2 - x1) / 4.0 * (count / total_count) if total_count > 0 else (x2 - x1) / 4.0
    if not above:
        height = -height
    control_y = max(y1, y2) + height if above else min(y1, y2) + height
    t = np.linspace(0, 1, 100)
    bezier_y = (1-t)**2 * y1 + 2*(1-t)*t*control_y + t**2 * y2
    bezier_x = (1-t)**2 * x1 + 2*(1-t)*t*mid_x + t**2 * x2
    ax.plot(bezier_x, bezier_y, color=color, linewidth=1)

def plot_sashimi(aggregated_data, args, label_dict, color_dict, chrom, plot_range, annotation=None):
    num_tracks = len(aggregated_data)
    total_height = args.height * num_tracks
    if annotation:
        total_height += args.ann_height
    fig, axs = plt.subplots(num_tracks + (1 if annotation else 0), 1,
                            figsize=(args.width, total_height), sharex=True)
    if num_tracks == 1 and not annotation:
        axs = [axs]
    else:
        axs = np.atleast_1d(axs)
    global_max = 0
    for track, data in aggregated_data.items():
        x, y, dons, accs, yd, ya, counts = data
        global_max = max(global_max, max(y))
    track_index = 0
    for track, data in aggregated_data.items():
        ax = axs[track_index]
        x, y, dons, accs, yd, ya, counts = data
        color = color_dict.get(track, "grey")
        ax.bar(x, y, width=1, color=color, alpha=args.alpha)
        ax.set_ylabel(label_dict.get(track, track), fontsize=args.base_size)
        ax.set_xlim(plot_range)
        if args.fix_y_scale:
            ax.set_ylim(0, global_max * 1.2)
        total_junction = sum(counts) if counts else 0
        for i, (don, acc, ydon, yacc, count) in enumerate(zip(dons, accs, yd, ya, counts)):
            above = (i % 2 == 0)
            draw_arc(ax, don, acc, ydon, yacc, count, total_junction, color, above=above)
            mid_x = (don + acc) / 2.0
            mid_y = (ydon + yacc) / 2.0
            ax.text(mid_x, mid_y, str(count), fontsize=args.base_size * 0.7, color=color, ha='center')
        track_index += 1

    # Annotation panel for transcripts and exons
    if annotation:
        ax_ann = axs[-1]
        transcripts = annotation['transcripts']
        exons = annotation['exons']
        y_pos = 0.5
        for tx, (tx_start, tx_end, strand) in transcripts.items():
            ax_ann.hlines(y_pos, tx_start, tx_end, color="black", linewidth=2)
            # Place transcript label above the line
            ax_ann.text(tx_start, y_pos + 0.1, tx, fontsize=args.base_size * 0.8,
                        verticalalignment='bottom', horizontalalignment='left', color="blue")
            for exon in exons.get(tx, []):
                exon_start, exon_end, _ = exon
                ax_ann.add_patch(Rectangle((exon_start, y_pos - 0.1), exon_end - exon_start, 0.2, color="black"))
            y_pos += 0.5
        ax_ann.set_ylim(0, y_pos)
        ax_ann.set_xlim(plot_range)
        ax_ann.set_yticks([])
        ax_ann.set_title("Annotation", fontsize=args.base_size)
        # Remove spines for a cleaner look
        for spine in ax_ann.spines.values():
            spine.set_visible(False)
        ax_ann.set_xlabel(f"{chrom} genomic position", fontsize=args.base_size)
    else:
        axs[-1].set_xlabel(f"{chrom} genomic position", fontsize=args.base_size)
    plt.tight_layout()
    return fig

def print_debug_info():
    info = OrderedDict()
    system = platform.system()
    info["OS"] = "{}-{}".format(system, platform.machine())
    info["Python"] = platform.python_version()
    info["ggsashimi"] = __version__
    print(get_version())
    print('')
    maxlen = max(len(k) for k in info.keys())
    for k, v in info.items():
        print("{:{width}}: {}".format(k, v, width=maxlen))
    print('')

if __name__ == "__main__":
    strand_dict = {"plus": "+", "minus": "-"}
    parser = define_options()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    if args.debug_info:
        print_debug_info()
        sys.exit(0)
    if args.aggr and not args.overlay:
        print("ERROR: Cannot apply aggregate function if overlay is not selected.")
        sys.exit(1)
    palette = read_palette(args.palette)
    bam_dict = {"+": OrderedDict()}
    overlay_dict = OrderedDict()
    color_dict = OrderedDict()
    id_list = []
    label_dict = OrderedDict()
    if args.strand != "NONE":
        bam_dict["-"] = OrderedDict()
    junctions_list = []
    for id, bam, overlay_level, color_level, label_text in read_bam_input(args.bam, args.overlay, args.color_factor, args.labels):
        if not os.path.isfile(bam):
            continue
        a, junctions = read_bam(bam, args.coordinates, args.strand)
        if list(a.keys()) == ["+"] and all(v == 0 for v in list(a.values())[0]):
            print("WARN: Sample {} has no reads in the specified area.".format(id))
            continue
        id_list.append(id)
        label_dict[id] = label_text
        for strand in a:
            if args.strand == "NONE" or args.out_strand == 'both' or strand == strand_dict.get(args.out_strand, strand):
                if args.junctions_bed:
                    for k, v in junctions[strand].items():
                        if v >= args.min_coverage:
                            junctions_list.append('\t'.join([args.coordinates.split(':')[0], str(k[0]), str(k[1]), id, str(v), strand]))
            bam_dict[strand][id] = prepare_for_plot(a[strand], junctions[strand], args.coordinates, args.min_coverage)
        if color_level is None:
            color_dict.setdefault(id, id)
        if overlay_level is not None:
            overlay_dict.setdefault(overlay_level, []).append(id)
            label_dict[overlay_level] = overlay_level
            color_dict.setdefault(overlay_level, overlay_level)
        if overlay_level is None:
            color_dict.setdefault(id, color_level)
    if not bam_dict["+"]:
        print("ERROR: No available bam files.")
        sys.exit(1)
    if args.junctions_bed:
        if not args.junctions_bed.endswith('.bed'):
            args.junctions_bed = args.junctions_bed + '.bed'
        with open(args.junctions_bed, 'w') as jbed:
            jbed.write('\n'.join(sorted(junctions_list)))
    annotation = None
    if args.gtf:
        transcripts, exons = read_gtf(args.gtf, args.coordinates)
        annotation = make_introns(transcripts, exons, None)
    chrom, start, end = parse_coordinates(args.coordinates)
    plot_range = (start, end)
    # Map each key in color_dict to a valid color from the palette
    if palette:
        new_color_dict = {}
        keys = sorted(color_dict.keys())
        for i, key in enumerate(keys):
            new_color_dict[key] = palette[i % len(palette)]
        color_dict = new_color_dict
    for strand in bam_dict:
        if args.out_prefix.endswith(('.pdf', '.png', '.svg', '.tiff', '.tif', '.jpeg', '.jpg')):
            base, ext = os.path.splitext(args.out_prefix)
            if (args.out_format == ext[1:] or
                (args.out_format == 'tiff' and ext in ('.tiff','.tif')) or
                (args.out_format == 'jpeg' and ext in ('.jpeg','.jpg'))):
                args.out_prefix = base
                out_suffix = ext[1:]
            else:
                out_suffix = args.out_format
        else:
            out_suffix = args.out_format
        out_prefix = args.out_prefix + "_" + strand if args.strand != "NONE" else args.out_prefix
        intersected_introns = None
        if args.shrink:
            intron_candidates = (v for vs in bam_dict[strand].values() for v in zip(vs[2], vs[3]))
            intersected_introns = list(intersect_introns(intron_candidates))
        aggregated_data = aggregate_tracks(bam_dict[strand], overlay_dict, args.aggr, intersected_introns)
        fig = plot_sashimi(aggregated_data, args, label_dict, color_dict, chrom, plot_range, annotation)
        out_file = "{}.{}".format(out_prefix, out_suffix)
        plt.savefig(out_file, dpi=args.out_resolution)
        plt.close(fig)
    sys.exit(0)