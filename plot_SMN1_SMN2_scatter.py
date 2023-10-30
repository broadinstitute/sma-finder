import argparse
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import re
import seaborn as sns
import string

np.random.seed(1)

P_SMN1 = 0.2       # based on an individual having 1 copy of SMN1 and 4 copies of SMN2
P_ERROR = 0.001/3  # Q30
MIN_COVERAGE_THRESHOLD = 14
MAX_SMN1_READS_THRESHOLD = 2
PHRED_SCALE_CONFIDENCE_THRESHOLD = 10

UNKNOWN_SAMPLES_LABEL = "Rare Disease Cases"
NEGATIVE_CONTROL_LABEL = "Unaffected Relatives"
POSITIVE_CONTROL_LABEL = "Cases With Validated\nSMA Diagnosis"

PALETTE = {
    UNKNOWN_SAMPLES_LABEL: "tab:gray",
    NEGATIVE_CONTROL_LABEL: "tab:blue",
    POSITIVE_CONTROL_LABEL: "tab:red",
}

C840_SMN1_READS_COLUMN = "c840_reads_with_smn1_base_C"
C840_TOTAL_READS_COLUMN = "c840_total_reads"
SMA_STATUS_COLUMN = "sma_status"


def compute_decision_boundary_slope_and_intercept(p, p_error, confidence_threshold):
    """Compute the slope of the decision boundary line when calling SMN1 copy number (0 or greater than 0)
    using likelihood estimation. Returns m from N = m * r where
    r = number of reads that support the presence of SMN1 at the c.840 position in SMN exon 7
    N = total number of reads (ie. coverage) at the c.840 position in SMN exon 7

    Args:
        p (float): smallest possible probability of seeing a read from SMN1 given the individual has at least 1 copy of
            SMN1 (The smaller this probability, the less likely we call a sample with 1 or 2 reads as SMA positive)
        p_error (float): this is the base error rate - ie. the probability of a base being sequenced incorrectly.
            This is used for the reference confidence model. [Poplin 2018]
        confidence_threshold (float): the PHRED scale confidence threshold - defined as the difference in
            log-likelihoods between the two possible outcomes (0 copies of SMA, or more than 0 copies).
            This is similar to a difference in raw PL's defined by HaplotypeCaller [Poplin 2018]
    Return:
        float, float: Returns 2 numbers: the slope and y-intercept of the decision boundary line
    """
    numerator = math.log10((p * (1 - p_error))/(p_error * (1 - p)))
    denominator = math.log10((1 - p_error)/(1 - p))

    confidence_threshold /= 10
    return numerator / denominator, confidence_threshold / denominator


def plot_figure(df_all, title=None):

    # compute decision boundary x,y points.
    decision_boundary_x = np.linspace(0, 1000, 100000)

    # decision boundary 1 = without a confidence threshold (shown as a solid line)
    slope, intercept = compute_decision_boundary_slope_and_intercept(P_SMN1, P_ERROR, 0)
    decision_boundary_y = slope * decision_boundary_x + intercept
    decision_boundary_y[decision_boundary_y < MIN_COVERAGE_THRESHOLD - 0.5] = MIN_COVERAGE_THRESHOLD - 0.5

    # decision boundary 2 = with a confidence threshold (shown as a dashed line)
    #slope2, intercept2 = compute_decision_boundary_slope_and_intercept(P_SMN1, P_ERROR, PHRED_SCALE_CONFIDENCE_THRESHOLD)
    #decision_boundary2_y = slope2 * decision_boundary_x + intercept2
    #decision_boundary2_y[decision_boundary2_y < MIN_COVERAGE_THRESHOLD - 0.5] = MIN_COVERAGE_THRESHOLD - 0.5

    # initialize the figure
    fig, ax = plt.subplots(1, 1, figsize=(9.2, 6.5))
    plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05, hspace=0.3)
    fig.set_tight_layout(True)
    if title: fig.suptitle(title, fontsize=17)

    x_axis_label = "# of reads that contain a 'C' at the c.840 position"
    y_axis_label = "Total # of reads at the c.840 position"

    # add jitter to points below the minimum coverage threshold to make them visible
    jitter = (np.random.rand(len(df_all)) - 0.5)/4
    df_all.loc[
        df_all[C840_TOTAL_READS_COLUMN] < MIN_COVERAGE_THRESHOLD,
        C840_TOTAL_READS_COLUMN
    ] = df_all[C840_TOTAL_READS_COLUMN] + jitter

    # render data points
    sns.scatterplot(
        ax=ax,
        data=df_all,
        x=C840_SMN1_READS_COLUMN,
        y=C840_TOTAL_READS_COLUMN,
        hue=SMA_STATUS_COLUMN,
        palette=PALETTE,
        s=15)

    # text labels for points that should be labeled
    df_samples_to_label = df_all[df_all["show_label"]]

    if len(df_samples_to_label) > 0:
        for _, row in df_samples_to_label.iterrows():
            x_offset = 2
            y_offset = 50
            ax.text(
                row[C840_SMN1_READS_COLUMN] + x_offset, row[C840_TOTAL_READS_COLUMN] + y_offset,
                row["sample_id"], size=12)
            ax.arrow(
                row[C840_SMN1_READS_COLUMN] + x_offset - 0.5, row[C840_TOTAL_READS_COLUMN] + y_offset - 0.5,
                -x_offset + 0.5 + 0.3, -y_offset + 0.5 + 2,
                head_width=x_offset*0.1, head_length=y_offset*0.15)

    ax.set_xlim([-0.5, 1000])
    ax.set_ylim([-0.5, 1000])
    ax.set_xscale("symlog", linthresh=MIN_COVERAGE_THRESHOLD)
    ax.set_yscale("symlog", linthresh=MIN_COVERAGE_THRESHOLD)
    ax.tick_params(labelsize=14)
    ax.set_xlabel(x_axis_label, size=15, labelpad=15)
    ax.set_ylabel(y_axis_label, size=15, labelpad=15)

    # decision boundary line(s)
    decision_boundary_line_color = "red"
    ax.plot(
        decision_boundary_x,
        decision_boundary_y,
        '-',
        label='Decision Boundary',
        c=decision_boundary_line_color)

    # minimum coverage threshold line
    ax.plot(
        [0, MIN_COVERAGE_THRESHOLD-0.5],
        [MIN_COVERAGE_THRESHOLD-0.5, MIN_COVERAGE_THRESHOLD-0.5],
        linestyle=(0, (5, 1)),
        linewidth=0.8,
        color='black',
        label="Minimum Read\nCoverage Threshold")

    # legend
    ax.legend(
        loc='upper left',
        bbox_to_anchor=[1.05, 1],
        frameon=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-p", "--positive-control", action="append", help="Sample id of a positive control. "
                   "This argument can be specified more than once")
    p.add_argument("-P", "--positive-control-sample-ids-path", help="Path of text file containing a list of positive "
                   "control sample ids (one per line)")
    p.add_argument("-l", "--sample-id-to-label", action="append", help="Sample id that should be labeled in the plot. "
                   "This argument can be specified more than once")
    p.add_argument("-L", "--sample-ids-to-label-path", help="Path of text file containing a list of sample ids that "
                   "should be labeled in the plot")
    p.add_argument("--sample-id-column", default="sample_id", help="Name of the sample id column in the input table")
    p.add_argument("--affected-status-column", default="affected", help="Name of the input table column that contains "
                   "the affected status of each sample")
    p.add_argument("--not-affected-label", default="Not Affected", help="Value in the affected status column that "
                   "indicates that the sample is unaffected")
    p.add_argument("--title", help="Plot title")
    p.add_argument("-o", "--output-prefix", help="Filename prefix for the output plot image file.")
    p.add_argument("-f", "--output-format", default="png", choices=["png", "pdf", "svg"], help="Output image format")
    p.add_argument("sma_finder_combined_results_path", help="Path of tsv file containing the combined results of "
                   "SMA Finder for all samples")

    args = p.parse_args()

    for path in args.positive_control_sample_ids_path, args.sample_ids_to_label_path, args.sma_finder_combined_results_path:
        if path and not os.path.isfile(path):
            p.error(f"File not found: {args.positive_control_sample_ids_path}")

    df = pd.read_table(args.sma_finder_combined_results_path, dtype={
        "confidence_score": float,
        "c840_reads_with_smn1_base_C": float,
        "c840_total_reads": float,
    })

    for column_name in args.sample_id_column, args.affected_status_column, C840_SMN1_READS_COLUMN, C840_TOTAL_READS_COLUMN:
        if column_name not in df.columns:
            p.error(f"Column '{column_name}' not found in {args.sma_finder_combined_results_path}")

    all_sample_ids = set(df[args.sample_id_column])

    positive_controls = set()
    if args.positive_control_sample_ids_path:
        with open(args.positive_control_sample_ids_path, "rt") as f:
            positive_controls |= {s.strip() for s in f if s.strip()}

    sample_ids_to_label = set()
    if args.sample_ids_to_label_path:
        with open(args.sample_ids_to_label_path, "rt") as f:
            sample_ids_to_label |= {s.strip() for s in f if s.strip()}

    for source, sample_ids in [
        ("--positive-control", args.positive_control),
        ("--sample-id-to-label", args.sample_id_to_label),
        (args.positive_control_sample_ids_path, positive_controls),
        (args.sample_ids_to_label_path, sample_ids_to_label)]:
        if not sample_ids:
            continue
        unknown_sample_ids = set(sample_ids) - all_sample_ids
        if unknown_sample_ids:
            print(f"WARNING: Unknown sample id(s) specified in {source}: " + ", ".join(sorted(unknown_sample_ids)))

    if args.positive_control:
        positive_controls |= set(args.positive_control)

    if args.sample_id_to_label:
        sample_ids_to_label |= set(args.sample_id_to_label)

    df[SMA_STATUS_COLUMN] = UNKNOWN_SAMPLES_LABEL
    
    if args.not_affected_label not in set(df[args.affected_status_column]):
        print(f"WARNING: --not-affected-label arg value '{args.not_affected_label}' was not found in the "
              f"{args.affected_status_column} column of {args.sma_finder_combined_results_path}")
        
    df.loc[df[args.affected_status_column] == args.not_affected_label, SMA_STATUS_COLUMN] = NEGATIVE_CONTROL_LABEL
    
    if positive_controls:
        df.loc[df["sample_id"].isin(positive_controls), SMA_STATUS_COLUMN] = POSITIVE_CONTROL_LABEL
    df["show_label"] = df["sample_id"].isin(sample_ids_to_label)

    plot_figure(df, args.title)

    if not args.output_prefix:
        args.output_prefix = re.sub(".tsv(.gz)?$", "", os.path.basename(args.sma_finder_combined_results_path))

    ouput_path = f"{args.output_prefix}.SMN1_vs_SMN2_plot.{args.output_format}"
    plt.savefig(ouput_path)
    print(f"Saved plot to {ouput_path}")

    plt.close()


if __name__ == "__main__":
    main()

#%%


