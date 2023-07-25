
#%%
import argparse
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
POSITIVE_CONTROL_LABEL = "Cases With Confirmed\n" \
                         "SMA Diagnosis"

PALETTE = {
    UNKNOWN_SAMPLES_LABEL: "tab:gray",
    NEGATIVE_CONTROL_LABEL: "tab:blue",
    POSITIVE_CONTROL_LABEL: "tab:red",
}


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


def plot_figure(df_all):

    exon_coverage_columns = [c for c in df_all.columns if c.startswith("reads_at_exon")]
    df_all.loc[:, "Average reads at exons 1-6"] = df_all[exon_coverage_columns].astype(int).sum(axis=1)/len(
        exon_coverage_columns)

    fig, ax = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    plt.subplots_adjust(left=0.05, right=0.85, top=0.95, bottom=0.05, hspace=0.3)
    fig.set_tight_layout(True)
    #fig.suptitle('sma_finder results for WGS, WES, and RNA samples')

    df_wgs = df_all[df_all.sample_type == "WGS"]
    df_wes = df_all[df_all.sample_type == "WES"]
    n_missing_sample_type = sum(df_all.sample_type.isna())
    assert len(df_wgs) + len(df_wes) + n_missing_sample_type == len(df_all)

    if n_missing_sample_type > 0:
        print(f"WARNING: sample type is missing for {n_missing_sample_type} samples")

    # compute decision boundary x,y points.
    # Decision boundary 1 = without a confidence threshold (shown as a solid line)
    # Decision boundary 2 = with a confidence threshold (shown as a dashed line)
    slope, intercept = compute_decision_boundary_slope_and_intercept(P_SMN1, P_ERROR, 0)
    slope2, intercept2 = compute_decision_boundary_slope_and_intercept(P_SMN1, P_ERROR, PHRED_SCALE_CONFIDENCE_THRESHOLD)

    decision_boundary_x_for_all_samples = np.linspace(0, 1000, 100000)
    decision_boundary_y = slope * decision_boundary_x_for_all_samples + intercept
    decision_boundary_y[decision_boundary_y < MIN_COVERAGE_THRESHOLD - 0.5] = MIN_COVERAGE_THRESHOLD - 0.5

    decision_boundary2_y = slope2 * decision_boundary_x_for_all_samples + intercept2
    decision_boundary2_y[decision_boundary2_y < MIN_COVERAGE_THRESHOLD - 0.5] = MIN_COVERAGE_THRESHOLD - 0.5

    # generate the 4 panels of Figure 3
    x_axis_label = "# of reads that contain a 'C' at the c.840 position"
    y_axis_label, total_reads_column = "Total # of reads at the c.840 position", "c840_total_reads"
    decision_boundary_x = decision_boundary_x_for_all_samples

    for ax_j, (df_for_sample_type, sample_type) in enumerate([(df_wgs, "WGS"), (df_wes, "WES")]):

        # add jitter to points below the minimum coverage threshold
        jitter = (np.random.rand(len(df_for_sample_type)) - 0.5)/4
        df_for_sample_type.loc[
            df_for_sample_type["c840_total_reads"] < MIN_COVERAGE_THRESHOLD,
            total_reads_column
        ] = df_for_sample_type[total_reads_column] + jitter

        current_ax = ax[ax_j]

        # render data points
        sns.scatterplot(
            ax=current_ax,
            data=df_for_sample_type,
            x="c840_reads_with_smn1_base_C",
            y=total_reads_column,
            hue="sma_status",
            palette=PALETTE,
            s=15)

        # add text label for new points
        df_samples_to_label = df_for_sample_type[df_for_sample_type["show_label"]]

        if len(df_samples_to_label) > 0:
            for _, row in df_samples_to_label.iterrows():
                x_offset = 2
                y_offset = 50
                current_ax.text(
                    row["c840_reads_with_smn1_base_C"] + x_offset, row[total_reads_column] + y_offset,
                    row["sample_id"], size=12)
                current_ax.arrow(
                    row["c840_reads_with_smn1_base_C"] + x_offset - 0.5, row[total_reads_column] + y_offset - 0.5,
                    -x_offset + 0.5 + 0.3, -y_offset + 0.5 + 2,
                    head_width=x_offset*0.1, head_length=y_offset*0.15)

        current_ax.set_xlim([-0.5, 1000])
        current_ax.set_ylim([-0.5, 1000])
        current_ax.set_xscale("symlog", linthresh=MIN_COVERAGE_THRESHOLD)
        current_ax.set_yscale("symlog", linthresh=MIN_COVERAGE_THRESHOLD)
        current_ax.tick_params(labelsize=14)
        current_ax.set_xlabel(x_axis_label, size=15, labelpad=15)

        current_ax.set_ylabel(y_axis_label if ax_j == 0 else "", size=15, labelpad=15)


        # A, B labels above each panel
        panel_title = f"All {sample_type} samples"
        current_ax.text(-0.1, 1.07, string.ascii_uppercase[ax_j], transform=current_ax.transAxes, size=20, weight='bold')
        current_ax.text(-0.01, 1.07, panel_title, transform=current_ax.transAxes, size=18)

        # decision boundary line (s)
        decision_boundary_line_color = "red"
        current_ax.plot(
            decision_boundary_x,
            decision_boundary_y,
            '-',
            label='Decision Boundary',
            c=decision_boundary_line_color)

        # minimum coverage threshold line
        current_ax.plot(
            [0, MIN_COVERAGE_THRESHOLD-0.5],
            [MIN_COVERAGE_THRESHOLD-0.5, MIN_COVERAGE_THRESHOLD-0.5],
            linestyle=(0, (5, 1)),
            linewidth=0.8,
            color='black',
            label="Minimum Read\nCoverage Threshold")

        # legend
        if ax_j > 0:
            current_ax.legend(
                loc='upper left',
                bbox_to_anchor=[1.05, 1],
                frameon=False)
        else:
            current_ax.get_legend().remove()

    #fig.subplots_adjust(right=0.8, left=0.05)  # , top=0.9, bottom=0.1, wspace=0.2


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--positive-controls-path", default="positive_controls.txt",
                   help="Path of text file containing a list of positive control sample ids (one per line)")
    p.add_argument("--no-label-path", default="no_label_samples.txt",
                   help="Path of text file containing a list of sample ids that shouldn't be labeled on the plot")    
    p.add_argument("-o", "--output-prefix", default="figure3",
                   help="Output prefix")
    p.add_argument("sma_finder_combined_results_path",
                   help="Path of tsv file containing the combined results of SMA Finder for all WES and WGS samples")

    args, _ = p.parse_known_args()

    positive_controls = set()
    if args.positive_controls_path:
        with open(args.positive_controls_path, "rt") as f:
            positive_controls = {s.strip() for s in f if s.strip()}
        print(f"Loaded {len(positive_controls)} positive controls: ", positive_controls)

    no_label_sample_ids = set()
    if args.no_label_path:
        with open(args.no_label_path, "rt") as f:
            no_label_sample_ids = {s.strip() for s in f if s.strip()}
        print(f"Loaded {len(no_label_sample_ids)} no-label sample ids")

    df = pd.read_table(args.sma_finder_combined_results_path)
    df.loc[:, "sma_status"] = UNKNOWN_SAMPLES_LABEL
    df.loc[df.affected == "Not Affected", "sma_status"] = NEGATIVE_CONTROL_LABEL
    df.loc[df["sample_id"].isin(positive_controls), "sma_status"] = POSITIVE_CONTROL_LABEL
    df["show_label"] = ~df["sample_id"].isin(no_label_sample_ids) & (df["sma_status"] == UNKNOWN_SAMPLES_LABEL) & (
            df["c840_reads_with_smn1_base_C"] <= MAX_SMN1_READS_THRESHOLD) & (
            df["c840_total_reads"] >= MIN_COVERAGE_THRESHOLD)

    plot_figure(df)

    plt.savefig(f"SMN1_vs_SMN2_plot.png")

    plt.show()
    plt.close()


if __name__ == "__main__":
    main()

#%%


