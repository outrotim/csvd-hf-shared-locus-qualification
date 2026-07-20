#!/usr/bin/env python3
"""Redraw source-independent publication displays from aggregate results only.

Figures 1, 2 and 4 are redrawn from released aggregate values. For Figure 3,
the script reproduces the qualification-state and prior-sensitivity components;
exact regional association panels require provider summary data that are not
redistributed.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/publication_results.json"
OUT = Path(os.environ.get("OUTPUT_DIR", ROOT / "output"))
BLUE = "#1769AA"
ORANGE = "#D97706"
RED = "#B91C1C"
GREY = "#777777"
TEAL = "#009E73"


def as_bool(values: pd.Series) -> np.ndarray:
    return values.astype(str).str.lower().isin({"true", "1", "yes"}).to_numpy()


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure1(payload: dict) -> None:
    flow = payload["flow"]
    steps = [
        ("External cSVD\ncandidate mapping", "55 exposures → 56 assays\n51 proteins → 48 genes", BLUE),
        ("HF MR technically\nidentifiable", f"{flow['hermes_identifiable_proteins']}/51 proteins", BLUE),
        ("HERMES BH-FDR\n< 0.05", "3 proteins\n2 genes", ORANGE),
        ("Rule-opened regional\nevaluations", "3 evaluations\n2 genomic regions", ORANGE),
        ("Primary sharing\ngate", "0 supported\nevaluations", RED),
        ("Downstream molecular\nanalyses", "Not opened", GREY),
    ]
    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.set(xlim=(0, len(steps)), ylim=(0, 1))
    ax.axis("off")
    for index, (title, value, color) in enumerate(steps):
        x, width = index + 0.08, 0.84
        ax.add_patch(FancyBboxPatch((x, 0.26), width, 0.48,
            boxstyle="round,pad=0.02,rounding_size=0.04", linewidth=1.5,
            edgecolor=color, facecolor="white"))
        ax.text(x + width / 2, 0.61, title, ha="center", va="center",
                fontsize=8.2, weight="bold", linespacing=1.05)
        ax.text(x + width / 2, 0.40, value, ha="center", va="center",
                fontsize=6.9 if index == 0 else 9, color=color)
        if index < len(steps) - 1:
            ax.annotate("", xy=(index + 1.04, 0.50), xytext=(index + 0.94, 0.50),
                        arrowprops={"arrowstyle": "->", "color": "#4D4D4D", "lw": 1.3})
    ax.set_title("Prespecified evidence-qualification workflow and observed stopping point",
                 fontsize=13, pad=8)
    ax.text(0.08, 0.08, "Candidate mapping and exclusions were fixed before candidate-specific "
            "HERMES outcome access. No regional evaluation met the sharing criterion, so "
            "downstream molecular analyses were not opened as rescue.", fontsize=8.6)
    save(fig, "Figure_1_evidence_qualification_funnel")


def figure2(payload: dict) -> None:
    plot = pd.DataFrame(payload["proteins"]).sort_values(
        "simes_pvalue_hermes", ascending=False).reset_index(drop=True)
    y = np.arange(len(plot))
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 12.8), sharey=True)
    panels = [
        ("hermes", BLUE, "a  Primary HERMES no-UKB HF", 1.0, "white"),
        ("bbj", ORANGE, "b  Post hoc robustness: BBJ chronic HF", 0.72, "#FAFAFA"),
    ]
    labels = [f"{p} ({g})" for p, g in zip(plot["protein"], plot["gene"])]
    for ax, (source, color, title, alpha, facecolor) in zip(axes, panels):
        values = -np.log10(np.clip(plot[f"simes_pvalue_{source}"].astype(float), 1e-300, 1))
        passed = as_bool(plot[f"hf_fdr_pass_{source}"])
        identifiable = (plot[f"technical_state_{source}"] == "IDENTIFIABLE").to_numpy()
        ax.set_facecolor(facecolor)
        ax.scatter(values[identifiable], y[identifiable], s=28, color=color, alpha=alpha, zorder=3)
        ax.scatter(values[passed], y[passed], s=88, facecolors="none", edgecolors=color,
                   linewidths=1.5, alpha=alpha, zorder=4)
        ax.scatter(np.zeros((~identifiable).sum()), y[~identifiable], s=34, marker="x",
                   color=GREY, linewidths=1.2, zorder=4)
        ax.set_xlabel("−log10 protein-level Simes P value", labelpad=8)
        ax.set_title(title, loc="left", fontsize=11, weight="bold")
        ax.grid(axis="x", color="#EEEEEE", linewidth=0.7)
        ax.set_ylim(-1, len(plot))
    axes[0].set_yticks(y, labels, fontsize=7.4)
    axes[1].tick_params(axis="y", labelleft=False)
    fig.suptitle("Complete 51-protein MR family: primary HERMES and post hoc BBJ robustness",
                 fontsize=13, y=0.995)
    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    fig.text(0.5, 0.016, "Large open symbols indicate BH-FDR <0.05 within the corresponding "
             "51-protein family; grey crosses indicate technical ineligibility.\nPanel b is "
             "post hoc and could not rescue the primary result.", ha="center", fontsize=8)
    save(fig, "Figure_2_complete_51_protein_mr")


def figure3(payload: dict) -> None:
    loci = pd.DataFrame(payload["regional_evaluations"])
    prior = pd.DataFrame(payload["abf_prior_sensitivity"])
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 4.4),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    ax.axis("off")
    cells = []
    for row in loci.itertuples(index=False):
        no_pair = row.evidence_state == "INCONCLUSIVE_NO_CREDIBLE_SET_PAIR"
        resolution = ("No credible-set pair\nPP3/PP4 not estimable" if no_pair
                      else f"PP4={row.primary_pp4:.3f}\nPP3={row.primary_pp3:.3f}")
        state = "Resolution-limited\nnot qualified" if no_pair else "Prior-sensitive\nnot qualified"
        cells.append([f"{row.protein}\n({row.gene})", "SuSiE" if row.locus_model == "SUSIE"
                      else row.locus_model, f"{int(row.shared_variants):,}", resolution, state])
    table = ax.table(cellText=cells, colLabels=["Candidate", "Model", "Shared\nvariants",
                     "Regional resolution", "Evidence state"], cellLoc="center", loc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.8)
    for (r, c), cell in table.get_celld().items():
        if r == 0:
            cell.set_facecolor("#F0F0F0")
            cell.set_text_props(weight="bold")
        elif c == 4:
            cell.set_facecolor("#FFF2CC" if r == 2 else "#EAF2F8")
    ax.set_title("a  Three assay-specific evaluations across two genomic regions",
                 loc="left", fontsize=11, weight="bold")
    ax = axes[1]
    ax.semilogx(prior["p12"], prior["pp4"], marker="o", color=BLUE, label="PP4 (shared signal)")
    ax.semilogx(prior["p12"], prior["pp3"], marker="s", color=ORANGE, label="PP3 (distinct signals)")
    ax.axhline(0.8, color="#555555", linestyle="--", lw=0.9, label="PP4 gate 0.8")
    ax.axvline(1e-5, color="#555555", linestyle=":", lw=0.9, label="Primary p12")
    ax.set(ylim=(-0.03, 1.03), xlabel="Prior probability p12", ylabel="Posterior probability")
    ax.set_title("b  Apo E2 ABF prior sensitivity", loc="left", fontsize=11, weight="bold")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.text(0.03, 0.12, "Apo E and BGAT: SuSiE formed no credible-set pair;\nPP3/PP4 are "
            "therefore not estimable, not zero.", transform=ax.transAxes, fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.92, "pad": 2.0})
    fig.suptitle("Regional sharing evidence remained resolution-limited or prior-sensitive",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    save(fig, "Figure_3_locus_evidence_and_prior_sensitivity")


def figure4(payload: dict) -> None:
    data = pd.DataFrame(payload["assay_precision"])
    fig = plt.figure(figsize=(10.8, 8.0))
    grid = fig.add_gridspec(
        2, 2, left=0.085, right=0.97, top=0.85, bottom=0.105,
        height_ratios=[1.0, 1.1], hspace=0.36, wspace=0.16,
    )
    fig.suptitle(
        "Conditional precision separates near-null estimates from information limits",
        x=0.085, y=0.965, ha="left", fontsize=15, weight="bold",
    )
    fig.text(
        0.085, 0.925,
        "Assay-level 80% MDEs are descriptive source-scale quantities—not equivalence margins or common concentration effects",
        ha="left", fontsize=9.2, color=GREY,
    )
    sources = [
        ("HERMES_NO_UKB", BLUE, "a  Primary HERMES no-UKB HF"),
        ("BBJ_ONLY_EAS", ORANGE, "b  Post hoc BBJ chronic HF"),
    ]
    for ax, (source, color, title) in zip(
        [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])], sources
    ):
        sub = data[data["outcome_source"] == source].dropna(
            subset=["mde_abs_log_or_80pct_nominal", "beta"]
        ).copy()
        sub["mde_abs_log_or_80pct_nominal"] = pd.to_numeric(
            sub["mde_abs_log_or_80pct_nominal"]
        )
        sub["beta"] = pd.to_numeric(sub["beta"])
        sub = sub.sort_values("mde_abs_log_or_80pct_nominal").reset_index(drop=True)
        x = np.arange(1, len(sub) + 1)
        mde = sub["mde_abs_log_or_80pct_nominal"].to_numpy(float)
        beta = np.abs(sub["beta"].to_numpy(float))
        ax.vlines(x, np.minimum(beta, mde), np.maximum(beta, mde), color="#B8B8B8", alpha=0.35, lw=0.8)
        ax.scatter(x, mde, marker="s", s=25, color=color, label="Conditional 80% MDE")
        ax.scatter(x, beta, marker="o", s=20, facecolor="white", edgecolor="#17324D", linewidth=0.9, label="Absolute MR estimate")
        ax.set_yscale("log")
        ax.set_title(f"{title}  (n={len(sub)} estimable assays)", loc="left", fontsize=10.5, weight="bold")
        ax.set_xlabel("Assays ordered by conditional MDE")
        ax.set_ylabel("Absolute log-odds magnitude\n(source-specific assay scale)")
        ax.grid(axis="y", color="#E5E7EB", linewidth=0.7)
        ax.legend(frameon=False, fontsize=7.4, loc="upper left")

    ax = fig.add_subplot(grid[1, :])
    keys = ["exposure_id", "soma_id"]
    h = data[data["outcome_source"] == "HERMES_NO_UKB"][
        keys + ["fluid", "protein", "gene", "mde_abs_log_or_80pct_nominal"]
    ].rename(columns={"mde_abs_log_or_80pct_nominal": "hermes_mde"})
    b = data[data["outcome_source"] == "BBJ_ONLY_EAS"][
        keys + ["mde_abs_log_or_80pct_nominal"]
    ].rename(columns={"mde_abs_log_or_80pct_nominal": "bbj_mde"})
    paired = h.merge(b, on=keys, validate="one_to_one").dropna(
        subset=["hermes_mde", "bbj_mde"]
    )
    paired["hermes_mde"] = pd.to_numeric(paired["hermes_mde"])
    paired["bbj_mde"] = pd.to_numeric(paired["bbj_mde"])
    fluid_colors = {"CSF": TEAL, "PLASMA": "#CC79A7"}
    for fluid, sub in paired.groupby("fluid", sort=False):
        ax.scatter(
            sub["hermes_mde"], sub["bbj_mde"], s=34, alpha=0.72,
            color=fluid_colors.get(str(fluid), GREY), edgecolor="white", linewidth=0.5,
            label=f"{fluid} assay",
        )
    low = min(paired["hermes_mde"].min(), paired["bbj_mde"].min()) * 0.8
    high = max(paired["hermes_mde"].max(), paired["bbj_mde"].max()) * 1.25
    ax.plot([low, high], [low, high], color=GREY, linestyle="--", lw=1.0, label="Equal precision")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_xlabel("HERMES conditional 80% MDE")
    ax.set_ylabel("BBJ conditional 80% MDE")
    ax.set_title("c  Matched-assay precision shift", loc="left", fontsize=10.5, weight="bold")
    ax.grid(color="#E5E7EB", linewidth=0.7)
    ratio = paired["bbj_mde"] / paired["hermes_mde"]
    ax.text(
        0.025, 0.95,
        f"{int((ratio > 1).sum())}/{len(paired)} matched assays above identity\nMedian BBJ/HERMES MDE ratio = {ratio.median():.2f}",
        transform=ax.transAxes, ha="left", va="top", fontsize=8.4, weight="bold",
    )
    ax.legend(frameon=False, loc="lower right", ncol=3, fontsize=7.6)
    fig.text(
        0.5, 0.025,
        "Comparisons are paired within assay across HF outcomes. MDEs describe observed precision and do not establish equivalence or absence.",
        ha="center", fontsize=8.4, color=GREY,
    )
    save(fig, "Figure_4_precision_interpretability")


def main() -> None:
    payload = json.loads(DATA.read_text(encoding="utf-8"))
    figure1(payload)
    figure2(payload)
    figure3(payload)
    figure4(payload)
    print(f"Wrote Figures 1–4 aggregate displays to {OUT}")


if __name__ == "__main__":
    main()
