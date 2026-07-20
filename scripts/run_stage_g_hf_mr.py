#!/usr/bin/env python3
"""Run the frozen Stage-G protein-to-heart-failure MR screen.

This program is deliberately written and checksummed before any
candidate-specific HERMES value is opened.  It requires both primary exposure
registries (updated public WUSTL CSF and official deCODE plasma); the WUSTL
plasma registry is analysed only as an independently labelled sensitivity.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import chi2, norm


ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "results/stage_g/csvd_candidate_universe.tsv"
WUSTL_INSTRUMENTS = ROOT / "results/stage_g/wustl_candidate_cis_instruments.tsv"
DECODE_INSTRUMENTS = ROOT / "results/stage_g/decode_candidate_cis_instruments.tsv"
HERMES = ROOT / "data/raw/genetics/HERMES_noUKB_harmonized.parquet"
OUTDIR = ROOT / "results/stage_g/hf_mr"

PRIMARY_ROLES = {
    "CSF": "CSF_PRIMARY_UPDATED_REESTIMATION",
    "PLASMA": "PLASMA_DECODE_PRIMARY",
}
SENSITIVITY_ROLE = "PLASMA_WUSTL_INDEPENDENT_SENSITIVITY"
PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
F_MINIMUM = 10.0
PALINDROMIC_MID_LOW = 0.42
PALINDROMIC_MID_HIGH = 0.58
EAF_TOLERANCE = 0.10


def sha256(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def clean_allele(value: object) -> str:
    return str(value).strip().upper()


def finite_float(value: object) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def other_exposure_allele(row: pd.Series) -> str | None:
    explicit = clean_allele(row.get("other_allele", ""))
    if explicit and explicit not in {"NAN", "NA", "."}:
        return explicit
    effect = clean_allele(row["effect_allele"])
    ref = clean_allele(row.get("ref", ""))
    alt = clean_allele(row.get("alt", ""))
    if effect == ref and alt:
        return alt
    if effect == alt and ref:
        return ref
    return None


def complement(allele: str) -> str | None:
    if len(allele) != 1:
        return None
    return COMPLEMENT.get(allele)


def harmonize_one(exposure: pd.Series, outcome: pd.Series) -> dict:
    """Harmonize one instrument; unresolved rows are retained with a reason."""

    result = exposure.to_dict()
    result.update(
        {
            "outcome_a1": clean_allele(outcome["A1"]),
            "outcome_a2": clean_allele(outcome["A2"]),
            "beta_outcome_raw": finite_float(outcome["BETA"]),
            "se_outcome": finite_float(outcome["SE"]),
            "p_outcome": finite_float(outcome["P"]),
            "eaf_outcome": finite_float(outcome.get("EAF")),
            "harmonization_status": "UNRESOLVED",
            "mr_keep": False,
            "beta_outcome_aligned": np.nan,
        }
    )
    ea = clean_allele(exposure["effect_allele"])
    oa = other_exposure_allele(exposure)
    a1 = result["outcome_a1"]
    a2 = result["outcome_a2"]
    bx = finite_float(exposure["beta_exposure"])
    sx = finite_float(exposure["se_exposure"])
    by = result["beta_outcome_raw"]
    sy = result["se_outcome"]
    if oa is None or bx is None or sx is None or by is None or sy is None or sx <= 0 or sy <= 0:
        result["harmonization_status"] = "INVALID_EFFECT_OR_ALLELE"
        return result

    options: list[tuple[int, str]] = []
    if (a1, a2) == (ea, oa):
        options.append((1, "EXACT"))
    if (a1, a2) == (oa, ea):
        options.append((-1, "SWAPPED"))
    ca1, ca2 = complement(a1), complement(a2)
    if ca1 is not None and ca2 is not None:
        if (ca1, ca2) == (ea, oa):
            options.append((1, "COMPLEMENT"))
        if (ca1, ca2) == (oa, ea):
            options.append((-1, "COMPLEMENT_SWAPPED"))
    options = list(dict.fromkeys(options))
    if not options:
        result["harmonization_status"] = "ALLELE_MISMATCH"
        return result

    is_palindromic = frozenset((ea, oa)) in PALINDROMIC
    if is_palindromic:
        ex = finite_float(exposure.get("effect_allele_frequency"))
        oy = result["eaf_outcome"]
        if ex is None or oy is None:
            result["harmonization_status"] = "PALINDROMIC_MISSING_FREQUENCY"
            return result
        if (
            PALINDROMIC_MID_LOW <= ex <= PALINDROMIC_MID_HIGH
            or PALINDROMIC_MID_LOW <= oy <= PALINDROMIC_MID_HIGH
        ):
            result["harmonization_status"] = "PALINDROMIC_MID_FREQUENCY"
            return result
        compatible = []
        for sign, label in options:
            expected_outcome_a1_frequency = ex if sign == 1 else 1.0 - ex
            if abs(oy - expected_outcome_a1_frequency) <= EAF_TOLERANCE:
                compatible.append((sign, label))
        signs = {sign for sign, _ in compatible}
        if len(signs) != 1:
            result["harmonization_status"] = "PALINDROMIC_FREQUENCY_UNRESOLVED"
            return result
        sign = next(iter(signs))
        label = next(label for candidate_sign, label in compatible if candidate_sign == sign)
    else:
        signs = {sign for sign, _ in options}
        if len(signs) != 1:
            result["harmonization_status"] = "NONPALINDROMIC_ORIENTATION_CONFLICT"
            return result
        sign = next(iter(signs))
        label = next(label for candidate_sign, label in options if candidate_sign == sign)

    f_stat = (bx / sx) ** 2
    result["f_statistic_recomputed"] = f_stat
    if f_stat < F_MINIMUM:
        result["harmonization_status"] = "WEAK_INSTRUMENT_F_LT_10"
        return result
    result["harmonization_status"] = label
    result["beta_outcome_aligned"] = sign * by
    result["mr_keep"] = True
    return result


def ivw_estimate(rows: pd.DataFrame) -> dict:
    """Wald ratio or multiplicative-random-effects IVW through the origin."""

    k = len(rows)
    if k == 0:
        return {
            "nsnp": 0,
            "method": "TECHNICALLY_INELIGIBLE",
            "beta": np.nan,
            "se": np.nan,
            "z": np.nan,
            "pvalue": 1.0,
            "or": np.nan,
            "or_lci95": np.nan,
            "or_uci95": np.nan,
            "q": np.nan,
            "q_df": 0,
            "q_pvalue": np.nan,
            "mre_scale": np.nan,
        }
    bx = rows["beta_exposure"].astype(float).to_numpy()
    by = rows["beta_outcome_aligned"].astype(float).to_numpy()
    sy = rows["se_outcome"].astype(float).to_numpy()
    if k == 1:
        beta = by[0] / bx[0]
        se = sy[0] / abs(bx[0])
        z = beta / se
        return {
            "nsnp": 1,
            "method": "WALD_RATIO",
            "beta": beta,
            "se": se,
            "z": z,
            "pvalue": 2 * norm.sf(abs(z)),
            "or": math.exp(beta),
            "or_lci95": math.exp(beta - 1.96 * se),
            "or_uci95": math.exp(beta + 1.96 * se),
            "q": np.nan,
            "q_df": 0,
            "q_pvalue": np.nan,
            "mre_scale": 1.0,
        }

    weights = 1.0 / np.square(sy)
    denominator = np.sum(weights * np.square(bx))
    beta = np.sum(weights * bx * by) / denominator
    fixed_se = math.sqrt(1.0 / denominator)
    q = float(np.sum(weights * np.square(by - beta * bx)))
    q_df = k - 1
    scale = math.sqrt(max(1.0, q / q_df))
    se = fixed_se * scale
    z = beta / se
    return {
        "nsnp": k,
        "method": "IVW_MULTIPLICATIVE_RANDOM_EFFECTS",
        "beta": beta,
        "se": se,
        "z": z,
        "pvalue": 2 * norm.sf(abs(z)),
        "or": math.exp(beta),
        "or_lci95": math.exp(beta - 1.96 * se),
        "or_uci95": math.exp(beta + 1.96 * se),
        "q": q,
        "q_df": q_df,
        "q_pvalue": chi2.sf(q, q_df),
        "mre_scale": scale,
    }


def simes(values: list[float]) -> float:
    ordered = np.sort(np.asarray(values, dtype=float))
    if len(ordered) == 0:
        return 1.0
    return float(min(1.0, np.min(ordered * len(ordered) / np.arange(1, len(ordered) + 1))))


def bh(values: list[float]) -> list[float]:
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    output = np.empty_like(adjusted)
    output[order] = adjusted
    return output.tolist()


def load_instruments() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if not WUSTL_INSTRUMENTS.is_file():
        raise FileNotFoundError(WUSTL_INSTRUMENTS)
    if not DECODE_INSTRUMENTS.is_file():
        raise RuntimeError(
            "Official deCODE primary instrument registry is absent; HERMES must remain unopened."
        )
    wustl = pd.read_csv(WUSTL_INSTRUMENTS, sep="\t", dtype={"rsid": str})
    decode = pd.read_csv(DECODE_INSTRUMENTS, sep="\t", dtype={"rsid": str})
    csf = wustl[wustl["source_role"] == PRIMARY_ROLES["CSF"]].copy()
    plasma_primary = decode[decode["source_role"] == PRIMARY_ROLES["PLASMA"]].copy()
    plasma_sensitivity = wustl[wustl["source_role"] == SENSITIVITY_ROLE].copy()
    if set(decode["source_role"]) != {PRIMARY_ROLES["PLASMA"]}:
        raise RuntimeError("deCODE registry contains a non-primary source role")
    return pd.concat([csf, plasma_primary], ignore_index=True), plasma_sensitivity, decode


def load_outcome(rsids: set[str]) -> pd.DataFrame:
    if not rsids:
        raise RuntimeError("No frozen rsID instruments")
    table = pq.read_table(
        HERMES,
        columns=["SNP", "A1", "A2", "BETA", "SE", "P", "EAF", "N"],
        filters=[("SNP", "in", sorted(rsids))],
    )
    outcome = table.to_pandas()
    if outcome["SNP"].duplicated().any():
        duplicates = outcome.loc[outcome["SNP"].duplicated(False), "SNP"].tolist()
        raise RuntimeError(f"Duplicate HERMES rows for frozen instruments: {duplicates[:10]}")
    return outcome


def expected_assays(candidates: pd.DataFrame, roles: dict[str, str]) -> pd.DataFrame:
    rows = []
    for record in candidates.to_dict("records"):
        if record["fluid"] not in roles:
            continue
        for soma_id in str(record["soma_ids"]).split(";"):
            rows.append({**record, "soma_id": soma_id, "source_role": roles[record["fluid"]]})
    return pd.DataFrame(rows)


def run_assays(
    expected: pd.DataFrame,
    instruments: pd.DataFrame,
    outcome: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    outcome_by_snp = outcome.set_index("SNP", drop=False)
    results, harmonized, loo = [], [], []
    grouped = {
        key: frame.copy()
        for key, frame in instruments.groupby(["exposure_id", "soma_id", "source_role"], dropna=False)
    }
    for item in expected.to_dict("records"):
        key = (item["exposure_id"], item["soma_id"], item["source_role"])
        exposure_rows = grouped.get(key, pd.DataFrame())
        kept = []
        if not exposure_rows.empty:
            for _, exposure in exposure_rows.iterrows():
                rsid = str(exposure["rsid"])
                if rsid not in outcome_by_snp.index:
                    row = exposure.to_dict()
                    row.update({"harmonization_status": "OUTCOME_RSID_ABSENT", "mr_keep": False})
                else:
                    row = harmonize_one(exposure, outcome_by_snp.loc[rsid])
                harmonized.append(row)
                if row.get("mr_keep"):
                    kept.append(row)
        estimate = ivw_estimate(pd.DataFrame(kept))
        results.append(
            {
                "exposure_id": item["exposure_id"],
                "protein_id": item["protein_id"],
                "fluid": item["fluid"],
                "protein": item["protein"],
                "gene": item["gene"],
                "soma_id": item["soma_id"],
                "source_role": item["source_role"],
                "frozen_instrument_count": len(exposure_rows),
                "harmonized_instrument_count": len(kept),
                **estimate,
            }
        )
        if len(kept) >= 3:
            kept_frame = pd.DataFrame(kept)
            for omitted in kept_frame["rsid"]:
                subset = kept_frame[kept_frame["rsid"] != omitted]
                loo.append(
                    {
                        "exposure_id": item["exposure_id"],
                        "soma_id": item["soma_id"],
                        "source_role": item["source_role"],
                        "omitted_rsid": omitted,
                        **ivw_estimate(subset),
                    }
                )
    return pd.DataFrame(results), pd.DataFrame(harmonized), pd.DataFrame(loo)


def aggregate_primary(assays: pd.DataFrame, candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    exposure_rows = []
    for exposure in candidates.to_dict("records"):
        subset = assays[assays["exposure_id"] == exposure["exposure_id"]]
        pvalues = subset["pvalue"].fillna(1.0).astype(float).tolist()
        exposure_rows.append(
            {
                "exposure_id": exposure["exposure_id"],
                "protein_id": exposure["protein_id"],
                "fluid": exposure["fluid"],
                "protein": exposure["protein"],
                "gene": exposure["gene"],
                "assay_count": len(subset),
                "eligible_assay_count": int((subset["harmonized_instrument_count"] > 0).sum()),
                "simes_pvalue": simes(pvalues),
            }
        )
    exposures = pd.DataFrame(exposure_rows)
    protein_rows = []
    for protein_id, subset in exposures.groupby("protein_id", sort=True):
        protein_rows.append(
            {
                "protein_id": protein_id,
                "protein": ";".join(sorted(set(subset["protein"]))),
                "gene": ";".join(sorted(set(subset["gene"]))),
                "fluid_specific_exposure_count": len(subset),
                "eligible_fluid_specific_exposure_count": int((subset["eligible_assay_count"] > 0).sum()),
                "simes_pvalue": simes(subset["simes_pvalue"].astype(float).tolist()),
            }
        )
    proteins = pd.DataFrame(protein_rows).sort_values("protein_id").reset_index(drop=True)
    if len(proteins) != 51:
        raise RuntimeError(f"Frozen primary family must contain 51 protein analytes, observed {len(proteins)}")
    proteins["bh_fdr"] = bh(proteins["simes_pvalue"].tolist())
    proteins["hf_fdr_pass"] = proteins["bh_fdr"] < 0.05
    return exposures, proteins


def write_table(frame: pd.DataFrame, name: str) -> Path:
    path = OUTDIR / name
    frame.to_csv(path, sep="\t", index=False, na_rep="NA")
    return path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--execute-frozen-outcome-lookup",
        action="store_true",
        help="required explicit gate; without it the program never opens HERMES",
    )
    parser.add_argument(
        "--self-check",
        action="store_true",
        help="check core MR estimators without accessing provider data",
    )
    args = parser.parse_args()
    if args.self_check:
        synthetic = pd.DataFrame(
            {
                "beta_exposure": [0.20, 0.30],
                "beta_outcome_aligned": [0.04, 0.06],
                "se_outcome": [0.02, 0.03],
            }
        )
        estimate = ivw_estimate(synthetic)
        adjusted = np.asarray(bh([0.001, 0.01, 0.20, 1.0]), dtype=float)
        if estimate["method"] != "IVW_MULTIPLICATIVE_RANDOM_EFFECTS":
            raise RuntimeError("IVW self-check failed")
        if not np.all(np.diff(adjusted) >= 0) or not np.all((0 <= adjusted) & (adjusted <= 1)):
            raise RuntimeError("BH self-check failed")
        print(json.dumps({"state": "PASS", "provider_data_accessed": False}, indent=2))
        return
    candidates = pd.read_csv(CANDIDATES, sep="\t", dtype=str)
    if len(candidates) != 55 or candidates["protein_id"].nunique() != 51:
        raise RuntimeError("Frozen candidate universe changed")
    if not args.execute_frozen_outcome_lookup:
        print(
            json.dumps(
                {
                    "state": "IMPLEMENTATION_VALIDATED_HERMES_UNOPENED",
                    "candidate_exposures": len(candidates),
                    "protein_analytes": candidates["protein_id"].nunique(),
                    "hermes_opened": False,
                },
                indent=2,
            )
        )
        return

    primary, sensitivity, decode = load_instruments()
    all_rsids = set(primary["rsid"].dropna().astype(str)) | set(sensitivity["rsid"].dropna().astype(str))
    outcome = load_outcome(all_rsids)
    primary_expected = expected_assays(candidates, PRIMARY_ROLES)
    sensitivity_expected = expected_assays(
        candidates[candidates["fluid"] == "PLASMA"],
        {"PLASMA": SENSITIVITY_ROLE},
    )
    primary_assays, primary_harmonized, primary_loo = run_assays(primary_expected, primary, outcome)
    sensitivity_assays, sensitivity_harmonized, sensitivity_loo = run_assays(
        sensitivity_expected, sensitivity, outcome
    )
    exposure_results, protein_results = aggregate_primary(primary_assays, candidates)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    outputs = [
        write_table(primary_assays, "primary_assay_mr.tsv"),
        write_table(primary_harmonized, "primary_harmonized_instruments.tsv"),
        write_table(primary_loo, "primary_leave_one_out.tsv"),
        write_table(exposure_results, "primary_fluid_exposure_simes.tsv"),
        write_table(protein_results, "primary_protein_family_bh.tsv"),
        write_table(sensitivity_assays, "wustl_plasma_sensitivity_assay_mr.tsv"),
        write_table(sensitivity_harmonized, "wustl_plasma_sensitivity_harmonized.tsv"),
        write_table(sensitivity_loo, "wustl_plasma_sensitivity_leave_one_out.tsv"),
    ]
    summary = {
        "state": "STAGE_G_HF_MR_COMPLETE_COLOCALIZATION_NOT_OPENED",
        "primary_candidate_exposures": 55,
        "primary_protein_family": 51,
        "primary_fdr_pass_proteins": int(protein_results["hf_fdr_pass"].sum()),
        "hermes_candidate_values_opened": True,
        "colocalization_opened": False,
        "input_sha256": {
            "candidates": sha256(CANDIDATES),
            "wustl_instruments": sha256(WUSTL_INSTRUMENTS),
            "decode_instruments": sha256(DECODE_INSTRUMENTS),
        },
        "output_sha256": {path.name: sha256(path) for path in outputs},
    }
    summary_path = OUTDIR / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
