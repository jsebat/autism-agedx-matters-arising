#!/usr/bin/env python3
"""Rare-variant cells for ALL legacy SPARK trio probands (sibling not required), for the age-at-diagnosis Matters Arising.
Source: the lab's diagnosis-dnm-expanded-20261014/legacy-collected/{spark_wes,spark_wgs}: trio-proven autistic children
({cohort}.parquet), tier-qualifying carrier rows for children and parents (tier-carriers.parquet) and de novo calls (DNMs.parquet).
Cells: 6 tiers x {dnm, inh} x {ndd, nonndd}, same definitions as build_dnm.py / out_dnm_v3:
  dnm = tier-qualifying allele of the child that is in the child's de novo call set; inh = every other tier-qualifying allele of the child
  NDD = SFARI score 1 + DDG2P confident (the lab's gene_sets); allele is NDD if any annotated SYMBOL is in the set. SNV/indel only (no CNV).
SSC and the new WGS release are not rebuilt here: their cells stay those of ~/rvh2/out_dnm_v3 (sibling families; new-WGS v3 strand filter).
Validation: children present in out_dnm_v3 must match cell for cell.
Output (participant-level, Expanse only): ~/rvh2/out_trio/{cohort}.parquet; receipt ~/rvh2/out_trio/receipt.json"""
import json, os
import numpy as np, pandas as pd
A = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08/analyses"; X = A + "/diagnosis-dnm-expanded-20261014/legacy-collected"
GS = A + "/family-genesets-meta-20261014/gene_sets"; V3 = os.path.expanduser("~/rvh2/out_dnm_v3"); OUT = os.path.expanduser("~/rvh2/out_trio"); os.makedirs(OUT, exist_ok=True)
TIERS = ["lof_t1", "lof_t2", "miss_t1", "miss_t2", "miss_t3", "miss_t4"]; CELLS = [f"{t}__{m}__{g}" for t in TIERS for m in ("dnm", "inh") for g in ("ndd", "nonndd")]
gs = lambda n: {l.strip() for l in open(f"{GS}/{n}.txt") if l.strip()}
NDD = gs("sfari_s1") | gs("ddg2p_confident")
rec = {}
for ds, cohorts in [("spark_wes", ["spark_wes_idt", "spark_wes_twist"]), ("spark_wgs", ["spark_wgs"])]:
    tc = pd.read_parquet(f"{X}/{ds}/tier-carriers.parquet", columns=["canonical_allele", "sample", "SYMBOL", "tier"])
    dn = pd.read_parquet(f"{X}/{ds}/DNMs.parquet", columns=["canonical_allele", "sample"]).drop_duplicates()
    dn["is_dnm"] = True
    for cohort in cohorts:
        kids = pd.read_parquet(f"{X}/{ds}/{cohort}.parquet")
        if cohort == cohorts[0]: rec[f"{ds}_columns"] = list(kids.columns)
        ks = set(kids["sample"])
        t = tc[tc["sample"].isin(ks)].copy(); t["ndd"] = t.SYMBOL.isin(NDD)
        tiers_per = t.groupby(["canonical_allele", "sample"]).tier.nunique()
        al = t.groupby(["canonical_allele", "sample"]).agg(tier=("tier", "first"), ndd=("ndd", "max")).reset_index()
        al = al.merge(dn, on=["canonical_allele", "sample"], how="left"); al["mode"] = np.where(al.is_dnm.fillna(False).astype(bool), "dnm", "inh")
        out = pd.DataFrame(index=pd.Index(sorted(ks), name="sample"))
        for c in CELLS:
            tt, md, g = c.split("__"); y = al[(al.tier == tt) & (al["mode"] == md) & (al.ndd == (g == "ndd"))]
            out[c] = y.groupby("sample").size()
        out = out.fillna(0).astype(int)
        keep = [c for c in ["sample", "ind_id", "family_id", "male", "diagnosis_age_years", "registration_age_years", "registration_year", "mother_sample", "father_sample"] + [f"PC{i}" for i in range(1, 11)] if c in kids.columns]
        x = kids[keep].drop_duplicates("sample").set_index("sample").join(out).reset_index(); x["cohort"] = cohort; x["case"] = 1
        x.to_parquet(f"{OUT}/{cohort}.parquet", index=False)
        # validation against the sibling-family cells
        v = pd.read_parquet(f"{V3}/{cohort}.parquet").set_index("sample")[CELLS]
        both = sorted(set(v.index) & set(out.index))
        diff = (out.loc[both, CELLS] - v.loc[both, CELLS]).abs()
        rec[cohort] = dict(children=len(x), alleles=len(al), alleles_multi_tier=int((tiers_per > 1).sum()), dnm_alleles=int((al["mode"] == "dnm").sum()),
                           dnm_calls_not_tier_qualifying=int((~dn[dn["sample"].isin(ks)].set_index(["canonical_allele", "sample"]).index.isin(al.set_index(["canonical_allele", "sample"]).index)).sum()),
                           validate_children=len(both), validate_cells_diff_total=int(diff.to_numpy().sum()),
                           validate_by_cell={c: int(diff[c].sum()) for c in CELLS if diff[c].sum()},
                           mean_counts={c: round(float(out[c].mean()), 4) for c in CELLS[:8]})
        print(cohort, json.dumps(rec[cohort])[:1500], flush=True)
json.dump(rec, open(f"{OUT}/receipt.json", "w"), indent=1)
