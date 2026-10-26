#!/usr/bin/env python3
"""Supplementary Table: sample sizes by subcohort and role for every analysis .
Rows = subcohort (genotyping platform / sequencing cohort); columns = autistic participants (all ages; registered before 22), their parents and
unaffected siblings with polygenic scores, and the rare-variant trio probands with their sequenced parents and the sibling reference.
Aggregate output only: ~/rvh2/agedx/supp_sample_table.tsv"""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
A = pd.read_parquet(f"{OUT}/full_probands.parquet"); A = A[A.cohort.eq("SPARK")]
reg = pd.read_csv(f"{P}/spark_registration.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
A["u22"] = pd.to_numeric(A.sample_id.map(reg.age_at_registration_years), errors="coerce").lt(22)
eur = A.predicted_ancestry.eq("EUR")
ped = pd.read_csv(f"{P}/spark_pedigree_ids.tsv", sep="\t", dtype=str)
parents = set(ped.biomother_sp_id.dropna()) | set(ped.biofather_sp_id.dropna())
A["role"] = np.where(A.case, "autistic", np.where(A.sample_id.isin(parents), "parent", "sibling or other relative"))
rows = []
for src, g in A[eur].groupby("source"):
    c = g[g.role.eq("autistic")]
    rows.append(dict(subcohort=src, autistic_all_ages=len(c), autistic_registered_before_22=int(c.u22.sum()),
                     autistic_u22_0_parents=int((c.u22 & c.n_parents_genotyped.eq(0)).sum()), autistic_u22_1_parent=int((c.u22 & c.n_parents_genotyped.eq(1)).sum()),
                     autistic_u22_2_parents=int((c.u22 & c.n_parents_genotyped.eq(2)).sum()),
                     parents_with_PGS=int(g.role.eq("parent").sum()), non_autistic_siblings_relatives_with_PGS=int(g.role.eq("sibling or other relative").sum())))
T = pd.DataFrame(rows)
R = pd.read_parquet(f"{OUT}/rare_full_cases.parquet")
sib = {c: int((pd.read_parquet(os.path.expanduser(f"~/rvh2/out_dnm_v3/{c}.parquet")).case == 0).sum()) for c in ["spark_wes_idt", "spark_wes_twist", "spark_wgs", "spark_wgs_new", "ssc"]}
rr = R.groupby("cohort").agg(trio_probands_all_ages=("sp", "size"), trio_probands_u22=("u22", "sum")).reset_index().rename(columns={"cohort": "sequencing_cohort"})
rr["sequenced_parents_u22"] = 2 * rr.trio_probands_u22; rr["unaffected_sibling_reference"] = rr.sequencing_cohort.map(sib)
T.to_csv(f"{OUT}/supp_sample_table_pgs.tsv", sep="\t", index=False); rr.to_csv(f"{OUT}/supp_sample_table_rare.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); print(T.to_string(index=False)); print(T.drop(columns="subcohort").sum().to_dict()); print(rr.to_string(index=False)); print(rr.drop(columns="sequencing_cohort").sum().to_dict())
