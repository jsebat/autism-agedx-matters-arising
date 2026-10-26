#!/usr/bin/env python3
"""Coverage of the 34 characteristics in the main polygenic CCA sample (registered before 22). Aggregate output only."""
import os
import numpy as np, pandas as pd
OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"
SUB = os.environ.get("SUBSET", "all")
A = pd.read_parquet(f"{OUT}/full_probands.parquet")
pro = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR")].copy()
regage = pd.to_numeric(pro.sample_id.map(pd.read_csv(f"{P}/spark_registration.tsv", sep="	").drop_duplicates("subject_sp_id").set_index("subject_sp_id").age_at_registration_years), errors="coerce")
u22 = regage.lt(22)   # Zhang et al. restrict their within-case age-at-diagnosis analyses to age_at_registration_years < 22
pro = {"all": pro, "zhang": pro[pro.zhang_era], "new": pro[~pro.zhang_era], "p0": pro[pro.n_parents_genotyped.eq(0)],
       "p1": pro[pro.n_parents_genotyped.eq(1)], "p2": pro[pro.n_parents_genotyped.eq(2)], "u22": pro[u22], "p0u22": pro[u22 & pro.n_parents_genotyped.eq(0)],
       "zhangu22": pro[u22 & pro.zhang_era], **{f"src_{k}": pro[pro.source.eq(k)] for k in ["spark_gsa_v1", "spark_gsa_v2", "spark_twist", "spark_wgs", "spark_wgs_2026_08"]}}[SUB].reset_index(drop=True)
PG = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
WITHIN = os.environ.get("WITHIN", "0") == "1"   # within families: pTDT = child score minus mean of both parents' scores (both parents genotyped)
if WITHIN:
    allA = pd.read_parquet(f"{OUT}/full_probands.parquet").drop_duplicates("sample_id").set_index("sample_id")
    pro = pro[pro.n_parents_genotyped.eq(2)].copy()
    for p_ in PG:
        mid = (pro.mother.map(allA[p_]) + pro.father.map(allA[p_])) / 2
        pro[p_] = pro[p_] - mid
    pro = pro.reset_index(drop=True); SUB = SUB + "_within"
pro = pro.rename(columns={"asd": "asd_2026"}); PG = ["asd_2026"] + PG[1:]
pro["sample"] = pro.sample_id; pro["family_key"] = pro.family.fillna(pro.sample_id)
rd = lambda f: pd.read_csv(f"{P}/{f}", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
dx, sev, ho, env, reg, sv = rd("spark_dx_derived.tsv"), rd("spark_severity.tsv"), rd("spark_heldout.tsv"), rd("spark_env.tsv"), rd("spark_registration.tsv"), rd("spark_survey.tsv")
s = pro["sample"]; num = lambda ser: pd.to_numeric(s.map(ser), errors="coerce")
XC = {"age at diagnosis": num(dx.dx_age_y), "age at registration": num(reg.age_at_registration_years), "registration year": num(reg.registration_year),
      "birth year": num(env.birth_year), "maternal age": num(env.maternal_age), "paternal age": num(env.paternal_age), "deprivation (ADI)": num(env.adi_pct),
      "FSIQ": num(sev.fsiq), "SCQ": num(sev.scq_res), "RBS-R": num(sev.rbsr_res), "SRS": num(sev.srs_res), "Vineland ABC": num(sev.vabc_res),
      "CBCL internalizing": num(sev.cbcl_int_res), "CBCL externalizing": num(sev.cbcl_ext_res), "DCDQ": num(sev.dcdq_res),
      "age walked": num(sev.walked_age_mos), "age first words": num(sev.used_words_age_mos), "prenatal PM2.5": num(env.prenatal_pm25)}
XB = {"female": num(dx.sex.map({"Female": 1, "Male": 0})), "intellectual disability": num(dx.ID), "minimal language": num(dx.minimal_language),
      "language regression": num(env.regression), "ADHD": num(ho.behav_adhd), "mood/anxiety": num(ho.mood_or_anx), "seizures": num(ho.seizures),
      "sleep disorder": num(ho.sleep_dx), "feeding disorder": num(ho.feeding_dx), "multiplex family": num(ho.multiplex), "specialist diagnosis": num(ho.dx_specialist),
      "validity flag": num(dx.validity_flag), "low household income (third)": num(sv.s09),
      "mother without degree": num(sv.s10), "ODD/conduct": num(sv.s29), "macrocephaly": num(sv.s37)}
X = pd.DataFrame({**XC, **XB}, index=pro.index)
cov = X.notna().sum().sort_values()
print("n participants:", len(X)); print((cov.to_frame("n_recorded").assign(pct=lambda d: (100 * d.n_recorded / len(X)).round(1))).to_string())
cov.to_csv(f"{OUT}/characteristic_coverage_u22.tsv", sep="	")
