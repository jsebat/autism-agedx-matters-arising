#!/usr/bin/env python3
"""Broad survey of case characteristics (Matters Arising on Zhang, Grove et al. 2025): build split codes for every SPARK autistic
participant from the canonical SPARK release (2026-09-28) plus the derived tables used throughout (Expanse pheno/). Each split is
1 (high / yes), 0 (low / no) or missing; continuous measures split top vs bottom third over all SPARK autistic participants.
Writes spark_survey.tsv (participant-level; copied to Expanse pgs_rare_h2_drift/pheno only) and, to the report folder, an
aggregate table survey_splits.tsv (n per half, phi with late diagnosis (>10 vs <6 y) and with adult registration)."""
import sys, numpy as np, pandas as pd
D = "C:/Users/Jonathan Sebat/sebatlab Dropbox/Autism/SPARK Phenotype Data/SPARKDataRelease_2026-09-28"; DER = sys.argv[1]; OUTP = sys.argv[2]; OUTA = sys.argv[3]
rd = lambda t, cols: pd.read_csv(f"{D}/{t}-2026-09-28.csv", usecols=lambda c: c in cols, low_memory=False).drop_duplicates("subject_sp_id").set_index("subject_sp_id")
reg = rd("individuals_registration", ["subject_sp_id", "asd", "sex", "age_at_registration_years", "registration_year", "diagnosis_age"])
reg = reg[reg.asd.astype(str).str.upper().eq("TRUE")]
core = rd("core_descriptive_variables", ["subject_sp_id", "cognitive_impairment_latest", "language_level_latest", "age_onset_mos", "current_depend_adult"])
bms = rd("basic_medical_screening", ["subject_sp_id", "birth_prem", "growth_macroceph", "growth_microceph", "behav_odd", "behav_conduct", "tics", "dev_motor", "dev_ld", "growth_obes"])
bh = rd("background_history_child", ["subject_sp_id", "annual_household_income", "mother_highest_education"])
der = {k: pd.read_csv(f"{DER}/spark_{k}.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id") for k in ["dx_derived", "env", "heldout", "severity"]}
ix = reg.index; g = lambda df, c: pd.to_numeric(df[c].reindex(ix), errors="coerce") if c in df.columns else pd.Series(np.nan, index=ix)
def tert(v): lo, hi = v.quantile(1 / 3), v.quantile(2 / 3); return pd.Series(np.where(v >= hi, 1.0, np.where(v <= lo, 0.0, np.nan)), index=ix)
def binm(df, c):  # medical screening items: 1 if endorsed; 0 if the participant completed the form
    v = pd.to_numeric(df[c], errors="coerce").reindex(ix); done = ix.isin(df.index); return pd.Series(np.where(v.fillna(0) > 0, 1.0, np.where(done, 0.0, np.nan)), index=ix)
dxm = g(reg, "diagnosis_age"); areg = g(reg, "age_at_registration_years"); dy = g(der["dx_derived"], "dx_year"); ry = g(reg, "registration_year")
inc = pd.to_numeric(bh.annual_household_income.reindex(ix).astype(str).str.extract(r"(\d[\d,]*)")[0].str.replace(",", ""), errors="coerce")
medu = bh.mother_highest_education.reindex(ix).astype(str).str.lower()
S = {
 "agedx: diagnosed >10 vs <6 y": pd.Series(np.where(dxm > 120, 1.0, np.where(dxm < 72, 0.0, np.nan)), index=ix),
 "design: registered as adult": (areg >= 18).astype(float).where(areg.notna()),
 "design: registered 2021+ vs before 2019": pd.Series(np.where(ry >= 2021, 1.0, np.where(ry < 2019, 0.0, np.nan)), index=ix),
 "design: diagnosed 2019+ vs before 2010": pd.Series(np.where(dy >= 2019, 1.0, np.where(dy < 2010, 0.0, np.nan)), index=ix),
 "design: multiplex family": g(der["heldout"], "multiplex"),
 "design: diagnosed by specialist": g(der["heldout"], "dx_specialist"),
 "design: SPARK validity flag": g(der["dx_derived"], "validity_flag"),
 "demog: female": pd.Series(np.where(reg.sex.eq("Female"), 1.0, np.where(reg.sex.eq("Male"), 0.0, np.nan)), index=ix),
 "SES: area deprivation (top vs bottom third)": tert(g(der["env"], "adi_pct")),
 "SES: household income (bottom vs top third)": 1 - tert(inc),
 "SES: mother without degree vs bachelor+": medu.map({"baccalaureate_degree": 0.0, "graduate_or_professional_degree": 0.0, "some_college": 1.0, "associate_degree": 1.0, "high_school_graduate": 1.0, "trade_school": 1.0, "ged_diploma": 1.0, "some_high_school": 1.0, "did_not_attend_high_school": 1.0}),
 "SES: maternal age (top vs bottom third)": tert(g(der["env"], "maternal_age")),
 "SES: paternal age (top vs bottom third)": tert(g(der["env"], "paternal_age")),
 "cog: intellectual disability": g(der["dx_derived"], "ID"),
 "cog: cognitive impairment (latest)": core.cognitive_impairment_latest.reindex(ix).astype(str).map({"True": 1.0, "False": 0.0}),
 "cog: minimal language": g(der["severity"], "minimal_language"),
 "dev: walked late (top vs bottom third)": tert(g(der["severity"], "walked_age_mos")),
 "dev: first words late (top vs bottom third)": tert(g(der["severity"], "used_words_age_mos")),
 "dev: language regression": g(der["env"], "regression"),
 "dev: early onset concern (bottom vs top third age)": 1 - tert(g(core, "age_onset_mos")),
 "dev: motor coordination problems (DCDQ)": 1 - tert(g(der["severity"], "dcdq_res")),
 "dev: motor delay": binm(bms, "dev_motor"),
 "dev: learning disability": binm(bms, "dev_ld"),
 "sev: SCQ (top vs bottom third)": tert(g(der["severity"], "scq_res")),
 "sev: RBS-R (top vs bottom third)": tert(g(der["severity"], "rbsr_res")),
 "sev: SRS (top vs bottom third)": tert(g(der["severity"], "srs_res")),
 "sev: Vineland ABC (bottom vs top third)": 1 - tert(g(der["severity"], "vabc_res")),
 "psy: ADHD": g(der["heldout"], "behav_adhd"),
 "psy: mood or anxiety disorder": g(der["heldout"], "mood_or_anx"),
 "psy: ODD or conduct disorder": pd.Series(np.fmax(binm(bms, "behav_odd"), binm(bms, "behav_conduct")), index=ix),
 "psy: tics": binm(bms, "tics"),
 "psy: CBCL internalizing (top vs bottom third)": tert(g(der["severity"], "cbcl_int_res")),
 "psy: CBCL externalizing (top vs bottom third)": tert(g(der["severity"], "cbcl_ext_res")),
 "med: seizures": g(der["heldout"], "seizures"),
 "med: sleep disorder": g(der["heldout"], "sleep_dx"),
 "med: feeding disorder": g(der["heldout"], "feeding_dx"),
 "med: preterm birth": binm(bms, "birth_prem"),
 "med: macrocephaly": binm(bms, "growth_macroceph"),
 "med: obesity": binm(bms, "growth_obes"),
 "med: genetic diagnosis reported": g(der["heldout"], "genetic_dx_reported"),
 "env: prenatal PM2.5 (top vs bottom third)": tert(g(der["env"], "prenatal_pm25")),
 "env: prenatal NO2 (top vs bottom third)": tert(g(der["env"], "prenatal_no2")),
 "design: born before 2005 vs 2012+": pd.Series(np.where(g(der["env"], "birth_year") < 2005, 1.0, np.where(g(der["env"], "birth_year") >= 2012, 0.0, np.nan)), index=ix),
}
P = pd.DataFrame(S); P.columns = [f"s{i:02d}" for i in range(len(S))]; P.index.name = "subject_sp_id"
P.reset_index().to_csv(OUTP, sep="\t", index=False)
late = S["agedx: diagnosed >10 vs <6 y"]; adult = S["design: registered as adult"]
rows = []
for i, (k, v) in enumerate(S.items()):
    ok = v.notna() & late.notna(); ok2 = v.notna() & adult.notna()
    rows.append(dict(code=f"s{i:02d}", split=k, n_high=int((v == 1).sum()), n_low=int((v == 0).sum()),
                     phi_late_dx=float(np.corrcoef(v[ok], late[ok])[0, 1]) if ok.sum() > 100 else np.nan,
                     phi_adult_reg=float(np.corrcoef(v[ok2], adult[ok2])[0, 1]) if ok2.sum() > 100 else np.nan))
A = pd.DataFrame(rows); A.to_csv(OUTA, sep="\t", index=False); pd.set_option("display.width", 200); print(A.round(3).to_string(index=False))
