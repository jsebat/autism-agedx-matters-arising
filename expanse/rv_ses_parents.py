#!/usr/bin/env python3
"""Tier-1 (high-confidence, constrained-gene) loss-of-function variants in NDD/ASD genes (SFARI score 1 + DDG2P confident):
(A) carrier PARENTS by their own education, family income and neighbourhood deprivation (SPARK legacy trios: WES IDT, WES Twist,
    iWGS v1.1); parent-level logistic regression, family-clustered SEs, adjusted for parent sex and cohort;
(B) per-gene counts of carrier parents (transmitted or not), by parent sex and education;
(C) the same variant pool inherited by autistic children vs their unaffected siblings (sibling families: SPARK WES IDT, Twist, iWGS v1.1,
    2026 WGS, SSC), conditional logistic regression within families; de novo NDD LoF as a positive control.
Inputs: tier-carriers / kids tables (legacy-collected), ~/rvh2/out_dnm_v3 (sibling families), ~/rvh2/agedx/spark_parent_edu.tsv
(SPARK background_history_child 2026-09-28: parents' education, household income), pheno/spark_env.tsv (ADI).
Outputs (aggregate): ~/rvh2/agedx/rv_parent_ses_t1.tsv, rv_parent_genes_t1.tsv, rv_sibling_assoc_t1.tsv"""
import os
import numpy as np, pandas as pd
import statsmodels.api as sm, statsmodels.formula.api as smf
from statsmodels.discrete.conditional_models import ConditionalLogit
A = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08/analyses"; X = A + "/diagnosis-dnm-expanded-20261014/legacy-collected"
GS = A + "/family-genesets-meta-20261014/gene_sets"; P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"; OUT = os.path.expanduser("~/rvh2/agedx")
gs = lambda n: {l.strip() for l in open(f"{GS}/{n}.txt") if l.strip()}
SF1, DDG = gs("sfari_s1"), gs("ddg2p_confident"); NDD = SF1 | DDG
env = pd.read_csv(f"{P}/spark_env.tsv", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")
bh = pd.read_csv(f"{OUT}/spark_parent_edu.tsv", sep="\t", dtype=str).drop_duplicates("subject_sp_id").set_index("subject_sp_id")
DEG = {"baccalaureate_degree", "graduate_or_professional_degree"}
INC = {"less_than_20000": 10, "21000_35000": 28, "36000_50000": 43, "51000_65000": 58, "66000_80000": 73, "81000_100000": 90,
       "101000_130000": 115, "131000_160000": 145, "over_161000": 180}

# ---------------- (A) + (B): parents in SPARK legacy trios
par_rows, var_rows, kid_inh = [], [], []
for ds, cohorts in [("spark_wes", ["spark_wes_idt", "spark_wes_twist"]), ("spark_wgs", ["spark_wgs"])]:
    tc = pd.read_parquet(f"{X}/{ds}/tier-carriers.parquet", columns=["canonical_allele", "sample", "SYMBOL", "tier"])
    tc = tc[tc.tier.eq("lof_t1") & tc.SYMBOL.isin(NDD)].drop_duplicates(["canonical_allele", "sample"])
    for cohort in cohorts:
        k = pd.read_parquet(f"{X}/{ds}/{cohort}.parquet", columns=["sample", "ind_id", "family_id", "mother_sample", "father_sample"]).drop_duplicates("sample")
        kid_alle = set(map(tuple, tc[tc["sample"].isin(set(k["sample"]))][["canonical_allele", "sample"]].to_numpy()))
        for role, col, ecol in [("mother", "mother_sample", "mother_highest_education"), ("father", "father_sample", "father_highest_education")]:
            q = k[["ind_id", "family_id", "sample", col]].rename(columns={col: "parent", "sample": "child"}).dropna(subset=["parent"])
            q["role"] = role; q["cohort"] = cohort
            q["edu"] = q.ind_id.map(bh[ecol]) if ecol in bh else np.nan
            q["inc"] = q.ind_id.map(bh.annual_household_income).map(INC)
            q["adi"] = pd.to_numeric(q.ind_id.map(env.adi_pct), errors="coerce")
            par_rows.append(q)
            v = tc[tc["sample"].isin(set(q.parent))].merge(q[["parent", "child", "role", "cohort", "edu", "family_id"]], left_on="sample", right_on="parent")
            v["transmitted"] = [(a, c) in kid_alle for a, c in zip(v.canonical_allele, v.child)]
            var_rows.append(v)
Pa = pd.concat(par_rows, ignore_index=True)
# one row per parent (a parent of several probands keeps the first proband's family data)
Pa = Pa.sort_values(["parent", "child"]).drop_duplicates("parent").reset_index(drop=True)
V = pd.concat(var_rows, ignore_index=True).drop_duplicates(["canonical_allele", "parent"])
Pa["carrier_lof_t1"] = Pa.parent.isin(set(V["parent"])).astype(int)
Pa["no_degree"] = np.where(Pa.edu.isna() | Pa.edu.eq(""), np.nan, (~Pa.edu.isin(DEG)).astype(float))
q3 = lambda s: pd.qcut(s, 3, labels=False, duplicates="drop")
fam = Pa.drop_duplicates("family_id")
inc_cut = fam.inc.quantile([1 / 3, 2 / 3]).to_numpy(); adi_cut = fam.adi.quantile([1 / 3, 2 / 3]).to_numpy()
Pa["low_income"] = np.where(Pa.inc.isna(), np.nan, np.where(Pa.inc <= inc_cut[0], 1.0, np.where(Pa.inc > inc_cut[1], 0.0, np.nan)))
Pa["high_adi"] = np.where(Pa.adi.isna(), np.nan, np.where(Pa.adi > adi_cut[1], 1.0, np.where(Pa.adi <= adi_cut[0], 0.0, np.nan)))
Pa["income_third"] = np.where(Pa.inc.isna(), np.nan, np.digitize(Pa.inc, inc_cut, right=True))       # 0 low .. 2 high
Pa["adi_third"] = np.where(Pa.adi.isna(), np.nan, np.digitize(Pa.adi, adi_cut, right=True))          # 0 least .. 2 most deprived
Pa["fam"] = pd.factorize(Pa.family_id.astype(str))[0]; Pa["male"] = Pa.role.eq("father").astype(int)
print("parents:", len(Pa), Pa.role.value_counts().to_dict(), "| education recorded:", int(Pa.no_degree.notna().sum()),
      "| income:", int(Pa.inc.notna().sum()), "| ADI:", int(Pa.adi.notna().sum()), "| cuts income", inc_cut, "ADI", adi_cut, flush=True)

rows = []
def logit(outcome, term, data, label, extra=""):
    x = data.dropna(subset=[term]); f = f"{outcome} ~ {term} + male + C(cohort){extra}"
    m = smf.glm(f, data=x, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": x.fam.to_numpy()})
    b, se = m.params[term], m.bse[term]
    g1, g0 = x[x[term] == 1], x[x[term] == 0]
    rows.append(dict(variants=outcome.replace("carrier_", ""), predictor=label, model="alone" if not extra else "joint (all three)", n_parents=int(m.nobs),
                     carriers_exposed=f"{int(g1[outcome].sum())}/{len(g1)} ({100 * g1[outcome].mean():.2f}%)" if len(g1) else "",
                     carriers_reference=f"{int(g0[outcome].sum())}/{len(g0)} ({100 * g0[outcome].mean():.2f}%)" if len(g0) else "",
                     OR=np.exp(b), lo95=np.exp(b - 1.96 * se), hi95=np.exp(b + 1.96 * se), p=m.pvalues[term]))
for outcome in ["carrier_lof_t1"]:
    logit(outcome, "no_degree", Pa, "own education: no degree vs degree")
    logit(outcome, "low_income", Pa, "family income: bottom vs top third")
    logit(outcome, "high_adi", Pa, "deprivation (ADI): most vs least deprived third")
    # trend per third, and a joint model with all three (thirds as ordinal, education binary)
    for term, lab in [("income_third", "family income, per third (higher)"), ("adi_third", "deprivation, per third (more deprived)")]:
        x = Pa.dropna(subset=[term]); m = smf.glm(f"{outcome} ~ {term} + male + C(cohort)", data=x, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": x.fam.to_numpy()})
        rows.append(dict(variants=outcome.replace("carrier_", ""), predictor=lab, model="trend", n_parents=int(m.nobs), OR=np.exp(m.params[term]),
                         lo95=np.exp(m.params[term] - 1.96 * m.bse[term]), hi95=np.exp(m.params[term] + 1.96 * m.bse[term]), p=m.pvalues[term]))
    x = Pa.dropna(subset=["no_degree", "income_third", "adi_third"])
    m = smf.glm(f"{outcome} ~ no_degree + income_third + adi_third + male + C(cohort)", data=x, family=sm.families.Binomial()).fit(cov_type="cluster", cov_kwds={"groups": x.fam.to_numpy()})
    for term, lab in [("no_degree", "own education: no degree"), ("income_third", "family income, per third (higher)"), ("adi_third", "deprivation, per third (more deprived)")]:
        rows.append(dict(variants=outcome.replace("carrier_", ""), predictor=lab, model="joint (all three)", n_parents=int(m.nobs), OR=np.exp(m.params[term]),
                         lo95=np.exp(m.params[term] - 1.96 * m.bse[term]), hi95=np.exp(m.params[term] + 1.96 * m.bse[term]), p=m.pvalues[term]))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/rv_parent_ses_t1.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(R.round(4).to_string(index=False), flush=True)

# (B) genes: carrier parents
V = V.merge(Pa[["parent", "no_degree"]], on="parent", how="left")
G = V.groupby(["SYMBOL"]).agg(carrier_parents=("parent", "nunique"),
                              mothers=("role", lambda s: int((s == "mother").sum())), fathers=("role", lambda s: int((s == "father").sum())),
                              transmitted=("transmitted", "sum"),
                              parent_no_degree=("no_degree", lambda s: int((s == 1).sum())), parent_degree=("no_degree", lambda s: int((s == 0).sum())))
G["sfari_s1"] = G.index.isin(SF1); G["ddg2p_confident"] = G.index.isin(DDG)
G = G.sort_values("carrier_parents", ascending=False); G.to_csv(f"{OUT}/rv_parent_genes_t1.tsv", sep="\t")
print(G.head(40).to_string(), flush=True)

# ---------------- (C) autistic children vs unaffected siblings, same families
S = pd.concat([pd.read_parquet(os.path.expanduser(f"~/rvh2/out_dnm_v3/{c}.parquet")) for c in ["spark_wes_idt", "spark_wes_twist", "spark_wgs", "spark_wgs_new", "ssc"]], ignore_index=True)
S["fam"] = pd.factorize(S.cohort.astype(str) + "_" + S.family_key.astype(str))[0]
crow = []
for lab, cols in [("inherited LoF tier 1, NDD genes", ["lof_t1__inh__ndd"]), ("de novo LoF tier 1, NDD genes (control)", ["lof_t1__dnm__ndd"])]:
    S["carrier"] = (S[cols].sum(axis=1) > 0).astype(int)
    for coh, x in [("pooled", S)] + list(S.groupby("cohort")):
        info = x.groupby("fam").agg(c=("case", "sum"), n=("case", "size")); keep = info[(info.c > 0) & (info.c < info.n)].index
        y = x[x.fam.isin(keep)]
        try:
            m = ConditionalLogit(y.case.to_numpy(), y[["carrier"]].to_numpy(float), groups=y.fam.to_numpy()).fit(disp=0)
            b, se, p = m.params[0], m.bse[0], m.pvalues[0]
        except Exception as e:
            b = se = p = np.nan
        crow.append(dict(variants=lab, cohort=coh, families=len(keep), autistic=int(y.case.sum()), siblings=int((y.case == 0).sum()),
                         carriers_autistic=f"{int(y.loc[y.case == 1, 'carrier'].sum())} ({100 * y.loc[y.case == 1, 'carrier'].mean():.2f}%)",
                         carriers_siblings=f"{int(y.loc[y.case == 0, 'carrier'].sum())} ({100 * y.loc[y.case == 0, 'carrier'].mean():.2f}%)",
                         OR=np.exp(b), lo95=np.exp(b - 1.96 * se), hi95=np.exp(b + 1.96 * se), p=p))
C = pd.DataFrame(crow); C.to_csv(f"{OUT}/rv_sibling_assoc_t1.tsv", sep="\t", index=False)
print(C.round(4).to_string(index=False), flush=True)
