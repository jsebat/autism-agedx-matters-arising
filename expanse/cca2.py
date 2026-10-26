#!/usr/bin/env python3
"""Canonical correlation analyses for the Matters Arising, with fully held-out evaluation (ccalib.py). Replaces the evaluation part of cca_full.py,
cca_trio.py and cca_ssc.py (same samples, characteristics and adjustment; see those scripts' docstrings) and writes the same output names.

env ANALYSIS = full | trio | ssc
  full: SUBSET (u22 | all | zhangu22 | p0 | p1 | p2 | p0u22 | src_<source> ...), PARTIAL ("" | source | source_family), WITHIN (1 = pTDT differences),
        PCS (1 = also regress genotype PCs 1-10)
  trio: YSET (all | pgs | rare | all_nolof), U22 (1 | 0)
  ssc : COHORT (ssc | spark), YSET (pgs | all)
env TASK = main | null | boot | combine ; CHUNK, NCHUNK (null/boot arrays); NPERM, NBOOT (totals)
Outputs (aggregate only), ~/rvh2/agedx/<tag>_{cv,structure,single,ratio,power,missbias}.tsv; chunks in ~/rvh2/agedx/cca2_chunks/."""
import os, sys, glob, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.expanduser("~/rvh2")); import ccalib as CL

OUT = os.path.expanduser("~/rvh2/agedx"); P = "/expanse/projects/sebat1/jsebat/pgs_rare_h2_drift/pheno"; CH = f"{OUT}/cca2_chunks"; os.makedirs(CH, exist_ok=True)
E = os.environ.get; AN = E("ANALYSIS", "full"); TASK = E("TASK", "main"); K = 4
PG = ["asd", "adhd_2022", "mdd_2025", "bip_2024", "scz_2022", "edu_2022", "iq_2018", "afb_2018", "height_2022"]
RC = ["de novo, NDD genes", "de novo, other genes", "inherited, NDD genes", "inherited, other genes", "inherited LoF tier 1"]
rd = lambda f: pd.read_csv(f"{P}/{f}", sep="\t").drop_duplicates("subject_sp_id").set_index("subject_sp_id")


def spark_chars(s):
    num = lambda ser: pd.to_numeric(s.map(ser), errors="coerce")
    dx, sev, ho, env, reg, sv = rd("spark_dx_derived.tsv"), rd("spark_severity.tsv"), rd("spark_heldout.tsv"), rd("spark_env.tsv"), rd("spark_registration.tsv"), rd("spark_survey.tsv")
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
    X = pd.DataFrame({**XC, **XB}); return X, [c in XC for c in X.columns]


def dummies(df, cols):
    Z = pd.get_dummies(df[cols].astype(str), drop_first=True).to_numpy(float); return np.column_stack([np.ones(len(df)), Z])


def load():
    A = pd.read_parquet(f"{OUT}/full_probands.parquet")
    if AN == "full":
        SUB, PARTIAL, WITHIN, PCS = E("SUBSET", "u22"), E("PARTIAL", ""), E("WITHIN", "0") == "1", E("PCS", "0") == "1"
        pro = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR")].copy()
        regage = pd.to_numeric(pro.sample_id.map(rd("spark_registration.tsv").age_at_registration_years), errors="coerce"); u22 = regage.lt(22)
        sel = {"all": pro.index == pro.index, "zhang": pro.zhang_era, "new": ~pro.zhang_era, "p0": pro.n_parents_genotyped.eq(0), "p1": pro.n_parents_genotyped.eq(1),
               "p2": pro.n_parents_genotyped.eq(2), "u22": u22, "p0u22": u22 & pro.n_parents_genotyped.eq(0), "zhangu22": u22 & pro.zhang_era,
               "arraysu22": u22 & pro.source.isin(["spark_gsa_v1", "spark_gsa_v2", "spark_twist"])}
        sel.update({f"src_{k}": pro.source.eq(k) for k in pro.source.unique()}); sel.update({f"srcu22_{k}": u22 & pro.source.eq(k) for k in pro.source.unique()})
        pro = pro[np.asarray(sel[SUB])].reset_index(drop=True); tag = f"cca_full_{SUB}"
        if WITHIN:
            allA = A.drop_duplicates("sample_id").set_index("sample_id"); pro = pro[pro.n_parents_genotyped.eq(2)].copy()
            for p_ in PG: pro[p_] = pro[p_] - (pro.mother.map(allA[p_]) + pro.father.map(allA[p_])) / 2
            pro = pro.reset_index(drop=True); tag += "_within"
        X, cont = spark_chars(pro.sample_id)
        Yd = pro[PG].rename(columns={"asd": "asd_2026"}); ok = Yd.notna().all(axis=1).to_numpy()
        zc = {"": [], "source": ["source"], "source_family": ["source", "n_parents_genotyped"]}[PARTIAL]
        Z = dummies(pro, zc) if zc else np.ones((len(pro), 1))
        tag += {"": "", "source": "_partial", "source_family": "_partialfam"}[PARTIAL]
        if PCS:
            pcs = pro[[f"PC{i}" for i in range(1, 11)]].apply(pd.to_numeric, errors="coerce")
            # SPARK-WGS-2026_08 has no PCs in full_probands: use the PC-AiR scores computed for its split-GWAS stratum (allstrata_prep.sbatch), SDSM -> SP
            J = "/expanse/projects/sebat1/j3guevar/SPARK-WGS-2026_08"
            meta = pd.read_csv(f"{J}/rare/inherited-benchmark-20261010/SPARK.WGS.2026_08.sample_metadata.tsv", sep="\t", dtype=str).set_index("sample_id").ind_id
            pa = pd.read_csv("/expanse/projects/sebat1/jsebat/cnv_drift/warrier_rep/spark_wgs_2026_08/pcair.sscore", sep="\t", dtype={"#IID": str})
            pa["sp"] = pa["#IID"].map(meta); pa = pa.dropna(subset=["sp"]).drop_duplicates("sp").set_index("sp")
            isnew = pro.source.eq("spark_wgs_2026_08")
            for i in range(1, 11): pcs.loc[isnew, f"PC{i}"] = pro.loc[isnew, "sample_id"].map(pa[f"PC{i}_AVG"]).to_numpy()
            print("PC coverage by source:", pcs.PC1.notna().groupby(pro.source).mean().round(3).to_dict(), flush=True)
            ok = ok & pcs.notna().all(axis=1).to_numpy()
            # PCs were computed per genotype set: regress them within platform (platform x PC interactions)
            src = pd.get_dummies(pro.source.astype(str)).to_numpy(float)
            Z = np.column_stack([Z] + [src * pcs[c].fillna(0).to_numpy()[:, None] for c in pcs.columns]); tag += "_pcs"
        singles = ["age at diagnosis", "deprivation (ADI)", "mother without degree", "low household income (third)", "maternal age", "minimal language", "ADHD", "birth year"] if not WITHIN else []
        D = dict(pro=pro, X=X, cont=cont, Yd=Yd, Z=Z, ok=ok, fam=pd.factorize(pro.family.fillna(pro.sample_id))[0], tag=tag, singles=singles, binfill=0.5)
    elif AN == "trio":
        YSET, U22 = E("YSET", "all"), E("U22", "1") == "1"
        R = pd.read_parquet(f"{OUT}/rare_full_cases.parquet"); R = R[~R.cohort.eq("ssc")]
        if U22: R = R[R.u22]
        Aa = A[A.cohort.eq("SPARK") & A.predicted_ancestry.eq("EUR")].drop_duplicates("sample_id")
        pro = R[["sp", "cohort", "family_key"] + RC].rename(columns={"cohort": "seq_cohort"}).merge(Aa, left_on="sp", right_on="sample_id", how="inner").reset_index(drop=True)
        X, cont = spark_chars(pro.sample_id)
        G = [p.replace("asd", "asd_2026") if p == "asd" else p for p in PG]; pro = pro.rename(columns={"asd": "asd_2026"})
        cols = {"all": G + RC, "pgs": G, "rare": RC, "all_nolof": G + RC[:4], "rare_nolof": RC[:4]}[YSET]
        Yd = pro[cols]; ok = Yd.notna().all(axis=1).to_numpy()
        Z = dummies(pro, ["source", "seq_cohort", "n_parents_genotyped"])
        D = dict(pro=pro, X=X, cont=cont, Yd=Yd, Z=Z, ok=ok, fam=pd.factorize(pro.family_key)[0], tag=f"cca_trio_{YSET}" + ("" if U22 else "_all"),
                 singles=["age at diagnosis"], binfill=0.5)
    else:
        COH, YSET = E("COHORT", "ssc"), E("YSET", "pgs"); num = lambda s: pd.to_numeric(s, errors="coerce"); A = A.drop_duplicates("sample_id")
        if COH == "ssc":
            pro = A[A.cohort.eq("SSC") & A.case & A.predicted_ancestry.eq("EUR")].copy().reset_index(drop=True)
            S = pd.read_csv(f"{P}/ssc_cca_chars.tsv", sep="\t", index_col=0).drop_duplicates("SS_RUID").set_index("SS_RUID")
            hz = pd.read_csv("/expanse/projects/sebat1/s3/data/sebat/SFARI_phenotypes/spark_ssc_phenotypes_child.csv", low_memory=False)
            hz = hz[hz.cohort.str.upper().eq("SSC")].dropna(subset=["repository_id"]).drop_duplicates("repository_id").set_index("repository_id")
            s = pro.sample_id; X = pd.DataFrame({c: num(s.map(S[c])) for c in S.columns})
            X["maternal age"] = num(s.map(hz.mother_age_at_birth)); X["paternal age"] = num(s.map(hz.father_age_at_birth))
            me = s.map(hz.mother_education).astype(str).str.lower()
            X["mother without degree"] = (~me.isin(["baccalaureate", "graduate", "post-graduate", "masters", "doctorate", "professional", "graduate-degree"])).astype(float).where(s.map(hz.mother_education).notna())
            inc = num(s.map(hz.household_income)); X["low household income (third)"] = (inc <= inc.quantile(1 / 3)).astype(float).where(inc.notna())
            Z = np.ones((len(pro), 1))
            if YSET == "all":
                R = pd.read_parquet(f"{OUT}/rare_full_cases.parquet"); R = R[R.cohort.eq("ssc")].set_index("sp")
                for c in RC: pro[c] = pro.sample_id.map(R[c])
        else:
            reg = rd("spark_registration.tsv"); pro = A[A.cohort.eq("SPARK") & A.case & A.predicted_ancestry.eq("EUR")].copy()
            pro = pro[num(pro.sample_id.map(reg.age_at_registration_years)).lt(22)].reset_index(drop=True)
            dx, sev, ho, env, sv = rd("spark_dx_derived.tsv"), rd("spark_severity.tsv"), rd("spark_heldout.tsv"), rd("spark_env.tsv"), rd("spark_survey.tsv")
            s = pro.sample_id; m_ = lambda ser: num(s.map(ser))
            X = pd.DataFrame({"female": m_(dx.sex.map({"Female": 1, "Male": 0})), "FSIQ": m_(sev.fsiq), "intellectual disability": m_(dx.ID), "SRS": m_(sev.srs_res), "RBS-R": m_(sev.rbsr_res),
                              "Vineland ABC": m_(sev.vabc_res), "CBCL internalizing": m_(sev.cbcl_int_res), "CBCL externalizing": m_(sev.cbcl_ext_res), "language regression": m_(env.regression),
                              "seizures": m_(ho.seizures), "age walked": m_(sev.walked_age_mos), "age first words": m_(sev.used_words_age_mos), "maternal age": m_(env.maternal_age),
                              "paternal age": m_(env.paternal_age), "mother without degree": m_(sv.s10), "low household income (third)": m_(sv.s09)})
            Z = dummies(pro, ["source", "n_parents_genotyped"])
        Yd = pro[PG + (RC if YSET == "all" else [])].rename(columns={"asd": "asd_2026"}); ok = Yd.notna().all(axis=1).to_numpy()
        D = dict(pro=pro, X=X, cont=[True] * X.shape[1], Yd=Yd, Z=Z, ok=ok, fam=pd.factorize(pro.family.fillna(pro.sample_id))[0], tag=f"cca_harm_{COH}_{YSET}",
                 singles=[], binfill=None)
    ok = D["ok"]
    D.update(Xn=D["X"].to_numpy(float)[ok], Y=D["Yd"].to_numpy(float)[ok], Z=D["Z"][ok], fam=D["fam"][ok], names=list(D["X"].columns),
             pro=D["pro"][ok].reset_index(drop=True), cont=np.array(D["cont"]))
    return D


D = load()
if E("NOSINGLE", "0") == "1": D["singles"] = []
tag = D["tag"]; X, cont, Y, Z, fam, bf = D["Xn"], D["cont"], D["Y"], D["Z"], D["fam"], D["binfill"]
print(f"{tag}: n = {len(Y)}, families = {len(np.unique(fam))}, X {X.shape}, Y {Y.shape}, Z {Z.shape}", flush=True)
REF = f"{CH}/{tag}_ref.npz"


def single_data(c):
    j = D["names"].index(c); m = np.isfinite(X[:, j])
    return X[m][:, [j]], Y[m], Z[m], fam[m], m


if TASK == "main":
    ff = CL.full_fit(X, cont, Y, Z, fam, K, bf, D["names"])
    ref = ff["sy"][:, :K]; np.savez(REF, ref=ref, V=ff["V"][:, :K])
    r, lams, U, V, _ = CL.nested(X, cont, Y, Z, fam, K, ref=ref, seed=0, binfill=bf)
    st = [dict(component=j + 1, block="characteristics", variable=c, r=float(ff["sx"][i, j])) for j in range(min(K, ff["sx"].shape[1])) for i, c in enumerate(ff["cols"])]
    st += [dict(component=j + 1, block="genetics", variable=c, r=float(ff["sy"][i, j])) for j in range(min(K, ff["sy"].shape[1])) for i, c in enumerate(D["Yd"].columns)]
    pd.DataFrame(st).to_csv(f"{OUT}/{tag}_structure.tsv", sep="\t", index=False)
    json.dump(dict(r=[float(x) for x in r], lams=[float(x) for x in lams], lam_full=ff["lam"], n=int(len(Y)), n_families=int(len(np.unique(fam)))), open(f"{CH}/{tag}_obs.json", "w"))
    sr = []
    for c in D["singles"]:
        xs, ys, zs_, fs_, _ = single_data(c); rr, ll, *_ = CL.nested(xs, np.array([True]), ys, zs_, fs_, 1, seed=0, binfill=bf)
        sr.append(dict(characteristic=c, n=int(len(ys)), n_families=int(len(np.unique(fs_))), cv_r=float(rr[0]), lams=ll)); print(sr[-1], flush=True)
    json.dump(sr, open(f"{CH}/{tag}_single_obs.json", "w"))
    # held-out genetic variates for power / missingness (no axis estimation on the evaluated person)
    if E("POWER", "0") == "1" or E("MISSBIAS", "0") == "1":
        from scipy.stats import norm, chi2 as chi2d
        Xraw = D["X"][D["ok"]].reset_index(drop=True); zcrit = norm.ppf(1 - 0.05 / 2 / Xraw.shape[1]) + norm.ppf(0.80)
        resid = lambda y, z: y - z @ np.linalg.lstsq(z, y, rcond=None)[0]; rng = np.random.default_rng(5)
        # family bootstrap as row weights (one draw = multinomial family counts), precomputed once and shared by all correlations
        _ui, _inv = np.unique(fam, return_inverse=True); WB = np.stack([np.bincount(rng.integers(0, len(_ui), len(_ui)), minlength=len(_ui))[_inv] for _ in range(300)]).astype(float)
        def cboot(x, v, fm, B=300, m=None):
            W = WB[:, m] if m is not None else WB; sw = W.sum(1)
            mx, mv = W @ x / sw, W @ v / sw; cxv = W @ (x * v) / sw - mx * mv; cxx = W @ (x * x) / sw - mx ** 2; cvv = W @ (v * v) / sw - mv ** 2
            z = np.arctanh(np.clip(cxv / np.sqrt(cxx * cvv), -0.999, 0.999)); return float(np.tanh(np.percentile(z, 2.5))), float(np.tanh(np.percentile(z, 97.5))), float(z.std())
        if E("POWER", "0") == "1":
            rows = []
            for c in Xraw.columns:
                m = Xraw[c].notna().to_numpy(); n = int(m.sum())
                if n < 50: continue
                xr = resid(Xraw[c].to_numpy(float)[m], Z[m]); nf = int(len(np.unique(fam[m])))
                for j in range(K):
                    vr = resid(V[m, j], Z[m]); r_ = float(np.corrcoef(xr, vr)[0, 1]); lo, hi, sez = cboot(xr, vr, fam[m], 300, m)
                    rows.append(dict(characteristic=c, axis=j + 1, n=n, n_families=nf, r=r_, lo95=lo, hi95=hi, mde_r_80pct_bonf=float(np.tanh(zcrit * sez))))
            pd.DataFrame(rows).to_csv(f"{OUT}/{tag}_power.tsv", sep="\t", index=False)
        if E("MISSBIAS", "0") == "1":
            fs = D["pro"].n_parents_genotyped.to_numpy(); Zs = dummies(D["pro"], ["source"]); rows = []
            for c in Xraw.columns:
                rec = Xraw[c].notna().to_numpy(); row = dict(characteristic=c, pct_recorded=100 * rec.mean(), **{f"pct_recorded_{k}par": 100 * rec[fs == k].mean() for k in (0, 1, 2)})
                for j in (0, 1): row[f"r_missing_axis{j + 1}"] = float(np.corrcoef(resid((~rec).astype(float), Zs), resid(V[:, j], Zs))[0, 1])
                est = []
                for k in (0, 1, 2):
                    m = rec & (fs == k); n = int(m.sum())
                    if n < 100: row[f"r_axis1_{k}par"] = np.nan; continue
                    xx = resid(Xraw[c].to_numpy(float)[m], Zs[m]); vv = resid(V[m, 0], Zs[m]); r_ = float(np.corrcoef(xx, vv)[0, 1])
                    _, _, sez = cboot(xx, vv, fam[m], 300, m); row[f"r_axis1_{k}par"] = r_; row[f"n_{k}par"] = n; est.append((np.arctanh(r_), sez ** 2))
                if len(est) >= 2:
                    z = np.array([e[0] for e in est]); w = 1 / np.array([e[1] for e in est]); zb = (w * z).sum() / w.sum(); Q = float((w * (z - zb) ** 2).sum())
                    row["Q_het"] = Q; row["p_het"] = float(chi2d.sf(Q, len(est) - 1))
                rows.append(row)
            pd.DataFrame(rows).to_csv(f"{OUT}/{tag}_missbias.tsv", sep="\t", index=False)
    print(tag, "observed r:", np.round(r, 3), "lambdas:", lams, flush=True)

elif TASK in ("null", "boot"):
    ch, nch = int(E("CHUNK", "0")), int(E("NCHUNK", "1")); tot = int(E("NPERM" if TASK == "null" else "NBOOT", "200"))
    reps = [i for i in range(tot) if i % nch == ch]; z = np.load(REF); ref, Vfull = z["ref"], z["V"]; out = []
    for i in reps:
        rng = np.random.default_rng(10_000 * (TASK == "boot") + i + 1)
        if TASK == "null":
            rec = dict(rep=i)
            Yp = CL.fl_perm(Y, Z, fam, rng); r, *_ = CL.nested(X, cont, Yp, Z, fam, K, ref=None, seed=i + 1, binfill=bf); rec["global"] = [float(x) for x in r]
            seq = []
            for j in range(1, K):   # axis j+1 after deflating on the first j full-sample genetic variates
                Za = np.column_stack([Z, Vfull[:, :j]]); Yp = CL.fl_perm(Y, Za, fam, rng); rr, *_ = CL.nested(X, cont, Yp, Za, fam, 1, seed=i + 1, binfill=bf); seq.append(float(rr[0]))
            rec["seq"] = seq; rec["single"] = {}
            for c in D["singles"]:
                xs, ys, zs_, fs_, _ = single_data(c); yp = CL.fl_perm(ys, zs_, fs_, rng)
                rr, *_ = CL.nested(xs, np.array([True]), yp, zs_, fs_, 1, seed=i + 1, binfill=bf); rec["single"][c] = float(rr[0])
        else:
            idx, nf_ = CL.boot_rows(fam, rng); r, *_ = CL.nested(X[idx], cont, Y[idx], Z[idx], nf_, K, ref=ref, seed=i + 1, binfill=bf)
            rec = dict(rep=i, axes=[float(x) for x in r], single={})
            for c in D["singles"]:
                j = D["names"].index(c); m = np.isfinite(X[idx, j])
                rr, *_ = CL.nested(X[idx][m][:, [j]], np.array([True]), Y[idx][m], Z[idx][m], nf_[m], 1, seed=i + 1, binfill=bf); rec["single"][c] = float(rr[0])
        out.append(rec); print(TASK, i, rec, flush=True)
    json.dump(out, open(f"{CH}/{tag}_{TASK}_{ch}.json", "w"))

elif TASK == "combine":
    ob = json.load(open(f"{CH}/{tag}_obs.json")); r = np.array(ob["r"])
    nul = sum([json.load(open(f)) for f in sorted(glob.glob(f"{CH}/{tag}_null_*.json"))], [])
    bt = sum([json.load(open(f)) for f in sorted(glob.glob(f"{CH}/{tag}_boot_*.json"))], [])
    rows = []
    for j in range(K):
        row = dict(component=j + 1, cv_canonical_r=r[j], shared_variance=r[j] ** 2 if np.isfinite(r[j]) else np.nan)
        if nul:
            g = np.array([x["global"][j] for x in nul]); row.update(null_mean=g.mean(), null_95=np.percentile(g, 95))
            if j == 0: row.update(p_perm=(1 + (g >= r[0]).sum()) / (1 + len(g)), test="global")
            else:
                obs_seq = json.load(open(f"{CH}/{tag}_seqobs.json"))[j - 1] if os.path.exists(f"{CH}/{tag}_seqobs.json") else r[j]
                s = np.array([x["seq"][j - 1] for x in nul]); row.update(seq_null_95=np.percentile(s, 95), p_perm=(1 + (s >= obs_seq).sum()) / (1 + len(s)), test="sequential", n_perm=len(s))
        if bt:
            b = np.array([x["axes"][j] for x in bt]); row.update(boot_lo=np.percentile(b, 2.5), boot_hi=np.percentile(b, 97.5), n_boot=len(b))
        rows.append(row)
    CV = pd.DataFrame(rows); CV["lambda_folds"] = ",".join(map(str, ob["lams"])); CV["n"] = ob["n"]; CV["n_families"] = ob["n_families"]
    CV.to_csv(f"{OUT}/{tag}_cv.tsv", sep="\t", index=False); print(CV.round(4).to_string(index=False))
    if os.path.exists(f"{CH}/{tag}_single_obs.json") and D["singles"]:
        so = json.load(open(f"{CH}/{tag}_single_obs.json")); S = []
        for x in so:
            c = x["characteristic"]; row = dict(characteristic=c, n=x["n"], n_families=x["n_families"], cv_r=x["cv_r"], shared_variance=x["cv_r"] ** 2)
            if nul:
                g = np.array([y["single"][c] for y in nul]); row.update(null_95=np.percentile(g, 95), p_perm=(1 + (g >= x["cv_r"]).sum()) / (1 + len(g)))
            if bt:
                b = np.array([y["single"][c] for y in bt]); row.update(boot_lo=np.percentile(b, 2.5), boot_hi=np.percentile(b, 97.5))
            S.append(row)
        S = pd.DataFrame(S); S.to_csv(f"{OUT}/{tag}_single.tsv", sep="\t", index=False); print(S.round(4).to_string(index=False))
        if bt:   # paired fold differences: shared variance of each axis / of each single characteristic, relative to age at diagnosis alone
            ref_c = "age at diagnosis"; RR = []
            aad = np.array([y["single"][ref_c] for y in bt]) ** 2; aad0 = S.set_index("characteristic").loc[ref_c, "shared_variance"]
            for j in range(K):
                b = np.array([y["axes"][j] for y in bt]) ** 2 / aad
                RR.append(dict(item=f"axis {j + 1}", fold=r[j] ** 2 / aad0, lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5)))
            for x in so:
                b = np.array([y["single"][x["characteristic"]] for y in bt]) ** 2 / aad
                RR.append(dict(item=x["characteristic"], fold=x["cv_r"] ** 2 / aad0, lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5)))
            pd.DataFrame(RR).to_csv(f"{OUT}/{tag}_ratio.tsv", sep="\t", index=False); print(pd.DataFrame(RR).round(2).to_string(index=False))

if TASK == "main" and os.environ.get("SEQOBS", "1") == "1":
    # observed sequential statistics (axis j+1 after deflating on the first j full-sample genetic variates)
    Vfull = np.load(REF)["V"]; so = []
    for j in range(1, K):
        Za = np.column_stack([Z, Vfull[:, :j]]); rr, *_ = CL.nested(X, cont, Y, Za, fam, 1, seed=0, binfill=bf); so.append(float(rr[0]))
    json.dump(so, open(f"{CH}/{tag}_seqobs.json", "w")); print("sequential observed:", np.round(so, 3), flush=True)
