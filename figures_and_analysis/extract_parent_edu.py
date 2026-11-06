"""Extract parental education and household income from the SPARK background_history_child table (release 2026-09-28), one row per
child, for rv_ses_parents.py (copy the output to ~/rvh2/agedx/spark_parent_edu.tsv on the cluster). Standard library only.
usage: python3 extract_parent_edu.py background_history_child-2026-09-28.csv spark_parent_edu.tsv"""
import csv, sys
src, dst = sys.argv[1:3]
keep = ["subject_sp_id", "biomother_sp_id", "biofather_sp_id", "mother_highest_education", "father_highest_education", "annual_household_income"]
seen = set(); n = 0
with open(src, newline="", encoding="utf-8-sig") as f, open(dst, "w", newline="") as g:
    r = csv.DictReader(f); w = csv.writer(g, delimiter="\t"); w.writerow(keep)
    for row in r:
        if row["subject_sp_id"] in seen: continue
        seen.add(row["subject_sp_id"]); w.writerow([row.get(k, "") for k in keep]); n += 1
print("rows", n)
