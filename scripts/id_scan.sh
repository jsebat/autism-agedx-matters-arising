#!/bin/bash
# Run before every push. Fails (exit 1) on study IDs, personal emails, credentials, tracked data files or personal names.
# Study IDs are not PHI, but a public repo stays ID-free (lab rule).
cd "$(dirname "$0")/.."
fail=0
chk() { local lab="$1" pat="$2"; local hits; hits=$(grep -rInE "$pat" --exclude-dir=.git --exclude=id_scan.sh . 2>/dev/null)
  if [ -n "$hits" ]; then echo "== $lab"; echo "$hits" | head -20; fail=1; fi; }
chk "SPARK subject IDs"        '\bSP0[0-9]{6}\b'
chk "SPARK family IDs"         '\bSF0[0-9]{6}\b'
chk "SSC RUIDs"                '\bSS(C)?0[0-9]{4,6}\b'
chk "SSC individual IDs"       '\b1[1-4][0-9]{3}\.(p1|s1|s2|fa|mo)\b'
chk "SDSM barcodes"            'SDSM-[0-9A-Z]{3,}'
chk "REACH IDs"                '\bREACH[0-9]{3,}\b'
chk "email addresses"          '[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.(edu|com|dk|uk|org)'
chk "credentials"              '(api[_-]?key|secret|password|token)[[:space:]]*[:=]'
chk "personal names"           '\b(tran|gendron|stephen|peter|joyce|vinh)\b'
data=$(git ls-files 2>/dev/null | grep -iE '\.(tsv|csv|parquet|txt|xlsx|bed|bim|fam|pgen|psam|pvar|gz|png|pdf)$' | grep -vE '(requirements|README)')
[ -n "$data" ] && { echo "== tracked data files"; echo "$data"; fail=1; }
[ $fail -eq 0 ] && echo "id_scan: clean" || { echo "id_scan: FAILED"; exit 1; }
