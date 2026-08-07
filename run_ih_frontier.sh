#!/bin/bash
# In-house frontier, phase 2 onward: chunked greedy deletion test over the
# level-20 extension corpus, C escalation of survivors, then the level-21
# iteration. Run after ext20_ih_raw.txt is complete.
set -e
N=${1:-8}
EXT=./ext6b.exe
split -d -n l/$N ext20_ih_raw.txt ihchunk_
for f in ihchunk_??; do
  $EXT frontier "$f" residue_inhouse.txt > "$f.out" 2> "$f.log" &
done
wait
cat ihchunk_??.out > frontier20_ih.txt
grep -h "^SURVIVES" frontier20_ih.txt | sed 's/^SURVIVES [0-9]* //' > f20_ih_surv.txt
echo "greedy survivors: $(wc -l < f20_ih_surv.txt)"
$EXT escalate f20_ih_surv.txt residue_inhouse.txt > f20_ih_esc_out.txt 2> f20_ih_esc.log
tail -1 f20_ih_esc.log
grep "^SURVIVES" f20_ih_esc_out.txt | sed 's/^SURVIVES [0-9]* //' > frontier20_ih_true.txt
echo "in-house level-20 true frontier: $(wc -l < frontier20_ih_true.txt)"
$EXT extend frontier20_ih_true.txt - prov21_ih.txt > ext21_ih_raw.txt 2> ext21_ih.log
tail -1 ext21_ih.log
$EXT frontier ext21_ih_raw.txt frontier20_ih_true.txt > frontier21_ih.txt 2> frontier21_ih.log
tail -1 frontier21_ih.log
grep "^SURVIVES" frontier21_ih.txt | sed 's/^SURVIVES [0-9]* //' > f21_ih_surv.txt
if [ -s f21_ih_surv.txt ]; then
  $EXT escalate f21_ih_surv.txt frontier20_ih_true.txt > f21_ih_esc_out.txt 2> f21_ih_esc.log
  tail -1 f21_ih_esc.log
  grep "^SURVIVES" f21_ih_esc_out.txt | sed 's/^SURVIVES [0-9]* //' > frontier21_ih_true.txt
  echo "in-house level-21 true frontier: $(wc -l < frontier21_ih_true.txt)"
else
  echo "in-house level-21 frontier EMPTY (greedy alone)"
fi
rm -f ihchunk_??
