clear all
use "temp/df.dta"
file open tsv using "temp/stata_out.tsv", write replace

timer clear 1
timer on 1
quietly robreg m y x z i.t, ivar(i) cluster(i) eff(95)
timer off 1
timer list 1
scalar run_time = r(t1)

file write tsv "run_time" _tab "est" _tab "se" _tab "scale" _n
file write tsv %21.0g (run_time) _tab %21.0g (_b[x]) _tab %21.0g (_se[x]) _tab %21.0g (e(scale)) _n
file close tsv
