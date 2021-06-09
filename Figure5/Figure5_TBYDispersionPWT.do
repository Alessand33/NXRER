
use PWT90.dta, clear
/*
gen nxy = csh_x + csh_m
gen try = csh_x - csh_m
gen xm = nxy/(try)
*/
tabstat xm nxy, by(year) stats(v)
/*gen lnpop = ln(pop)*/

kdensity lnpop
tabstat xm nxy if pop>1, by(year) stats(v)
tabstat xm nxy if pop>2, by(year) stats(v)
tabstat xm nxy, by(year) stats(iqr)
tabstat xm nxy, by(year) stats(iqr count)
tabstat xm nxy, by(year) stats(var)
tabstat xm nxy if abs(nxy)<1, by(year) stats(var)
/* sum xm nxy if abs(nxy)>1, by(year) stats(count) */

tabstat xm nxy if abs(nxy)>1, by(year) stats(count)
tabstat xm nxy if abs(nxy)>0.5, by(year) stats(count)
list year country xm nxy if abs(nxy)>0.5
list year country xm nxy if abs(nxy)>0.5
list year country pop xm nxy if abs(nxy)>0.5
list year country pop xm nxy if abs(nxy)>0.5 & nxy!=.
tabstat xm nxy if pop>5, by(year) stats(v)
tabstat try if pop>5, by(year) stats(mean count)
tabstat try if pop>5, by(year) stats(median mean count)
/* gen dy = ln(rgdpna/l.rgdpna) */

/* egen id = group(country)*/
tsset year id
tsset id year
/* gen dy = ln(rgdpna/l.rgdpna)*/
tabstat xm nxy dy if pop>5, by(year) stats(v)
drop if try<0

egen NXYvar5 = sd(nxy) if pop>5, by(year)
egen ct5 = count(nxy) if pop>5, by(year)
egen XMsd5 = sd(xm) if pop>5, by(year)
egen try5 = mean(try) if pop>5, by(year)
egen try5median = median(try) if pop>5, by(year)
egen dysd5 = sd(dy) if pop>5, by(year)
/*
twoway scatter NXYvar5 try5median if year>=1970||lfit NXYvar5 try5median if year>=1970
twoway scatter XMsd5 try5median if year>=1970, yscale(log) xscale(log)||lfit XMsd5 try5median if year>=1970
*/
twoway (scatter NXYvar5 try5median if year>=1970)|| /// 
(scatter XMsd5 try5median if year>=1970, yscale(log) xscale(log))|| ///
(lfit NXYvar5 try5median if year>=1970) || ///
(lfit XMsd5 try5median if year>=1970, yscale(log) xscale(log))
 
