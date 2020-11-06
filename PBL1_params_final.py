# all parameters should be in nM for conc and hr for time
kb = 30
kd = 8.6e-11
karna = 0.055
kdrna = 1.38e-8
drna = 8.4
dkcat = 104.4
kmrnase = 163
en = 2.44e3
bt = 10**4
bs1 = 1
bs2 = 1
bs3 = 1
bs12 = 1
bs13 = 1
bs23 = 1
ls1 = 2e-3
ls2 = 2e-3
ls3 = 2e-3
ks3 = 114
ks1 = 430
ks2 = 30
dt = 6000
kr = 50
rmrna = 5.7
rlux = 66
rcin = 22
rain = 142.2
ds1 = 0.108
ds2 = 0.017
ds3 = 0.53
ltac = 1.5e-3
lsal = 2.1e-4
lara = 2.1e-4
Ktac = 1.4e+5
Ksal = 4.3e+4
Kara = 7.3e+3
D = 600
Q = 3e+4
F = 0.3 #0.5 for continuous, 0 for batch
Kcat = 79812
KM = 2.28e9
AHLlac = 2.28e8
nutf = 1000
ktox = 5
kc = 0.6
d = 0.01
k1nut = 8
k2nut = 8
k3nut = 8
y1 = 2
y2 = 2
y3 = 2
dc = 0.8
# was missing
mat10 = 0
a10 = 0
at10 = 0
mat20 = 0
a20 = 0 
at20 = 0 
mat30 = 0
a30 = 0
at30 = 0
sal = 15000
i = 500000
ara= 8000

c2growth =1
c3growth = 1

## initial values 3 member mRNA
mt10 = 0
t10 = 0
ma10=0
mc10 = 0
mt20 = 0
t20 = 0
ma20 = 0
mc20 = 0
mt30 = 0
t30 = 0
ma30 = 0
mc30 = 0
s10 = 0
s20 = 0
s30 = 0
cini0 = 0
luxi0 = 0
aini0 = 0
c10 = 500
c20 = 500
c30 = 500
nut0 = 20000

# initial values
y0_two_protein = [mt10, t10, mat10, a10, at10, mt20, t20, mat20, a20, at20, s10, s20, cini0, luxi0, c10, c20, nut0]

y0_two_mRNA = [mt10, t10, ma10, mc10, mt20, t20, ma20, mc20, s10, s20, cini0, luxi0, c10, c20, nut0]

y0_three_protein = [mt10, t10, mat10, a10, at10, mt20, t20, mat20, a20, at20, mt30, t30, mat30, a30, at30, s10, s20, s30, cini0, luxi0, aini0, c10, c20, c30, nut0]

y0_three_mRNA = [mt10, t10, ma10, mc10, mt20, t20, ma20, mc20, mt30, t30, ma30, mc30, s10, s20, s30, cini0, luxi0, aini0, c10, c20, c30, nut0]

# parameters
param_two_protein = [bs1, ls1, ks1, drna, bt, kd, kb, dt, bs2, ls2, ks2, kr, rmrna, lsal, Ksal, ltac, Ktac, ds1, ds2, D, Q, rlux, rcin, F, nutf, ktox, kc, d, k1nut, k2nut, sal, i, ara, y1, y2, dc,AHLlac,KM,Kcat]

param_two_mRNA = [bs1, ls1, ks1, kdrna, karna, drna, bs2, ls2, ks2, karna, dkcat, en, kmrnase, bt, dt, rmrna, lsal, Ksal, ltac, Ktac, ds1, ds2, D, Q, rlux, rcin, F, nutf, ktox, kc, d, k1nut, k2nut, sal, i, ara, y1, y2,AHLlac,KM,Kcat]

param_three_protein = [bs1, ls1, ks1, drna, bt, kd, kb, dt, bs2, ls2, ks2, kr, bs23,ls3,ks3,bs13,bs3,bs12, rmrna, lsal, Ksal, ltac, Ktac, lara, Kara, ds1, ds2, ds3, D, Q, rlux, rcin, rain, F, nutf, ktox, kc, d, k1nut, k2nut, k3nut, sal, i, ara, y1, y2, y3,AHLlac,KM,Kcat]

param_three_mRNA = [bs1, ls1, ks1, kdrna, karna, drna, bs2, ls2, ks2, karna, dkcat, en, kmrnase, bt, dt,bs23,ls3,ks3,bs13, bs3, bs12, rmrna, lsal, Ksal, ltac, Ktac, lara, Kara, ds1, ds2, ds3, D, Q, rlux, rcin, rain, F, nutf, ktox, kc, d, k1nut, k2nut, k3nut, sal, i, ara, y1, y2, y3,AHLlac,KM,Kcat,c2growth,c3growth]
