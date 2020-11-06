import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import odeint
from PBL1_params_final import *

def two_protein(y, t, params): #2 member protein
    mt1, t1, mat1, a1, at1, mt2, t2, mat2, a2, at2, s1, s2, cini, luxi, c1, c2, nut = y
    bs1, ls1, ks1, drna, bt, kd, kb, dt, bs2, ls2, ks2, kr, rmrna, lsal, Ksal, ltac, Ktac, ds1, ds2, D, Q, rlux, rcin, F, nutf, ktox, kc, d, k1nut, k2nut, sal, i, ara, y1, y2, dc,AHLlac,KM,Kcat = params
    derivs = [bs1*(ls1+s1**2/(ks1+s1**2))-drna*mt1, #mt1
        bt*mt1+kd*at1-kb*a1*t1-dt*t1, #t1
        bs2*(ls2+s2**2/(ks2+s2**2))-drna*mat1,# mat1
        kr*bt*mt1+kd*at1-kb*a1*t1-dt*a1,#a1
        kb*a1*t1 - kd*at1 - dt*at1, #at1
        bs2*(ls2+s2**2/(ks2+s2**2))-drna*mt2, #mt2
        bt*mt2+kd*at2-kb*a2*t2-dt*t2, #t2
        bs1*(ls1+s1**2/(ks1+s1**2))-drna*mat2, #mat2
        kr*bt*mt2+kd*at2-kb*a2*t2-dt*a2, #a2
        kb*a2*t2 - kd*at2 - dt*at2, #at2
        rlux*luxi*c1-ds1*s1-Kcat*(c1+c2)*AHLlac*s1/(KM+s1)-F*s1, #s1
        rcin*cini*c2-ds2*s2-Kcat*(c1+c2)*AHLlac*s2/(KM+s2)-F*s2, #s2
        rmrna*(lsal+sal**2/(Ksal+sal**2))-D*cini/(Q+cini), #cini
        rmrna*(ltac+i**2/(Ktac+i**2))-D*luxi/(Q+luxi), #luxi
        kc*(nut/(k1nut+nut))*c1-dc*c1*t1/(ktox+t1)-d*c1-F*c1, #c1
        kc*(nut/(k2nut+nut))*c2-dc*c2*t2/(ktox+t2)-d*c2-F*c2, #c2
        nutf*F-nut*F-kc*(nut/(k1nut+nut))*c1*y1-kc*(nut/(k2nut+nut))*c2*y2 #nutrient
    ]
    return derivs

def two_mRNA(y, t, params): #2 member mRNA
    mt1, t1, ma1, mc1, mt2, t2, ma2, mc2, s1, s2, cini, luxi, c1, c2, nut = y
    bs1, ls1, ks1, kdrna, karna, drna, bs2, ls2, ks2, karna, dkcat, en, kmrnase, bt, dt, rmrna, lsal, Ksal, ltac, Ktac, ds1, ds2, D, Q, rlux, rcin, F, nutf, ktox, kc, d, k1nut, k2nut, sal, i, ara, y1, y2,AHLlac,KM,Kcat = params
    derivs = [bs1*(ls1+s1**2/(ks1+s1**2))+kdrna*mc1-karna*ma1*mt1-drna*mt1,  #mt1
              bt*mt1-dt*t1,  #t1
              bs2 * (ls2+s2**2/(ks2+s2**2)) + kdrna * mc1 - karna * ma1 * mt1 - drna * ma1,  # ma1
              karna * ma1 * mt1 - kdrna * mc1 - dkcat * en * mc1 / (kmrnase+mc1),  #mc1
              bs2 * (ls2+s2**2/(ks2+s2**2)) + kdrna * mc2 - karna * ma2 * mt2 - drna * mt2,  #mt2
              bt * mt2 - dt * t2,  #t2
              bs1 * (ls1+s1**2/(ks1+s1**2)) + kdrna * mc2 - karna * ma2 * mt2 - drna * ma2,  #ma2
              karna * ma2 * mt2 - kdrna * mc2 - dkcat * en * mc2 / (kmrnase+mc2),  #mc2
              rlux * luxi * c1 - ds1 * s1 - Kcat * (c1+c2) * AHLlac * s1 / (KM+s1) - F * s1,  #s1
              rcin * cini * c2 - ds2 * s2 - Kcat * (c1+c2) * AHLlac * s2 / (KM+s2) - F * s2,  #s2
              rmrna * (lsal+sal**2/(Ksal+sal**2)) - D * cini / (Q+cini),  #cini
              rmrna * (ltac+i**2/(Ktac+i**2)) - D * luxi / (Q+luxi), #luxi
              kc * (nut / (k1nut + nut)) * c1 - dc * c1 * t1 / (ktox + t1) - d * c1 - F * c1,  # c1
              kc * (nut / (k2nut + nut)) * c2 - dc * c2 * t2 / (ktox + t2) - d * c2 - F * c2,  # c2
              nutf * F - nut * F - kc * (nut / (k1nut + nut)) * c1 * y1 - kc * (nut / (k2nut + nut)) * c2 * y2  # nutrient
              ]
    return derivs

def three_protein(y, t, params): #3 member protein
    mt1, t1, mat1, a1, at1, mt2, t2, mat2, a2, at2, mt3, t3, mat3, a3,at3, s1, s2, s3, cini, luxi, aini, c1, c2, c3, nut = y
    bs1, ls1, ks1, drna, bt, kd, kb, dt, bs2, ls2, ks2, kr, bs23,ls3,ks3,bs13,bs3,bs12, rmrna, lsal, Ksal, ltac, Ktac, lara, Kara, ds1, ds2, ds3, D, Q, rlux, rcin, rain, F, nutf, ktox, kc, d, k1nut, k2nut, k3nut, sal, i, ara, y1, y2, y3,AHLlac,KM,Kcat = params
    derivs = [bs1*(ls1+s1**2/(ks1+s1**2))-drna*mt1,  #mt1
              bt*mt1+kd*at1-kb*a1*t1-dt*t1,  #t1
              bs23*(ls2*ls3+(s2**2/(ks2+s2**2)*(s3**2/(ks3+s3**2))))-drna*mat1,  # mat1
              kr*bt*mt1+kd*at1-kb*a1*t1-dt*a1,  #a1
              kb*a1*t1 - kd*at1 - dt*at1,  #at1
              bs2*(ls2+s2**2/(ks2+s2**2))-drna*mt2,  #mt2
              bt*mt2+kd*at2-kb*a2*t2-dt*t2,  #t2
              bs13*(ls1*ls3+(s1**2/(ks1+s1**2)*(s3**2/(ks3+s3**2))))-drna*mat2,  #mat2
              kr*bt*mt2+kd*at2-kb*a2*t2-dt*a2,  #a2
              kb*a2*t2 - kd*at2 - dt*at2,  #at2
              bs3*(ls2+s3**2/(ks3+s3**2))-drna*mt3,  #mt3
              bt * mt3 + kd * at3 - kb * a3 * t3 - dt * t3,  # t3
              bs12 * (ls1*ls2+(s1**2/(ks1+s1**2))*(s2**2/(ks2+s2**2))) - drna * mat3,  # mat3
              kr * bt * mt3 + kd * at3 - kb * a3 * t3 - dt * a3,  #a3
              kb * a3 * t3 - kd * at3 - dt * at3,  #at3
              rlux * luxi * c1 - ds1 * s1 - Kcat * (c1+c2+c3) * AHLlac * s1 / (KM+s1) - F * s1,  #s1
              rcin * cini * c2 - ds2 * s2 - Kcat * (c1+c2+c3) * AHLlac * s2 / (KM+s2) - F * s2,  #s2
              rain * aini * c3 - ds3 * s3 - Kcat * (c1+c2+c3) * AHLlac * s3 / (KM+s3) - F * s3,  #s3
              rmrna * (lsal+sal**2/(Ksal+sal**2)) - D * cini / (Q+cini),  #cini
              rmrna * (ltac+i**2/(Ktac+i**2)) - D * luxi / (Q+luxi),  #luxi
              rmrna * (lara+ara**2/(Kara+ara**2)) - D * aini / (Q+aini), #aini
              kc * (nut / (k1nut + nut)) * c1 - dc * c1 * t1 / (ktox + t1) - d * c1 - F * c1,  # c1
              kc * (nut / (k2nut + nut)) * c2 - dc * c2 * t2 / (ktox + t2) - d * c2 - F * c2,  # c2
              kc * (nut / (k2nut + nut)) * c3 - dc * c3 * t3 / (ktox + t3) - d * c3 - F * c3,  # c3
              nutf * F - nut * F - kc * (nut / (k1nut + nut)) * c1 * y1 - kc * (nut / (k2nut + nut)) * c2 * y2 - kc * (nut / (k3nut + nut)) * c3 * y3 # nutrient
              ]
    return derivs

def three_mRNA(y, t, params): #3 member mRNA
    mt1, t1, ma1, mc1, mt2, t2, ma2, mc2, mt3, t3, ma3, mc3, s1, s2, s3, cini, luxi, aini, c1, c2,c3, nut = y
    bs1, ls1, ks1, kdrna, karna, drna, bs2, ls2, ks2, karna, dkcat, en, kmrnase, bt, dt,bs23,ls3,ks3,bs13, bs3, bs12, rmrna, lsal, Ksal, ltac, Ktac, lara, Kara, ds1, ds2, ds3, D, Q, rlux, rcin, rain, F, nutf, ktox, kc, d, k1nut, k2nut, k3nut, sal, i, ara, y1, y2, y3,AHLlac,KM,Kcat,c2growth,c3growth = params
    derivs = [bs1*(ls1+s1**2/(ks1+s1**2))+kdrna*mc1-karna*ma1*mt1-drna*mt1,  #mt1
              bt*mt1-dt*t1,  #t1
              bs23*(ls2*ls3+(s2**2/(ks2+s2**2)*(s3**2/(ks3+s3**2))))+kdrna*mc1-karna*ma1*mt1-drna*ma1,  # ma1
              karna*ma1*mt1-kdrna*mc1-dkcat*en*mc1/(kmrnase+mc1),  #mc1
              bs2*(ls2+s2**2/(ks2+s2**2))+kdrna*mc2-karna*ma2*mt2-drna*mt2,  #mt2
              bt*mt2-dt*t2,  #t2
              bs13*(ls1*ls3+(s1**2/(ks1+s1**2)*(s3**2/(ks3+s3**2))))+kdrna*mc2-karna*ma2*mt2-drna*ma2,  #ma2
              karna*ma2*mt2-kdrna*mc2-dkcat*en*mc2/(kmrnase+mc2),  #mc2
              bs3 * (ls3+s3**2/(ks3+s3**2)) + kdrna * mc3 - karna * ma3 * mt3 - drna * mt3,  # mt3
              bt * mt3 - dt * t3,  #t3
              bs12 * (ls1*ls2+(s1**2/(ks1+s1**2))*(s2**2/(ks2+s2**2))) + kdrna * mc3 - karna * ma3 * mt3 - drna * ma3,  #ma3
              karna * ma3 * mt3 - kdrna * mc3 - dkcat * en * mc3 / (kmrnase+mc3),  #mc3
              rlux * luxi * c1 - ds1 * s1 - Kcat * (c1+c2+c3) * AHLlac * s1 / (KM+s1) - F * s1,  #s1
              rcin * cini * c2 - ds2 * s2 - Kcat * (c1+c2+c3) * AHLlac * s2 / (KM+s2) - F * s2,  #s2
              rain * aini * c3 - ds3 * s3 - Kcat * (c1+c2+c3) * AHLlac * s3 / (KM+s3) - F * s3,  #s3
              rmrna * (lsal+sal**2/(Ksal+sal**2)) - D * cini / (Q+cini),  #cini
              rmrna * (ltac+i**2/(Ktac+i**2)) - D * luxi / (Q+luxi),  #luxi
              rmrna * (lara+ara**2/(Kara+ara**2)) - D * aini / (Q+aini), #aini
              kc * (nut / (k1nut + nut)) * c1 - dc * c1 * t1 / (ktox + t1) - d * c1 - F * c1,  # c1
              c2growth* kc * (nut / (k2nut + nut)) * c2 - dc * c2 * t2 / (ktox + t2) - d * c2 - F * c2,  # c2
              c3growth* kc * (nut / (k2nut + nut)) * c3 - dc * c3 * t3 / (ktox + t3) - d * c3 - F * c3,  # c3
              nutf * F - nut * F - kc * (nut / (k1nut + nut)) * c1 * y1 - kc * (nut / (k2nut + nut)) * c2 * y2 - kc * (nut / (k3nut + nut)) * c3 * y3# nutrient
              ]
    return derivs



def run_simulation(func, params, y0,tstop):

    tInc = 0.5
    t = np.arange(0.,tstop,tInc)
    y0 = np.array(y0)

    psoln, infodict = odeint(func, y0, t, args=(params,),full_output=True)

    return t, psoln

def plot_two(t, psoln, i_1, i_2):
    plt.plot(t,psoln[:,i_1])
    plt.plot(t, psoln[:, i_2])
    plt.legend(['c1', 'c2'])
    plt.tight_layout()
    plt.show()

def plot_three(t, psoln, i_1, i_2, i_3):
    plt.plot(t, psoln[:, i_1])
    plt.plot(t, psoln[:, i_2])
    plt.plot(t,psoln[:,i_3])
    plt.legend(['c1', 'c2', 'c3'])
    plt.tight_layout()
    plt.show()


def different_initial(init_cond, change, param_number_2mrna, param_number_3mrna, perturbation):
    # init_conditions - a list of initial conditions to test. If more than one changing, making a list of lists - [[1,2],[3,4]], no more than 3 can be tested using this code
    # change = enter 1 for parameters, 0 for initial conditions
    # param_number - the number of the changed parameters in the list of init conditions or parameters, as a list in the same order as initial conditions
    # perturbation - 1 if add 100 cells, 0 if no perturbation
    fig, axs = plt.subplots(2,4,sharex='col', sharey='row',gridspec_kw={'hspace': 0, 'wspace': 0})
    for F in [0, 0.3]:  # for both batch and chemostat
        if F == 0:
            location = 0
            culture = 'Batch Culture'
            tstop = 1000
        else:
            location = 1
            culture = 'Chemostat'
            tstop = 20000
        #2 member
        params = param_two_mRNA.copy()
        params[26] = F
        y0 = y0_two_mRNA.copy()
        for ii in range(len(init_cond)):
            if ii == 0:
                lines = ':'
            elif ii == 1:
                lines = '--'
            elif ii == 2:
                lines = '-'
            elif ii == 3:
                lines = '-.'
            if change == 1:
                for jj in range(len(param_number_2mrna)):
                    params[param_number_2mrna[jj]] = init_cond[ii][jj]
            else:
                for jj in range(len(param_number_2mrna)):
                    y0[param_number_2mrna[jj]] = init_cond[ii][jj]
            t, psoln = run_simulation(two_mRNA, params, y0,tstop)
            t = np.array(t)
            tpoints = len(t)
            yinter = y0_two_mRNA.copy()
            c1 = np.array(psoln[:, 12])
            c2 = np.array(psoln[:, 13])
            if perturbation == 1:
                for jj in range(len(y0_two_mRNA)):
                    yinter[jj] = psoln[tpoints - 1, jj]
                yinter[12] = yinter[12] + 100
                t2, psoln2 = run_simulation(two_mRNA, params, yinter,tstop)
                t2 = np.array(t2)
                t2 = t2 + t[tpoints - 1]
                t = np.concatenate((t, t2))
                c12 = np.array(psoln2[:, 12])
                c22 = np.array(psoln2[:, 13])
                c1 = np.concatenate((c1, c12))
                c2 = np.concatenate((c2, c22))
            ctotal = c1 + c2
            for kk in range(len(c1)):
                if c1[kk] < 1:
                    c1[kk] = 0
            for kk in range(len(c2)):
                if c2[kk] < 1:
                    c2[kk] = 0
            c1frac = np.divide(c1, ctotal, out=np.zeros_like(c1), where=ctotal != 0)
            c2frac = np.divide(c2, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
            axs[1, location].plot(t, c1frac, color='red', linestyle=lines)
            axs[1, location].plot(t, c2frac, color='blue', linestyle=lines)
            axs[0, location].plot(t, c1, color='red', linestyle=lines)
            axs[0, location].plot(t, c2, color='blue', linestyle=lines)
        # 3 member
        params = param_three_mRNA.copy()
        params[36] = F
        y0 = y0_three_mRNA.copy()
        for ii in range(len(init_cond)):
            if ii == 0:
                lines = ':'
            elif ii == 1:
                lines = '--'
            elif ii == 2:
                lines = '-'
            elif ii == 3:
                lines='-.'
            if change == 1:
                for jj in range(len(param_number_3mrna)):
                    params[param_number_3mrna[jj]] = init_cond[ii][jj]
            else:
                for jj in range(len(param_number_3mrna)):
                    y0[param_number_3mrna[jj]] = init_cond[ii][jj]
            t, psoln = run_simulation(three_mRNA, params, y0,tstop)
            t = np.array(t)
            tpoints = len(t)
            yinter = y0_three_mRNA.copy()
            c1 = np.array(psoln[:, 18])
            c2 = np.array(psoln[:, 19])
            c3 = np.array(psoln[:, 20])
            if perturbation == 1:
                for jj in range(len(y0_three_mRNA)):
                    yinter[jj] = psoln[tpoints - 1, jj]
                yinter[18] = yinter[18] + 100
                t2, psoln2 = run_simulation(three_mRNA, params, yinter,tstop)
                t2 = np.array(t2)
                t2 = t2 + t[tpoints - 1]
                t = np.concatenate((t, t2))
                c12 = np.array(psoln2[:, 18])
                c22 = np.array(psoln2[:, 19])
                c32 = np.array(psoln2[:, 20])
                c1 = np.concatenate((c1, c12))
                c2 = np.concatenate((c2, c22))
                c3 = np.concatenate((c3, c32))
            ctotal = c1 + c2 + c3
            for kk in range(len(c1)):
                if c1[kk] < 1:
                    c1[kk] = 0
            for kk in range(len(c2)):
                if c2[kk] < 1:
                    c2[kk] = 0
            for kk in range(len(c3)):
                if c3[kk] < 1:
                    c3[kk] = 0
            c1frac = np.divide(c1, ctotal, out=np.zeros_like(c1), where=ctotal != 0)
            c2frac = np.divide(c2, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
            c3frac = np.divide(c3, ctotal, out=np.zeros_like(c3), where=ctotal != 0)
            axs[1, location+2].plot(t, c1frac, color='red', linestyle=lines)
            axs[1, location+2].plot(t, c2frac, color='blue', linestyle=lines)
            axs[1, location+2].plot(t, c3frac, color='green', linestyle=lines)
            axs[0, location+2].plot(t, c1, color='red', linestyle=lines)
            axs[0, location+2].plot(t, c2, color='blue', linestyle=lines)
            axs[0, location+2].plot(t, c3, color='green', linestyle=lines)
    fig.suptitle('Population Dynamics with Different Initial Ratios')
    custom_lines = [Line2D([0], [0], color='red', linestyle='-'), Line2D([0], [0], color='blue', linestyle='-'),
                    Line2D([0], [0], color='green', linestyle='-'), Line2D([0], [0], color='black', linestyle=':'),
                    Line2D([0], [0], color='black', linestyle='--'), Line2D([0], [0], color='black', linestyle='-')]
    axs[0,3].legend(custom_lines, ['C1', 'C2', 'C3', 'C1:100, C2:100, C3:100','C1: 50, C2:75, C3: 175','C1: 210, C2: 75,C3: 15'])
    for ax in axs.flat:
        ax.set(xlabel='Time(Hours)')
        ax.label_outer()
    axs[0, 0].set(ylabel='Cell Population')
    axs[1, 0].set(ylim=(0, 1.1))
    axs[0, 0].set(ylim=(0, 1250))
    axs[1, 0].set(ylabel='Cell Proportion')
    axs[1, 1].set(ylim=(0, 1.1))
    axs[0, 1].set(ylim=(0, 2500))
    axs[0, 0].set_title('Two Members \nin Batch Culture')
    axs[0, 1].set_title('Two Members \nin Chemostat')
    axs[0, 2].set_title('Three Members\n in Batch Culture')
    axs[0, 3].set_title('Three Members \nin Chemostat')
    plt.show()

def timecourse_both_1condition(perturbation):
    #perturbation = 1 if perturbing system, 0 if not
    fig, axs = plt.subplots(2, 4, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
    for mode in ['protein', 'mrna']:  # to do both protein and mRNA circuits
        for F in [0, 0.3]:  # for both batch and chemostat
            if F == 0:
                location = 0
                tstop = 500
            else:
                location = 1
                tstop = 15000
            if mode == 'protein':
                params = param_two_protein.copy()
                params[23] = F
                y0 = y0_two_protein.copy()
                t, psoln = run_simulation(two_protein, params, y0,tstop)
                t = np.array(t)
                tpoints = len(t)
                yinter = y0_two_protein.copy()
                c1 = np.array(psoln[:, 14])
                c2 = np.array(psoln[:, 15])
                if perturbation == 1:
                    for jj in range(len(y0_two_protein)):
                        yinter[jj] = psoln[tpoints - 1, jj]
                    yinter[14] = yinter[14] + 100
                    t2, psoln2 = run_simulation(two_protein, params, yinter,tstop)
                    t2 = np.array(t2)
                    t2 = t2 + t[tpoints - 1]
                    t = np.concatenate((t, t2))
                    c12 = np.array(psoln2[:, 14])
                    c22 = np.array(psoln2[:, 15])
                    c1 = np.concatenate((c1, c12))
                    c2 = np.concatenate((c2, c22))
                for kk in range(len(c1)):
                    if c1[kk]<1:
                        c1[kk] = 0
                for kk in range(len(c2)):
                    if c2[kk] < 1:
                        c2[kk] = 0
                ctotal = c1 + c2
                c1frac = np.divide(c1, ctotal, out=np.zeros_like(c1), where=ctotal != 0)
                c2frac = np.divide(c2, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
                axs[1, location].plot(t, c1frac, color='red', linestyle=':')
                axs[1, location].plot(t, c2frac, color='blue', linestyle=':')
                axs[0, location].plot(t, c1, color='red', linestyle=':')
                axs[0, location].plot(t, c2, color='blue', linestyle=':')
                # onto 3 member
                params = param_three_protein.copy()
                params[33] = F
                y0 = y0_three_protein.copy()
                t, psoln = run_simulation(three_protein, params, y0,tstop)
                t = np.array(t)
                tpoints = len(t)
                yinter = y0_three_protein.copy()
                c1 = np.array(psoln[:, 21])
                c2 = np.array(psoln[:, 22])
                c3 = np.array(psoln[:, 23])
                if perturbation == 1:
                    for jj in range(len(y0_three_protein)):
                        yinter[jj] = psoln[tpoints - 1, jj]
                    yinter[21] = yinter[21] + 100
                    t2, psoln2 = run_simulation(three_protein, params, yinter,tstop)
                    t2 = np.array(t2)
                    t2 = t2 + t[tpoints - 1]
                    t = np.concatenate((t, t2))
                    c12 = np.array(psoln2[:, 21])
                    c22 = np.array(psoln2[:, 22])
                    c32 = np.array(psoln2[:, 23])
                    c1 = np.concatenate((c1, c12))
                    c2 = np.concatenate((c2, c22))
                    c3 = np.concatenate((c3, c32))
                for kk in range(len(c1)):
                    if c1[kk]<1:
                        c1[kk] = 0
                for kk in range(len(c2)):
                    if c2[kk] < 1:
                        c2[kk] = 0
                for kk in range(len(c3)):
                    if c3[kk] < 1:
                        c3[kk] = 0
                ctotal = c1 + c2 + c3
                c1frac = np.divide(c1, ctotal, out=np.zeros_like(c1), where=ctotal != 0)
                c2frac = np.divide(c2, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
                c3frac = np.divide(c3, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
                axs[1, location].plot(t, c1frac, color='red', linestyle='-')
                axs[1, location].plot(t, c2frac, color='blue', linestyle='-')
                axs[1, location].plot(t, c3frac, color='green', linestyle='-')
                axs[0, location].plot(t, c1, color='red', linestyle='-')
                axs[0, location].plot(t, c2, color='blue', linestyle='-')
                axs[0, location].plot(t, c3, color='green', linestyle='-')
            if mode == 'mrna':
                params = param_two_mRNA.copy()
                params[26] = F
                y0 = y0_two_mRNA.copy()
                t, psoln = run_simulation(two_mRNA, params, y0,tstop)
                t = np.array(t)
                tpoints = len(t)
                yinter = y0_two_mRNA.copy()
                c1 = np.array(psoln[:, 12])
                c2 = np.array(psoln[:, 13])
                if perturbation == 1:
                    for jj in range(len(y0_two_mRNA)):
                        yinter[jj] = psoln[tpoints - 1, jj]
                    yinter[12] = yinter[12] + 100
                    t2, psoln2 = run_simulation(two_mRNA, params, yinter,tstop)
                    t2 = np.array(t2)
                    t2 = t2 + t[tpoints - 1]
                    t = np.concatenate((t, t2))
                    c12 = np.array(psoln2[:, 12])
                    c22 = np.array(psoln2[:, 13])
                    c1 = np.concatenate((c1, c12))
                    c2 = np.concatenate((c2, c22))
                for kk in range(len(c1)):
                    if c1[kk]<1:
                        c1[kk] = 0
                for kk in range(len(c2)):
                    if c2[kk] < 1:
                        c2[kk] = 0
                ctotal = c1 + c2
                c1frac = np.divide(c1, ctotal, out=np.zeros_like(c1), where=ctotal != 0)
                c2frac = np.divide(c2, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
                axs[1, location+2].plot(t, c1frac, color='red', linestyle=':')
                axs[1, location+2].plot(t, c2frac, color='blue', linestyle=':')
                axs[0, location+2].plot(t, c1, color='red', linestyle=':')
                axs[0, location+2].plot(t, c2, color='blue', linestyle=':')
                # start 3 member
                params = param_three_mRNA.copy()
                params[36] = F
                y0 = y0_three_mRNA.copy()
                t, psoln = run_simulation(three_mRNA, params, y0,tstop)
                t = np.array(t)
                tpoints = len(t)
                yinter = y0_three_mRNA.copy()
                c1 = np.array(psoln[:, 18])
                c2 = np.array(psoln[:, 19])
                c3 = np.array(psoln[:, 20])
                if perturbation == 1:
                    for jj in range(len(y0_three_mRNA)):
                        yinter[jj] = psoln[tpoints - 1, jj]
                    yinter[18] = yinter[18] + 100
                    t2, psoln2 = run_simulation(three_mRNA, params, yinter,tstop)
                    t2 = np.array(t2)
                    t2 = t2 + t[tpoints - 1]
                    t = np.concatenate((t, t2))
                    c12 = np.array(psoln2[:, 18])
                    c22 = np.array(psoln2[:, 19])
                    c32 = np.array(psoln2[:, 20])
                    c1 = np.concatenate((c1, c12))
                    c2 = np.concatenate((c2, c22))
                    c3 = np.concatenate((c3, c32))
                for kk in range(len(c1)):
                    if c1[kk]<1:
                        c1[kk] = 0
                for kk in range(len(c2)):
                    if c2[kk] < 1:
                        c2[kk] = 0
                for kk in range(len(c3)):
                    if c3[kk] < 1:
                        c3[kk] = 0
                ctotal = c1 + c2 + c3
                c1frac = np.divide(c1, ctotal, out=np.zeros_like(c1), where=ctotal != 0)
                c2frac = np.divide(c2, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
                c3frac = np.divide(c3, ctotal, out=np.zeros_like(c2), where=ctotal != 0)
                axs[1, location+2].plot(t, c1frac, color='red', linestyle='-')
                axs[1, location+2].plot(t, c2frac, color='blue', linestyle='-')
                axs[1, location+2].plot(t, c3frac, color='green', linestyle='-')
                axs[0, location+2].plot(t, c1, color='red', linestyle='-')
                axs[0, location+2].plot(t, c2, color='blue', linestyle='-')
                axs[0, location+2].plot(t, c3, color='green', linestyle='-')
    fig.suptitle('System Response to Perturbation in Batch and Continuous Culture')
    custom_lines = [Line2D([0], [0], color='red', linestyle='-'), Line2D([0], [0], color='blue', linestyle='-'),
                    Line2D([0], [0], color='green', linestyle='-'),
                        Line2D([0], [0], color='black', linestyle=':'),
                        Line2D([0], [0], color='black', linestyle='-')]
    axs[0,3].legend(custom_lines, ['C1', 'C2','C3', '2 Member System', '3 Member System'])
    for ax in axs.flat:
        ax.set(xlabel='Time(Hours)')
        axs[0, 0].set(ylabel='Cell Population')
        axs[1,0].set(ylim=(0,1.1))
        axs[0, 0].set(ylim=(0, 1250))
        axs[1, 0].set(ylabel='Cell Proportion')
        axs[1, 1].set(ylim=(0, 1.1))
        ax.label_outer()
        axs[0, 1].set(ylim=(0, 2000))
        axs[0, 0].set_title('Protein Circuit in Batch Culture')
        axs[0, 1].set_title('Protein Circuit in Chemostat')
        axs[0, 2].set_title('mRNA Circuit in Batch Culture')
        axs[0, 3].set_title('mRNA Circuit in Chemostat')
    plt.show()

def two_parameter_scan_2m(parameter1_range,parameter1_listnum,parameter2_range,parameter2_listnum,membernumber,show,divisor):
    #parameter range - conditions to be tested in format of (min,max,inc)
    # paramater_listnum - number it appears on the mrna list, whether it be 2 or 3 member
    # membernumber - number of members - either 2 or 3
    # show = 1, 2 or 3 for cell type you want to show
    x1 = []
    x2=[]
    y = []
    if membernumber == 2:
        params = param_two_mRNA.copy()
        y0 = y0_two_mRNA.copy()
    else:
        params = param_three_mRNA.copy()
        y0 = y0_three_mRNA.copy()
    for ii in range(parameter1_range[0],parameter1_range[1],parameter1_range[2]):
        params[parameter1_listnum] = ii/divisor
        for jj in range(parameter2_range[0],parameter2_range[1],parameter2_range[2]):
            params[parameter2_listnum] = jj/divisor
            x1.append(ii/divisor)
            x2.append(jj/divisor)
            if membernumber==2:
                t1, psoln = run_simulation(two_mRNA, params, y0,20000)
            else:
                t1, psoln = run_simulation(three_mRNA, params, y0,20000)
            t1 = np.array(t1)
            tpoints = len(t1)
            if membernumber == 2:
                c1 = np.array(psoln[:, 12])
                c2 = np.array(psoln[:, 13])
                if c1[tpoints-1] <0:
                    c1[tpoints - 1]
                if c2[tpoints-1] <0:
                    c2[tpoints - 1]
                ctotal = c1[tpoints-1] + c2 [tpoints-1]
                if show == 1:
                    if ctotal != 0:
                        y.append(c1[tpoints-1]/ctotal)
                    else:
                        y.append(0)
                else:
                    if ctotal != 0:
                        y.append(c2[tpoints-1]/ctotal)
                    else:
                        y.append(0)
            else:
                c1 = np.array(psoln[:, 18])
                c2 = np.array(psoln[:, 19])
                c3 = np.array(psoln[:, 20])
                if c1[tpoints-1] <0:
                    c1[tpoints - 1]
                if c2[tpoints-1] <0:
                    c2[tpoints - 1]
                if c2[tpoints-1] <0:
                    c2[tpoints - 1]
                if c3[tpoints-1]<0:
                    c3[tpoints-1] = 0
                ctotal = c1[tpoints-1] + c2[tpoints-1] + c3[tpoints-1]
                if show == 1:
                    if ctotal != 0:
                        y.append(c1[tpoints-1]/ctotal)
                    else:
                        y.append(0)
                elif show ==2:
                    if ctotal != 0:
                        y.append(c2[tpoints-1]/ctotal)
                    else:
                        y.append(0)
                else:
                    if ctotal != 0:
                        y.append(c3[tpoints-1]/ctotal)
                    else:
                        y.append(0)
    return x1,x2,y

def one_parameter_scan(parameter_range,param_num,members,divisor):
    x = []
    c1ss = []
    c2ss = []
    c3ss = []
    for ii in range(parameter_range[0],parameter_range[1],parameter_range[2]):
        x.append(ii/divisor)
        if members == 2:
            params = param_two_mRNA.copy()
            params[param_num] = ii/divisor
            t1, psoln = run_simulation(two_mRNA, params, y0_two_mRNA,20000)
            t1 = np.array(t1)
            c1 = np.array(psoln[:, 12])
            c2 = np.array(psoln[:, 13])
            tpoints = len(t1)
            if c1[tpoints - 1] < 0:
                c1[tpoints - 1]
            if c2[tpoints - 1] < 0:
                c2[tpoints - 1]
            if c2[tpoints - 1] < 0:
                c2[tpoints - 1]
            if c3[tpoints - 1] < 0:
                c3[tpoints - 1] = 0
            ctotal = c1[tpoints - 1] + c2[tpoints - 1]
            c1ss.append(c1[tpoints - 1] / ctotal)
            c2ss.append(c2[tpoints - 1] / ctotal)
            c3ss.append(0)
        else:
            params = param_three_mRNA.copy()
            params[param_num] = ii/divisor
            t1, psoln = run_simulation(three_mRNA, params, y0_three_mRNA,20000)
            t1 = np.array(t1)
            c1 = np.array(psoln[:, 18])
            c2 = np.array(psoln[:, 19])
            c3 = np.array(psoln[:, 20])
            tpoints = len(t1)
            if c1[tpoints - 1] < 0:
                c1[tpoints - 1]
            if c2[tpoints - 1] < 0:
                c2[tpoints - 1]
            if c2[tpoints - 1] < 0:
                c2[tpoints - 1]
            if c3[tpoints - 1] < 0:
                c3[tpoints - 1] = 0
            ctotal = c1[tpoints - 1] + c2[tpoints - 1] + c3[tpoints - 1]
            c1ss.append(c1[tpoints - 1]/ctotal)
            c2ss.append(c2[tpoints - 1]/ctotal)
            c3ss.append(c3[tpoints - 1]/ctotal)
    return x,c1ss,c2ss,c3ss



#timecourse_both_1condition(1)

# t, psoln = run_simulation(two_protein, param_two_protein, y0_two_protein)
# plot_two(t, psoln, 14, 15)

# t, psoln = run_simulation(two_mRNA, param_two_mRNA, y0_two_mRNA)
# plot_two(t, psoln, 12, 13)

# t, psoln = run_simulation(three_protein, param_three_protein, y0_three_protein)
# plot_three(t, psoln, 21, 22, 23)

#t, psoln = run_simulation(three_mRNA,param_three_mRNA,y0_three_mRNA)
#plot_three(t, psoln,18,19,20)

#time_diff_init_3m([[1.08e7],[2.28e8],[3.48e9]],1,[47],[50],1)

different_initial([[100,100,100],[50,75,175],[210,75,15]],0,[12,13],[18,19,20],0)

#different_initial([[0],[2.28e6],[2.28e12]],1,[38],[50],1)

# the following is the code to generate the 3x3 parameter scan figure
fig, axs = plt.subplots(3,3)
#top left
iptg,sal,c1 = two_parameter_scan_2m([0,1100,100],34,[0,1100,100],33,2,1,1)
data1= axs[0,0].scatter(x=iptg,y=sal,c=c1,vmin=0.2,vmax=0.8)
fig.colorbar(data1,ax=axs[0,0])
axs[0,0].set_title('C1 in 2 Member System')
axs[0,0].set(xlabel='IPTG (nM)', ylabel='Sal (nM)')

#top center
iptg,ara,c1=two_parameter_scan_2m([0,2100,200],45,[0,210,20],46,3,1,1)
#sal held constant 500 nM
data2 = axs[0,1].scatter(x=iptg,y=ara,c=c1,vmin=0.2,vmax=0.7)
plt.colorbar(data2,ax=axs[0,1])
axs[0,1].set_title('C1 in 3 Member System')
axs[0,1].set(xlabel='IPTG (nM)', ylabel='Ara (nM)')

#top right
iptg,ara,c3=two_parameter_scan_2m([0,2100,200],45,[0,210,20],46,3,3,1)
#sal held constant 500 nM
data3=axs[0,2].scatter(x=iptg,y=ara,c=c3,vmin=0.05,vmax=0.55)
plt.colorbar(data3,ax=axs[0,2])
axs[0,2].set_title('C3 in 3 Member System')
axs[0,2].set(xlabel='IPTG (nM)',ylabel='Ara (nM)')

#center left
ahllac,c1,c2,c3 = one_parameter_scan([0,10001000000,1000000000],50,3,1)
axs[1,0].plot(ahllac,c1,color='red')
axs[1,0].plot(ahllac,c2,color='blue')
axs[1,0].plot(ahllac,c3,color='green')
initial = []
for ii in range(len(ahllac)):
    initial.append(1/3)
axs[1,0].plot(ahllac,initial,color='black',linestyle=':')
axs[1,0].legend(['C1','C2','C3','Initial Conditions'],loc='upper right')
axs[1,0].set(xlabel='AHL-Lactonase (nM)',ylabel='Proportion of cells',ylim=(0,1))

#center
rnase,c1,c2,c3 = one_parameter_scan([0,20010,1000],11,3,1)
axs[1,1].plot(rnase,c1,color='red')
axs[1,1].plot(rnase,c2,color='blue')
axs[1,1].plot(rnase,c3,color='green')
initial = []
for ii in range(len(rnase)):
    initial.append(1/3)
axs[1,1].plot(rnase,initial,color='black',linestyle=':')
axs[1,1].legend(['C1','C2','C3','Initial Conditions'])
axs[1,1].set(xlabel='RNAse III (nM)',ylabel='Proportion of cells',ylim=(0,1))

#center right
flow,c1,c2,c3 = one_parameter_scan([0,61,2],36,3,100)
# after 0.6, flow rate larger than growth rate so no steady state achieved
axs[1,2].plot(flow,c1,color='red')
axs[1,2].plot(flow,c2,color='blue')
axs[1,2].plot(flow,c3,color='green')
initial = []
for ii in range(len(flow)):
    initial.append(1/3)
axs[1,2].plot(flow,initial,color='black',linestyle=':')
axs[1,2].legend(['C1','C2','C3','Initial Conditions'])
axs[1,2].set(xlabel='Chemostat Flow (mL/hr)',ylabel='Proportion of cells',ylim=(0,1.1))

#bottom left
y2,y3,c1=two_parameter_scan_2m([95,106,1],53,[95,106,1],54,3,1,100)
data4=axs[2,0].scatter(x=y2,y=y3,c=c1,vmin=0,vmax=1)
plt.colorbar(data4,ax=axs[2,0])
axs[2,0].set_title('C1 in 3 Member System')
axs[2,0].set(xlabel='C2/C1 Growth Constant Ratio',ylabel='C3/C1 Growth Constant Ratio')

#bottom center
y2,y3,c3=two_parameter_scan_2m([95,106,1],53,[95,106,1],54,3,3,100)
data5=axs[2,1].scatter(x=y2,y=y3,c=c3,vmin=0,vmax=1)
fig.colorbar(data5, ax=axs[2,1])
axs[2,1].set_title('C3 in 3 Member System')
axs[2,1].set(xlabel='C2/C1 Growth Constant Ratio',ylabel='C3/C1 Growth Constant Ratio')

#bottom right
trans,c1,c2,c3 = one_parameter_scan([100,1000010,10000],13,3,1)
axs[2,2].plot(trans,c1,color='red')
axs[2,2].plot(trans,c2,color='blue')
axs[2,2].plot(trans,c3,color='green')
initial = []
for ii in range(len(trans)):
    initial.append(1/3)
axs[2,2].plot(trans,initial,color='black',linestyle=':')
axs[2,2].legend(['C1','C2','C3','Initial Conditions'],loc='upper right')
axs[2,2].set(xlabel='Toxin mRNA Translation Rate (nM/hr)',ylabel='Proportion of cells',ylim=(0,1.1),xscale='log')

fig.suptitle('Parameter Effects on Cell Proportions')
fig.subplots_adjust(wspace=0.5,hspace=0.7)
plt.show()