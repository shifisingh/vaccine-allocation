import math

import pandas as pd

# import functions and variables from other files
from alpha import findAlpha
from variables import *
from vaccinations import *
from death_rates import du, dh, dq

# initial data for each risk class
# S, E, I, UD, UR, HD, HR, QD, QR, R, D
class1 = [1274292 / population, 1274292 / population, 73164 / population,
          (60970/population) * du[0], (60970/population) * (1-du[0]),
          (14/population) * dh[0], (14/population) * (1-dh[0]),
          (12180/population) * (dq[0]), (12180/population) * (1-dq[0]),
          (95719/population), (11/population)]
class2 = [676765/population, 676765/population, 87102/population,
          72585/population * (du[1]), 72585/population * (1-du[1]),
          14/population * (dh[1]), 14/population * (1-dh[1]),
          14503/population * (dq[1]), 14503/population * (1-dq[1]),
          113953/population, 28/population]
class3 = [1278927/population, 1278927/population, 134136/population,
          111780/population * (du[2]),  111780/population * (1-du[2]),
          28/population * (dh[2]), 28/population * (1-dh[2]),
          22328/population * (dq[2]), 22328/population * (1-dq[2]),
          175487/population, 440/population]
class4 = [1440307/population, 1440307/population, 105396/population,
          87830/population * (du[3]), 87830/population * (1-du[3]),
          140/population * (dh[3]), 140/population * (1-dh[3]),
          17426/population * (dq[3]), 17426/population * (1-dq[3]),
          137886/population, 3066/population]
class5 = [632134/population, 632134/population, 35712/population,
          29760/population * (du[4]), 29760/population * (1-du[4]),
          234/population * (dh[4]), 234/population * (1-dh[4]),
          5718/population * (dq[4]), 5718/population * (1 - dq[4]),
          46721/population, 7896/population]

print('du ', du)
print('dh ', dh)
print('dq ', dq)

#find iud, ihd, iqd rates
#disagreggated rates in appendix page 26

# 2d array that stores all risk class' day 0 data
data = [class1, class2, class3, class4, class5]

def tPlusOne(k, t):
    #global s, e, i, ud, ur, hd, hr, qd, qr, d, r, m, v
    alpha = findAlpha(k, 0, math.ceil(t/7))  # nominal infection rate

    # compute ksum
    ksum = i[0][t-1] + i[1][t-1] + i[2][t-1] + i[3][t-1] + i[4][t-1]

    # Susceptible people, equation 16
    Skt = s[k][t-1] - v[k][t-1]
    s[k][t] = Skt - beta*v[k][t] - ((alpha[0]*gamma)*(Skt - beta*v[k][t])*ksum)

    # Exposed people, equation 17
    Ekt = e[k][t-1]
    e[k][t] = Ekt + ((alpha[0]*gamma*(Skt - beta*v[k][t]))*ksum - (ri*Ekt))

    # six rates
    rDetected = .15
    UDkt = ud[k][t - 1] * (1 - rr)
    rud = rDetected*mMin*(1-pDetected)
    URkt = ud[k][t - 1] * rr
    rur = rDetected*(1-mMin)*(1-pDetected)
    HDkt = hd[k][t - 1] * (1 - rh)  # set of people that are currently hospitalized + going to die
    rhd = rDetected*mMin*pDetected*pHospitalized  # rate that infected people get to a state where they are hospitalized + die
    HRkt = hr[k][t - 1] * rh
    rhr = rDetected*(1-mMin)*pDetected*pHospitalized
    QDkt = qd[k][t - 1] * (1 - rr)
    rqd = rDetected*mMin*pDetected*(1-pHospitalized)
    QRkt = qr[k][t - 1] * rr
    rqr = rDetected*(1-mMin)*pDetected*(1-pHospitalized)

    # people no longer in the infected category for any reason
    rDetected = rud + rur + rhd + rhr + rqd + rqr

    # Infected people, equation 18
    Ikt = i[k][t-1]
    i[k][t] = Ikt + (ri*Ekt - rDetected*Ikt)

    # Undiagnosed people that die, equation 19
    ud[k][t] = ud[k][t-1] + (rud*Ikt - rDeath*ud[k][t-1])
    # Undiagnosed people that recover
    ur[k][t] = URkt + (rur*Ikt - rr*URkt)

    # Hospitalized people that die, equation 20
    # the people hospitalized yesterday + (the people that are infected - the people that die)
    hd[k][t] = hd[k][t-1] + (rhd*Ikt - rDeath*hd[k][t-1])
    # Hospitalized people that recover
    hr[k][t] = HRkt + (rhr*Ikt - rr*HRkt)

    # Quarantined people that die, equation 21
    qd[k][t] = qd[k][t-1] + (rqd*Ikt - rDeath*qd[k][t-1])
    # Quarantined people that recover
    qr[k][t] = QRkt + (rqr*Ikt - rr*QRkt)

    # People that die, equation 22
    Dkt = d[k][t-1]
    d[k][t] = Dkt+rDeath*(UDkt + HDkt + QDkt)

    # People that recover, equation 13
    r[k][t] = rr*(ur[k][t]+qr[k][t])+rh*hr[k][t]

    # Vaccinated people
    Vkt = v[k][t - 1]
    vaxRate = s[k][t]*.01
    v[k][t] = Vkt + vaxRate

    # Mkt, equation 23
    Mkt = m[k][t-1]
    m[k][t] = Mkt + beta*v[k][t]

def populateUntilT(T):
    global data, s, e, i, ud, ur, hd, hr, qd, qr, d, r, m, v

    # populate the preliminary data
    # note: initializing instead of setting to size variable b/c all 2d arrays were set to each other
    # for every risk class
    # [[t0, t1, ...],[t0, t1, ...],[t0, t1, ...],[t0, t1, ...],[t0, t1, ...]]
    s = [[0 for i in range(T + 1)] for i in range(5)]
    e = [[0 for i in range(T + 1)] for i in range(5)]
    i = [[0 for i in range(T + 1)] for i in range(5)]
    ud = [[0 for i in range(T + 1)] for i in range(5)]
    ur = [[0 for i in range(T + 1)] for i in range(5)]
    hd = [[0 for i in range(T + 1)] for i in range(5)]
    hr = [[0 for i in range(T + 1)] for i in range(5)]
    qd = [[0 for i in range(T + 1)] for i in range(5)]
    qr = [[0 for i in range(T + 1)] for i in range(5)]
    r = [[0 for i in range(T + 1)] for i in range(5)]
    d = [[0 for i in range(T + 1)] for i in range(5)]
    m = [[0 for i in range(T + 1)] for i in range(5)]
    v = [[0 for i in range(T + 1)] for i in range(5)]

    #v = buildV(strat, T)
    #v = buildPredictedV()

    for k in range(0, 5):
        for t in range(0, T+1):
            if t == 0:
                s[k][0] = data[k][0]
                e[k][0] = data[k][1]
                i[k][0] = data[k][2]
                ud[k][0] = data[k][3]
                ur[k][0] = data[k][4]
                hd[k][0] = data[k][5]
                hr[k][0] = data[k][6]
                qd[k][0] = data[k][7]
                qr[k][0] = data[k][8]
                r[k][0] = data[k][9]
                d[k][0] = data[k][10]
                m[k][0] = 0
                v[k][0]
            else:
                tPlusOne(k, t)

def buildDF(k, T):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, r, m, v

    populateUntilT(T)

    # set the column names
    # S, E, I, UD, UR, HD, HR, QD, QR, D, R, M, sum
    columns = ['time since day 0', 'susceptible','exposed', 'infected',
               'UD', 'UR', 'HD', 'HR', 'QD', 'QR', 'died', 'recovered', 'immune state', 'vaccinated', 'sum']
    df = pd.DataFrame(columns=columns)

    # populateUntilT(T)
    #buildPredictedV()

    # add to dataframe risk class by risk class
    for t in range(0, T+1):
        sum = (s[k][t]+e[k][t]+i[k][t]+ud[k][t]+ur[k][t]+hd[k][t]+hr[k][t]+qd[k][t]+qr[k][t]+d[k][t]+r[k][t]+m[k][t])
        df.loc[t] = [t, s[k][t], e[k][t], i[k][t], ud[k][t], ur[k][t], hd[k][t], hr[k][t],
                     qd[k][t], qr[k][t], d[k][t], r[k][t], m[k][t], v[k][t], sum]
    return df

# build for 1 day for risk class 1
df = buildDF(1, 20)
# specify column 1 values are integers
df.to_excel(r'data.xlsx', index = False)
