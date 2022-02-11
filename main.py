import pandas as pd

# import functions and variables from other files
from variables import beta, gamma, alpha, ksum, ri, rDeath, rDetected, rr, rh, s, e, i, ud, ur, hd, hr, qd, qr, d, m, \
    mMin, pDetected, pHospitalized, rDetected
from variables import *
from vaccinations import v
from vaccinations import *
from death_rates import du, dh, dq

# initial data for each risk class
# S, E, I, UD, UR, HD, HR, QD, QR, R, D
class1 = [1274292 / 1516350, 1274292 / 1516350, 73164 / 1516350,
          (60970/1516350) * du[0], (60970/1516350) * (1-du[0]),
          (14/1516350) * dh[0], (14/1516350) * (1-dh[0]),
          (12180/1516350) * (dq[0]), (12180/1516350) * (1-dq[0]),
          (11/1516350), (95719/1516350)]
class2 = [676765/964950, 676765/964950, 87102/964950,
          72585/964950 * (du[1]), 72585/964950 * (1-du[1]),
          14/964950 * (dh[1]), 14/964950 * (1-dh[1]),
          14503/964950 * (dq[1]), 14503/964950 * (1-dq[1]),
          28/964950, 113953/964950]
class3 = [1278927/1723126, 1278927/1723126, 134136/1723126,
          111780/1723126 * (du[2]),  111780/1723126 * (1-du[2]),
          28/1723126 * (dh[2]), 28/1723126 * (1-dh[2]),
          22328/1723126 * (dq[2]), 22328/1723126 * (1-dq[2]),
          440/1723126, 175487/1723126]
class4 = [1440307/1792051, 1440307/1792051, 105396/1792051,
          87830/1792051 * (du[3]), 87830/1792051 * (1-du[3]),
          140/1792051 * (dh[3]), 140/1792051 * (1-dh[3]),
          17426/1792051 * (dq[3]), 17426/1792051 * (1-dq[3]),
          3066/1792051, 137886/1792051]
class5 = [632134/758175, 632134/758175, 35712/758175,
          29760/758175 * (du[4]), 29760/758175 * (1-du[4]),
          234/758175 * (dh[4]), 234/758175 * (1-dh[4]),
          5718/758175 * (dq[4]), 5718/758175 * (1 - dq[4]),
          7896/758175, 46721/758175]
rDetected = .15

print('du ', du)
print('dh ', dh)
print('dq ', dq)

#find iud, ihd, iqd rates
#disagreggated rates in appendix page 26

# 2d array that stores all risk class' day 0 data
data = [class1, class2, class3, class4, class5]

def tPlusOne(k, t):
    #global s, e, i, ud, ur, hd, hr, qd, qr, d, m, v

    # compute ksum
    ksum = i[0][t-1] + i[1][t-1] + i[2][t-1] + i[3][t-1] + i[4][t-1]

    # Susceptible people, equation 16
    Skt = s[k][t-1]
    s[k][t] = Skt - beta*v[k][t] - ((alpha*gamma)*(Skt - beta*v[k][t])*ksum)

    # Exposed people, equation 17
    Ekt = e[k][t-1]
    e[k][t] = Ekt + ((alpha*gamma*(Skt - beta*v[k][t]))*ksum - (ri*Ekt))

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

    # Mkt, equation 23
    Mkt = m[k][t-1]
    m[k][t] = Mkt + beta*v[k][t]

def populateUntilT(T):
    global data, s, e, i, ud, ur, hd, hr, qd, qr, d, m, v

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
                #r[k][0] = data[k][6]
                d[k][0] = data[k][10]
                m[k][0] = 0
            else:
                tPlusOne(k, t)

def buildDF(k, T):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m

    populateUntilT(T)

    # set the column names
    # S, R, I, UD, UR, HD, HR, QD, QR, D, M
    columns = ['time since day 0', 'susceptible', 'infected',
               'UD', 'UR', 'HD', 'HR', 'QD', 'QR', 'died', 'immune state']
    df = pd.DataFrame(columns=columns)

    # populateUntilT(T)
    #buildPredictedV()

    # add to dataframe risk class by risk class
    for t in range(0, T+1):
        df.loc[t] = [t, s[k][t], i[k][t], ud[k][t], ur[k][t], hd[k][t], hr[k][t],
                     qd[k][t], qr[k][t], d[k][t], m[k][t]]
    return df

# build for 1 day for risk class 1
print(buildDF(1, 1))
df = buildDF(1, 20)
# specify column 1 values are integers
df.to_excel(r'data.xlsx', index = False)
