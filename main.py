import pandas as pd
from sklearn.linear_model import LinearRegression

# for a risk class k, this may need to become a class
beta = .95 # vaccine efficacy
gamma = .5 # government intervention
alpha = .83 # nominal infection rate
# 2d array of number of vaccinations, vaccinated[k][t], where k represents the risk class and t represents the time
v = [[]]
# probability of getting into contact with an infected person (initialized to 0)
ksum = 0
# progression rate

ri = .06
# death rate, r^D
rDeath = 0.028
# detection rate, r^d
rDetected = .15
# recovery rate for non hospitalized
rr = .03
# recovery rate for hospitalized
rh = .02

# 2d array of states, state[k][t], where k represents the risk class and t represents the time
s = [[]]
e = [[]]
i = [[]]
ud = [[]]
ur = [[]]
hd = [[]]
hr = [[]]
qd = [[]]
qr = [[]]
#r = [[]]
d = [[]]
m = [[]]

def findDeathRates(k):
    # for each class
    # x = [[u0, h0, q0], [u1, h1, q1], [u2, h2, q2], [u3, h3, q3]]
    # which will be t arrays size 3 because each of the values will represent its states
    # y = death array

    df = pd.read_csv('risk_class_'+str(k)+'.csv')

    x_feat_list = ['undetected', 'hospitalized', 'quarantined']

    x = df.loc[:, x_feat_list].values
    y = df.loc[:,'total deaths'].values

    reg = LinearRegression()
    reg.fit(x, y)

    # array of the coefficients
    # parse them out
    return reg.coef_

# death rate for undetected
du = [findDeathRates(1)[0], findDeathRates(2)[0], findDeathRates(3)[0], findDeathRates(4)[0], findDeathRates(5)[0]]
# death rate for hospitalized
dh = [findDeathRates(1)[1], findDeathRates(2)[1], findDeathRates(3)[1], findDeathRates(4)[1], findDeathRates(5)[1]]
# death rate for quarantined
dq = [findDeathRates(1)[2], findDeathRates(2)[2], findDeathRates(3)[2], findDeathRates(4)[2], findDeathRates(5)[2]]

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
# 2d array that stores all risk class' data
data = [class1, class2, class3, class4, class5]

def singleDoseVaxStrat(t):
    # build the 2d array so v[k][t] is the total number of people vaccinated at time t in risk class k
    # t arrays of size 5
    # initialize values to 0
    v = [[0 for i in range(t + 1)] for i in range(5)]

    # vaccination hesitancy set to 30%
    hesitancy = .3

    for time in t:
        for k in ['risk_class_1', 'risk_class_2', 'risk_class_3', 'risk_class_4', 'risk_class_5']:
            if time > 0 & time <= 42:
                # rate =
                v[k][time] = (v[k][time-1] * (1-hesitancy)) * rate
            # phase 2
            if time >= 43 & time <= 62:
                # rate =
                v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate
            # phase 3
            if time >= 63 & time <= 76:
                # rate =
                v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate
            # phase 4
            if time >= 77:
                # rate =
                v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate
    return v

def flatVaxStrat(t):
    v = [[0 for i in range(t + 1)] for i in range(5)]

    hesitancy = .3
    rate = .02

    for time in range(1, t+1):
        for k in [0, 1, 2, 3, 4]:
            v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate

    return v

def buildV(strat,t):
    global v
    v = [[0 for i in range(t + 1)] for i in range(5)]

    if strat == 'flat vax':
        v = flatVaxStrat(t)
    if strat == 'single dose':
        v = singleDoseVaxStrat(t)

    return v

def tPlusOne(k, t):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m, v

    # compute ksum
    ksum = i[0][t] + i[1][t] + i[2][t] + i[3][t] + i[4][t]

    # Susceptible people, equation 16
    Skt = s[k][t-1]
    s[k][t] = Skt - beta*v[k][t] - ((alpha*gamma)*(Skt - beta*v[k][t])*ksum)

    # Exposed people, equation 17
    Ekt = e[k][t-1]
    e[k][t] = Ekt + ((alpha*gamma*(Skt - beta*v[k][t]))*ksum - (ri*Ekt))

    # Infected people, equation 18
    Ikt = i[k][t-1]
    i[k][t] = Ikt + (ri*Ekt - rDetected*Ikt)

    # Undiagnosed people that die, equation 19
    UDkt = ud[k][t-1] * (1 - rr)
    rud = UDkt / ud[k][t-1]
    ud[k][t] = UDkt + (rud*Ikt - rDeath*UDkt)
    # Undiagnosed people that recover
    URkt = ud[k][t-1] * rr
    rur = URkt / ud[k][t-1]
    ur[k][t] = URkt + (rur*Ikt - rr*URkt)

    # Hospitalized people that die, equation 20
    HDkt = hd[k][t-1] * (1 - rh)
    rhd = HDkt / hd[k][t-1]
    hd[k][t] = HDkt + (rhd*Ikt - rDeath*HDkt)
    # Hospitalized people that recover
    HRkt = hr[k][t-1] * rh
    rhr = HRkt / hr[k][t-1]
    hr[k][t] = HRkt + (rhr*Ikt - rr*HRkt)

    # Quarantined people that die, equation 21
    QDkt = qd[k][t-1] * (1 - rr)
    rqd = QDkt / qd[k][t-1]
    qd[k][t] = QDkt + (rqd*Ikt - rDeath*QDkt)
    # Quarantined people that recover
    QRkt = qr[k][t-1] * rr
    rqr = QRkt / qr[k][t-1]
    qr[k][t] = QRkt + (rqr*Ikt - rr*QRkt)

    # People that die, equation 22
    Dkt = d[k][t-1]
    d[k][t] = (Dkt+rDeath)*(UDkt + HDkt + QDkt)

    # Mkt, equation 23
    Mkt = m[k][t-1]
    m[k][t] = Mkt + beta*v[k][t]

def populateUntilT(T, strat):
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

    v = buildV(strat, T)

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

def buildDF(k, T, strat):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m

    # set the column names
    # S, R, I, UD, UR, HD, HR, QD, QR, D, M
    columns = ['time since day 0', 'susceptible', 'infected',
               'UD', 'UR', 'HD', 'HR', 'QD', 'QR', 'died', 'immune state']
    df = pd.DataFrame(columns=columns)

    populateUntilT(T, strat)

    # add to dataframe risk class by risk class
    for t in range(0, T+1):
        df.loc[t] = [t, s[k][t], i[k][t], ud[k][t], ur[k][t], hd[k][t], hr[k][t],
                     qd[k][t], qr[k][t], d[k][t], m[k][t]]
    return df

# build for 1 day for risk class 1
print(buildDF(1, 1, 'flat vax'))