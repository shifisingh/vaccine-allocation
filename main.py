import pandas as pd

# initial data for each risk class
# S, E, I, U, H, Q, R, D
class1 = [1274292/1516350, 1274292/1516350, 73164/1516350, 60970/1516350, 14/1516350, 12180/1516350,
          11/1516350, 95719/1516350]
class2 = [676765/964950, 676765/964950, 87102/964950, 72585/964950, 14/964950, 14503/964950,
          28/964950, 113953/964950]
class3 = [1278927/1723126, 1278927/1723126, 134136/1723126, 111780/1723126, 28/1723126, 22328/1723126,
          440/1723126, 175487/1723126]
class4 = [1440307/1792051, 1440307/1792051, 105396/1792051, 87830/1792051, 140/1792051, 17426/1792051,
          3066/1792051, 137886/1792051]
class5 = [632134/758175, 632134/758175, 35712/758175, 29760/758175, 234/758175, 5718/758175,
          7896/758175, 46721/758175]
# 2d array that stores all risk class' data
data = [class1, class2, class3, class4, class5]

# for a risk class k, this may need to become a class
beta = .95
gamma = .5
alpha = .83
# 2d array of number of vaccinations, state[k][t], where k represents the risk class and t represents the time
v = [0] * 1
# probability of getting into contact with an infected person
ksum = .1
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

def tPlusOne(k, t):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m

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
    ud[k][t] = URkt + (rur*Ikt - rr*URkt)

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

def populateUntilT(T):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m

    for k in range(0, 5):
        for t in range(0, T):
            if t == 0:
                # populate the preliminary data
                size = [[k for i in range(T+1)] for i in range(5)]
                s = e = i = ud = ur = hd = hr = qd = qr = d = m = size
                s[k][0] = data[k][0]
                e[k][0] = data[k][1]
                i[k][0] = data[k][2]
                ud[k][0] = 0
                ur[k][0] = 0
                hd[k][0] = 0
                hr[k][0] = 0
                qd[k][0] = 0
                qr[k][0] = 0
                #r[k][0] = data[k][6]
                d[k][0] = data[k][7]
                m[k][0] = 0
            else:
                tPlusOne(k, t)

def buildDF(k, T):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m

    # set the column names
    # S, R, I, UD, UR, HD, HR, QD, QR, D, M
    columns = ['time since day 0', 'susceptible', 'infected',
               'UD', 'UR', 'HD', 'HR', 'QD', 'QR', 'died', 'immune state']
    df = pd.DataFrame(columns=columns)

    populateUntilT(T)

    # add to dataframe risk class by risk class
    for t in range(0, T+1):
        df.loc[t] = [t, s[k][t], i[k][t], ud[k][t], ur[k][t], hd[k][t], hr[k][t],
                     qd[k][t], qr[k][t], d[k][t], m[k][t]]
    return df

# build for 1 day for risk class 1
print(buildDF(1, 1))