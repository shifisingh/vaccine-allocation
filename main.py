import pandas as pd
from sklearn.linear_model import LinearRegression

# for a risk class k, this may need to become a class
beta = .95
gamma = .5
alpha = .83
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

# initial data for each risk class
# S, E, I, UD, UR, HD, HR, QD, QR, R, D
class1 = [1274292 / 1516350, 1274292 / 1516350, 73164 / 1516350,
          (60970/1516350) * (1 - rr), (60970/1516350) * (rr),
          (14/1516350) * (1 - rh), (14/1516350) * (rh),
          (12180/1516350) * (1 - rr), (12180/1516350) * (rr),
          (11/1516350), (95719/1516350)]
class2 = [676765/964950, 676765/964950, 87102/964950,
          72585/964950 * (1 - rr), 72585/964950 * (rr),
          14/964950 * (1 - rh), 14/964950 * (rh),
          14503/964950 * (1 - rr), 14503/964950 * (rr),
          28/964950, 113953/964950]
class3 = [1278927/1723126, 1278927/1723126, 134136/1723126,
          111780/1723126 * (1 - rr),  111780/1723126 * (rr),
          28/1723126 * (1 - rh), 28/1723126 * (rh),
          22328/1723126 * (1 - rr), 22328/1723126 * (rr),
          440/1723126, 175487/1723126]
class4 = [1440307/1792051, 1440307/1792051, 105396/1792051,
          87830/1792051 * (1 - rr), 87830/1792051 * (rr),
          140/1792051 * (1 - rh), 140/1792051 * (rh),
          17426/1792051 * (1 - rr), 17426/1792051 * (rr),
          3066/1792051, 137886/1792051]
class5 = [632134/758175, 632134/758175, 35712/758175,
          29760/758175 * (1 - rr), 29760/758175 * (rr),
          234/758175 * (1 - rh), 234/758175 * (rh),
          5718/758175 * (1 - rr), 5718/758175 * (1 - rr),
          7896/758175, 46721/758175]
# 2d array that stores all risk class' data
data = [class1, class2, class3, class4, class5]

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

def findDeathRates():
    # for each class
    # x = [[u0, h0, q0], [u1, h1, q1], [u2, h2, q2], [u3, h3, q3]]
    # which will be t arrays size 3 because each of the values will represent its states
    # y = death array

    x1 = [[60970, 14, 12180]]
    x2 = [[72585, 14, 14503]]
    x3 = [[11780, 28, 22328]]
    x4 = [[87830, 140, 17426]]
    x5 = [[29760, 234, 5718]]
    y1 = [[11]]
    y2 = [[28]]
    y3 = [[440]]
    y4 = [[3066]]
    y5 = [[7896]]

    reg = LinearRegression()
    reg.fit(x1, y1)

    # array of the coefficients
    # parse them out
    return reg.coef_

def findVaxRates(df):
    # given a dataframe of vaccination data, figure out the rates for which people
    # are getting vaccinated during each phase

    # find the x & y
    x = df.loc[:, ['risk_class_1', 'risk_class_2', 'risk_class_3', 'risk_class_4', 'risk_class_5']].values.reshape(-1, 1)
    y = df.loc[:, 'total'].values

    # do regression on the given x & y values
    reg = LinearRegression()
    reg.fit(x, y)

    # returns an array holding the coefficient which is parsed
    return reg.coef_

def singleDoseVaxStrat(t):
    # build the 2d array so v[k][t] is the total number of people vaccinated at time t in risk class k
    # t arrays of size 5
    v = [[0 for i in range(t + 1)] for i in range(5)]
    p = [[0 for i in range(t + 1)] for i in range(5)]

    # vaccination hesitancy set to 30%
    hesitancy = .3

    # import the csv as a dataframe
    phase_1_df = pd.read_csv('phase_1.csv')

    for time in t:
        for k in ['risk_class_1', 'risk_class_2', 'risk_class_3', 'risk_class_4', 'risk_class_5']:
            if time > 0 & time <= 42:
                rate = findVaxRates(phase_1_df, k)
                v[k][time] = (v[k][time-1] * (1-hesitancy)) * rate
            # phase 2
            if time >= 43 & time <= 62:
                # rate = findVaxRates(df, k)
                v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate
            # phase 3
            if time >= 63 & time <= 76:
                # rate = findVaxRates(df, k)
                v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate
            # phase 4
            if time >= 77:
                # rate = findVaxRates(df, k)
                v[k][time] = (v[k][time - 1] * (1 - hesitancy)) * rate
    return v

phase_1_df = pd.read_csv('phase_1.csv')
# x = phase_1_df.loc[:, ('risk_class_1', 'risk_class_2', 'risk_class_3', 'risk_class_4', 'risk_class_5')].values
x = [[1274292, 676765, 1278927, 1440307, 632134]]
y = [35618]
#y = phase_1_df.loc[:, 'total'].values[0]
reg = LinearRegression()
reg.fit(x, y)
print(reg.coef_)

def tPlusOne(k, t):
    global s, e, i, ud, ur, hd, hr, qd, qr, d, m, v

    # compute ksum
    ksum = i[0][t] + i[1][t] + i[2][t] + i[3][t] + i[4][t]

    # Vax strategy
    v = [[0 for i in range(t + 1)] for i in range(5)]
    # add new equations to support our vax strat

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