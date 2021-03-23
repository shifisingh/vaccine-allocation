import pandas as pd

# data for each risk class
# S, E, I, U, H, Q, R, D
class1 = [.840, .840, .048, .0402, .000009233, .00803, .063, .000007254]
class2 = []
class3 = []
class4 = []
class5 = []
# 2d array that stores all risk class' data
data = [class1, class2, class3, class4, class5]

riskClasses = ['0-19', '20-29', '30-49', '50-69', '70+']

# for a risk class k, this may need to become a class
beta = .95
gamma = .5
alpha = .83
t = 30
vk = [10] * t
# S, E, I, U, H, Q, R, D
initValues = [.1]*8
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

def tPlusOne(initValues, t):
    # newvalues will be a list which in order produces : S, R, I, UD, UR, HD, HR,
    # QD, QR, D, M (not R though - this may need to be added in perhaps by
    # summing UR, HR, and QR)
    newValues = [0]*11

    # Susceptible people, equation 16
    Skt = initValues[0]
    newValues[0] = Skt - beta*vk[t] - ((alpha*gamma)*(Skt - beta*vk[t]))

    # Exposed people, equation 17
    Ekt = initValues[1]
    newValues[1] = Ekt + ((alpha*gamma*(Ekt - beta*vk[t]))*ksum - (ri*Ekt))

    # Infected people, equation 18
    Ikt = initValues[2]
    newValues[2] = Ikt + (ri*Ekt - rDetected*Ikt)

    # Undiagnosed people that die, equation 19
    UDkt = initValues[3] * (1 - rr)
    rud = UDkt / initValues[3]
    newValues[3] = UDkt + (rud*Ikt - rDeath*UDkt)
    # Undiagnosed people that recover
    URkt = initValues[3] * rr
    rur = URkt / initValues[3]
    newValues[4] = URkt + (rur*Ikt - rr*URkt)

    # Hospitalized people that die, equation 20
    HDkt = initValues[4] * (1 - rh)
    rhd = HDkt / initValues[4]
    newValues[5] = HDkt + (rhd*Ikt - rDeath*HDkt)
    # Hospitalized people that recover
    HRkt = initValues[4] * rh
    rhr = HRkt / initValues[4]
    newValues[6] = HRkt + (rhr*Ikt - rr*HRkt)

    # Quarantined people that die, equation 21
    QDkt = initValues[5] * .97
    rqd = QDkt / initValues[5]
    newValues[7] = QDkt + (rqd*Ikt - rDeath*QDkt)
    # Quarantined people that recover
    QRkt = initValues[5] * .03
    rqr = QRkt / initValues[5]
    newValues[8] = QRkt + (rqr*Ikt - rr*QRkt)

    # People that die, equation 22
    Dkt = initValues[6]
    newValues[9] = (Dkt+rDeath)*(UDkt + HDkt + QDkt)

    # Mkt, equation 23
    Mkt = initValues[7]
    newValues[10] = Mkt + beta*vk[t]

    return newValues

def populateUntilT(t,k):
    i = 0
    listOfList = [[] for i in range(t)]
    while i < t:
        if i == 0:
            listOfList[i] = tPlusOne(k, i)
        else:
            listOfList[i] = tPlusOne(k, i)
        i += 1
    return listOfList

def buildDF(t, k):
    # set the column names
    # S, R, I, UD, UR, HD, HR, QD, QR, D, M
    columns = ['time since day 0', 'susceptible', 'recovered', 'infected',
               'UD', 'UR', 'HD', 'HR', 'QD', 'QR', 'died', 'M']
    df = pd.DataFrame(columns = columns)

    # add to dataframe risk class by risk class
    data = populateUntilT(t, k)
    for i in range(len(data)):
        row = data[i]
        df.loc[i] = [i, row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9]]
    return df

# build for 10 days for risk class 1
buildDF(10, 1)




