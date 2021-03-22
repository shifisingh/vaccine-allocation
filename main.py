
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

def populateUntilT(t):
    i = 0
    listOfList = [[] for i in range(t)]
    while i < t:
        if i == 0:
            listOfList[i] = tPlusOne(initValues, i)
        else:
            listOfList[i] = tPlusOne(listOfList[i - 1], i)
        i += 1
    return listOfList

print(populateUntilT(30))






