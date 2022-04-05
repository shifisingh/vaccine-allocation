# population counts
class1population = 1516350
class2population = 964950
class3population = 1723126
class4population = 1792051
class5population = 758175
population = class1population + class2population + class3population + class4population + class5population
# for a risk class k, this may need to become a class
beta = .95 # vaccine efficacy
gamma = .5 # government intervention
# alpha = .83 # nominal infection rate
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
# minimum mortality rate (assume = mortality rate for now)
mMin = .01
# percentage of infectious cases detected
pDetected = .20
# percentage of detected cases hospitalized
pHospitalized = .15

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

