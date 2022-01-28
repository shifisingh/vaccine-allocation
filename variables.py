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

