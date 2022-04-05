from math import floor

import pandas as pd

v = [[]]

# FILL BASED ON HISTORICAL DATA
# Populate vaccinations for a flat vax rate
def populateVaxByWeek(week):
    global v

    for k in [0, 1, 2, 3, 4]:
        df = pd.read_csv('risk_class_' + str(k) + '.csv')
        # for week 1, start at day 7 end before day 15
        for i in range((week*7)+1, (week*7)+8):
            # the proportion of vaccinated people
            rate = (df['vaccinated'][week * 7])/(df['suceptible'][week * 7])
            v[k][i] = v[k][i-1] * (1+rate)


def singleDoseVaxStrat(t):
    global v

    # vaccination hesitancy set to 30%
    hesitancy = .3

    #populateInitVax(70)

    # for every risk class
    for k in [0, 1, 2, 3, 4]:
        df = pd.read_csv('risk_class_' + str(k) + '.csv')
        rate = (df['vaccinated'][10])/(df['suceptible'][10])
        # loop over days now
        for time in range(71, t):
            v[k][time] = (v[k][time-1] * (1-hesitancy)) * (1+rate)
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

# STRATEGY 2
# Filling in initial values for a simulation
# We have 15 weeks of data but here, we are filling in the 'v' array for the first 10 weeks
# and seeing if we can correctly predict the next 5 weeks of data

def fillKnownValues():
    global v

    for k in [0, 1, 2, 3, 4]:
        df = pd.read_csv('risk_class_' + str(k) + '.csv')
        # fill up 70 days from day 0, until 2/27
        # start populating from day 1 b/c day 0 should have 0 vaccinations
        for i in range(1, 105):
            if v[k][i] == 0:
                # data is not cumulative so take the vax count for each week, divide by 7, to split between the week.
                v[k][i] = df['vaccinated'][(floor(i/7) + 1)] / 7

def buildPredictedV():
    global v

    # initialize for the 15 weeks (we have all this data but predicting 10 weeks)
    v = [[0 for i in range(105)] for i in range(5)]

    # fill in the in-between values
    fillKnownValues()

    return v;

    #print(v[0])
    #print(len(v[0]))

# buildPredictedV()

# NEW VAX STRATS
# fixed number of vaccines per week
