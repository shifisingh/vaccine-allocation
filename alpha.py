import pandas as pd
from sklearn.linear_model import LinearRegression
from variables import *

# summation equation adding up all infected cases across all age groups on a particular week
def summation(time):
    df0 = pd.read_csv('risk_class_0.csv')
    df1 = pd.read_csv('risk_class_1.csv')
    df2 = pd.read_csv('risk_class_2.csv')
    df3 = pd.read_csv('risk_class_3.csv')
    df4 = pd.read_csv('risk_class_4.csv')

    infected = (df0.loc[:, 'active cases'][time] + df0.loc[:, 'undetected'][time]) \
               + (df1.loc[:, 'active cases'][time] + df1.loc[:, 'undetected'][time]) \
               + (df2.loc[:, 'active cases'][time] + df2.loc[:, 'undetected'][time]) \
               + (df3.loc[:, 'active cases'][time] + df3.loc[:, 'undetected'][time]) \
               + (df4.loc[:, 'active cases'][time] + df4.loc[:, 'undetected'][time])
    infected_sum = infected / population

    return infected_sum

def findAlpha(k, t):
    # for each class
    # y = [[Δs for each t]...]
    # which will be an array of size t-1 of arrays of size 1 because start on week 2 (calculating change)
    # x = alpha * gamma * s[t] * i[t]
    # incubation period is 7 days

    df = pd.read_csv('risk_class_' + str(k) + '.csv')
    infected_sum = summation(t)

    x = [0 for i in range(t)]
    y = [0 for i in range(t)]

    for time in range(0, t):
        infected_sum = summation(t)
        suceptible = df.loc[:, 'suceptible'][time]
        x[time] = [gamma * (suceptible / population) * infected_sum]

    # fill out first value
    # q for prof: y is a list of differences so len(x) = len(y) + 1
    # -> should we tack on a default val at the end or beg of y?
    y[0] = 0

    # start at week 2 to week t-1 (non inclusive)
    for time in range(1, t):
        prev_suceptible = df.loc[:, 'suceptible'].values[time - 1]
        current_suceptible = df.loc[:, 'suceptible'].values[time]
        y[time] = prev_suceptible / population - current_suceptible / population

    print(x)
    print(y)

    reg = LinearRegression(positive=True)
    reg.fit(x, y)

    # array of the coefficients
    # parse them out
    return reg.coef_

# nominal infection rate over all 15 weeks for risk class 4
alpha = findAlpha(4, 15)
print(alpha)