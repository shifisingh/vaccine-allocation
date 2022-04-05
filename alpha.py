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

def findAlpha(k, start, end):
    # for each class
    # y = [[Δs for each t]...]
    # which will be an array of size t-1 of arrays of size 1 because start on week after start (calculating change)
    # x = alpha * gamma * s[t] * i[t]
    # incubation period is 7 days

    df = pd.read_csv('risk_class_' + str(k) + '.csv')

    x = []
    y = []

    for time in range(start, end):
        infected_sum = summation(time)
        # print(infected_sum)
        suceptible = df.loc[:, 'suceptible'][time]
        x.append([gamma * (suceptible / population) * infected_sum])

    # start at the week after start to end (non inclusive)
    for time in range(start + 1, end):
        prev_suceptible = df.loc[:, 'suceptible'].values[time - 1]
        current_suceptible = df.loc[:, 'suceptible'].values[time]
        y.append(prev_suceptible / population - current_suceptible / population)

    # fill out first value
    # q for prof: y is a list of differences so len(x) = len(y) + 1
    # -> should we tack on a default val at the end or beg of y?
    y.append(0)

    # print(x)
    # print(y)

    reg = LinearRegression(positive=True)
    reg.fit(x, y)

    # array of the coefficients
    # parse them out
    return reg.coef_

# nominal infection rate over all 15 weeks for risk class 4
alpha = findAlpha(4, 9, 15)
# print(alpha)