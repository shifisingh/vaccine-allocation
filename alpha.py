import pandas as pd
from sklearn.linear_model import LinearRegression
from variables import gamma


def findAlpha(k, t):
    # for each class
    # y = [[Δs for each t]...]
    # which will be an array of size t-1 of arrays of size 1 because start on week 2 (calculating change)
    # x = alpha * gamma * s[t] * i[t]
    # incubation period is 7 days

    df = pd.read_csv('risk_class_' + str(k) + '.csv')

    x = [0 for i in range(t)]
    y = [0 for i in range(t)]

    print((df.loc[:, 'active cases']))

    for time in range(0, t):
        x[time] = [gamma * df.loc[:, 'suceptible'][time] * (df.loc[:, 'active cases'][time]
                                                            + df.loc[:, 'undetected'][time])]

    for time in range(0, t - 1):
        y[time] = df.loc[:, 'suceptible'].values[time + 1] - df.loc[:, 'suceptible'].values[time]

    reg = LinearRegression(positive=True)
    reg.fit(x, y)

    # array of the coefficients
    # parse them out
    return reg.coef_

.0000083728
alpha = findAlpha(1, 2)
print(alpha)