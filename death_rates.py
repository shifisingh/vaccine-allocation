import pandas as pd
from sklearn.linear_model import LinearRegression

def findDeathRates(k):
    # for each class
    # x = [[u0, h0, q0], [u1, h1, q1], [u2, h2, q2], [u3, h3, q3]]
    # which will be t arrays size 3 because each of the values will represent its states
    # y = death array

    df = pd.read_csv('risk_class_'+str(k)+'.csv')

    x_feat_list = ['undetected', 'hospitalized', 'quarantined']

    x = df.loc[:, x_feat_list].values
    y = df.loc[:,'total deaths'].values

    reg = LinearRegression(positive=True)
    reg.fit(x, y)

    # array of the coefficients
    # parse them out
    return reg.coef_

# death rate for undetected
du = [findDeathRates(0)[0], findDeathRates(1)[0], findDeathRates(2)[0],
      findDeathRates(3)[0], findDeathRates(4)[0]]
# death rate for hospitalized
dh = [findDeathRates(0)[1], findDeathRates(1)[1], findDeathRates(2)[1],
      findDeathRates(3)[1], findDeathRates(4)[1]]
# death rate for quarantined
dq = [findDeathRates(0)[2], findDeathRates(1)[2], findDeathRates(2)[2],
      findDeathRates(3)[2], findDeathRates(4)[2]]