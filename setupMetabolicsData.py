'''
Jon Stingel
10/09/2020
File to open and load all the metabolics data, experimental and simulation
into usable data structures for machine learning and other techniques. 

TODO: decide what form that file or structure should take. 
'''
import os
import pandas as pd
import glob
import math
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import statsmodels.api as sm

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score,mean_squared_error
import seaborn as sns


####################################################
## import the experimental values
# set paths
repobasedir = os.getcwd()
experimentalfile = os.path.join(repobasedir, 'experimentalMetabolics_all.csv')
# read in the file to dataframe
expmetcost_df = pd.read_csv(experimentalfile)
# print(expmetcost_df)
# print(expmetcost)
# f, ax = plt.subplots()
# expmetcost.plot(x="experiment [pi]",y="cost [W/kg]",title="plot of all exp costs",ax=ax)
# plt.show()

##################################################
## import all of the simulation results
# set all the paths
simresultspath = os.path.join(repobasedir,'..\\metabolicsResults\\')
muscleinversepath = os.path.join(simresultspath,'muscleInverse\\')
muscleInverseWithEMGpath = os.path.join(simresultspath,'muscleInverseWithEMG\\')

## first handle the values in the regular muscle driven inverse problem
# get all the filenames
musclefiles = glob.glob(os.path.join(muscleinversepath,'*.csv'))
# load them all into a single dataframe
df_from_each_file = (pd.read_csv(f) for f in musclefiles)
# print(df_from_each_file)
muscle_df = pd.concat(df_from_each_file, ignore_index=True)
# print(muscle_df)


test_df = muscle_df.groupby(['subjectname','condname','trialname']).agg({'metabolics_all_avg':['mean']})
test_df.columns = ['metabolics_all_avg_mean']
test_df = test_df.reset_index()
# print(test_df)

exp_df = expmetcost_df.groupby(['subjectname','condname','trialname']).agg({'metabolics_all_avg':['mean']})
exp_df.columns = ['metabolics_all_avg_mean']
exp_df = exp_df.reset_index()
# print(exp_df)


both_df = pd.merge(test_df, exp_df, how='right', on=['subjectname','condname','trialname'])
pd.set_option('display.max_rows',None) #,'display.max_columns',None)

bothtrim_df = both_df.dropna()
# print(bothtrim_df)

# figure out how to split them
# bothtrim_df.plot(kind='scatter',x='metabolics_all_avg_mean_x',y='metabolics_all_avg_mean_y')
# plt.show()


## testing
print('testing')
X = bothtrim_df['metabolics_all_avg_mean_x'].values.reshape(-1,1)
# print('x')
# print(X)
Y = bothtrim_df['metabolics_all_avg_mean_y'].values.reshape(-1,1)
# print('y')
# print(Y)

reg = LinearRegression()
reg.fit(X,Y)

print("The linear model is: Y = {:.5} + {:.5}X".format(reg.intercept_[0], reg.coef_[0][0]))


predictions = reg.predict(X)

# fig = plt.figure()
# ax = fig.add_subplot(111)

# plt.figure(figsize=(8,8))
# plt.scatter(bothtrim_df['metabolics_all_avg_mean_x'],
#     bothtrim_df['metabolics_all_avg_mean_y'],
#     c='black')
# plt.plot(bothtrim_df['metabolics_all_avg_mean_x'],
#     predictions,
#     c='blue',
#     linewidth=2)
# ax.set_aspect('equal')
# plt.xlim((3,10))
# plt.ylim((3,10))
# plt.show()

########################################################
## random forest implementation

# X = bothtrim_df.metabolics_all_avg_mean_x.values.reshape(-1,1)
# Y = bothtrim_df.metabolics_all_avg_mean_y.values.reshape(-1,1)

# x_train, x_test, y_train, y_test = train_test_split(X,Y,test_size=0.1,random_state=12)

# # RandomForestRegModel = RandomForestRegressor()
# # RandomForestRegModel.fit(x_train, y_train)

# print(x_train)
# print(y_train)