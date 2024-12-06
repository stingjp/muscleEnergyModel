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

# from sklearn.linear_model import LinearRegression
# from sklearn.metrics import r2_score
# import statsmodels.api as sm

# from sklearn.model_selection import train_test_split
# from sklearn.ensemble import RandomForestRegressor
# from sklearn.metrics import r2_score,mean_squared_error
import seaborn as sns
import pdb
import matplotlib.lines as mlines
import matplotlib.transforms as mtransforms
import seaborn as sns
from pylab import plot, show, savefig, xlim, figure, ylim, legend, boxplot, setp, axes #hold


####################################################
## import the experimental values
# set paths
repobasedir = os.getcwd()
experimentalfile = os.path.join(repobasedir, 'experimentalMetabolics_all.csv')
# read in the file to dataframe
expmetcost_df = pd.read_csv(experimentalfile)
print(expmetcost_df)
# print(expmetcost)
# f, ax = plt.subplots()
# expmetcost.plot(x="experiment [pi]",y="cost [W/kg]",title="plot of all exp costs",ax=ax)
# expmetcost_df.plot(y="metabolics_all_avg",title="plot of all exp costs",ax=ax, marker='o', alpha=0.7)
# plt.show()
# import pdb
# pdb.set_trace()


##################################################
## import all of the simulation results
# set all the paths
simresultspath = os.path.join(repobasedir,'..\\metabolicsResults\\')
muscleinversepath = os.path.join(simresultspath,'muscleInverse\\')
muscleInverseWithEMGpath = os.path.join(simresultspath,'muscleInverseWithEMG\\')
muscleMetabolicspath = os.path.join(simresultspath,'muscleInverse_100con\\')
# muscletestpath = 'G:\\Shared drives\\Exotendon\\muscleModel\\metabolicsResults\\muscleInverse_100con\\'
# musclestancepath = 'G:\\Shared drives\\Exotendon\\muscleModel\\metabolicsResults\\muscleStance\\'
# muscleswingpath = 'G:\\Shared drives\\Exotendon\\muscleModel\\metabolicsResults\\muscleSwing\\'
muscletestpath = os.path.join(simresultspath,'muscleInverse_100con\\')
musclestancepath = os.path.join(simresultspath,'muscleStance\\')
muscleswingpath = os.path.join(simresultspath,'muscleSwing\\')


## first handle the values in the regular muscle driven inverse problem
# get all the filenames
musclefiles = glob.glob(os.path.join(muscletestpath,'*.csv'))

# load them all into a single dataframe
df_from_each_file = (pd.read_csv(f) for f in musclefiles)
print(df_from_each_file)
muscle_df = pd.concat(df_from_each_file, ignore_index=True)
print(muscle_df)

# testcode to get the full set of stance and swing data all condensed into one dataframe each...
stancepath = glob.glob(os.path.join(musclestancepath,'*.csv'))
swingpath = glob.glob(os.path.join(muscleswingpath,'*.csv'))
# load all the stance files into a dataframe
df_stance_files = (pd.read_csv(f) for f in stancepath)
stance_df = pd.concat(df_stance_files, ignore_index=True)
print(stance_df)

# load all the swing files into a dataframe
df_swing_files = (pd.read_csv(f) for f in swingpath)
swing_df = pd.concat(df_swing_files, ignore_index=True)
print(swing_df)

## working on formatting the dataframe for an excel sheet where I do all the muscle combinations
stance1 = stance_df.groupby(['Var2','subjectname','condname','trialname']).agg({'Var1':['mean']})
stance1.columns = ['metabolics_stance_avg']
stance1 = stance1.reset_index()
# Filter out rows where 'Var2' does not contain 'TOTAL'
stance2 = stance1[stance1['Var2'].str.contains('metabolics_combined')]
stance2 = stance2.reset_index()
stance3 = stance2.groupby(['Var2','condname','subjectname']).agg({'metabolics_stance_avg':['mean']})
stance3.columns = ['metabolics_stance_avg']
stance3 = stance3.reset_index()
# Save stance3 to a CSV file
stance3.to_csv('stance3_output.csv', index=False)
# Pivot the stance3 dataframe to get the desired orientation
stance3_pivot = stance3.pivot(index='Var2', columns=['condname','subjectname'], values='metabolics_stance_avg')
# Save the pivoted dataframe to a CSV file
stance3_pivot.to_csv('stance3_pivot_output.csv')

# now do the same for the swing data
swing1 = swing_df.groupby(['Var2','subjectname','condname','trialname']).agg({'Var1':['mean']})
swing1.columns = ['metabolics_swing_avg']
swing1 = swing1.reset_index()
# Filter out rows where 'Var2' does not contain 'TOTAL'
swing2 = swing1[swing1['Var2'].str.contains('metabolics_combined')]
swing2 = swing2.reset_index()
swing3 = swing2.groupby(['Var2','condname','subjectname']).agg({'metabolics_swing_avg':['mean']})
swing3.columns = ['metabolics_swing_avg']
swing3 = swing3.reset_index()
# Save swing3 to a CSV file
swing3.to_csv('swing3_output.csv', index=False)
# Pivot the swing3 dataframe to get the desired orientation
swing3_pivot = swing3.pivot(index='Var2', columns=['condname','subjectname'], values='metabolics_swing_avg')
# Save the pivoted dataframe to a CSV file
swing3_pivot.to_csv('swing3_pivot_output.csv')





pdb.set_trace()
### working on organizing the data frame - test full first
print(muscle_df)

# get the full metabolics dataframe
test_df = muscle_df.groupby(['subjectname','condname','trialname']).agg({'metabolics_all_avg':['mean']})
print(test_df)
test_df.columns = ['metabolics_all_avg_mean']
print(test_df)
test_df = test_df.reset_index()
print('test_df: full metabolics df')
print(test_df)
pdb.set_trace()




# get the swing metabolics dataframe
swing_df = muscle_df.groupby(['subjectname','condname','trialname']).agg({'metabolics_swing_avg':['mean']})
swing_df.columns = ['metabolics_swing_avg_mean']
swing_df = swing_df.reset_index()
print('swing_df')
print(swing_df)
# get the stance metabolics dataframe
stance_df = muscle_df.groupby(['subjectname','condname','trialname']).agg({'metabolics_stance_avg':['mean']})
stance_df.columns = ['metabolics_stance_avg_mean']
stance_df = stance_df.reset_index()
print('stance_df')
print(stance_df)

# get the experimental metabolics dataset in the same form - averaged across subject conditions and trials. 
exp_df = expmetcost_df.groupby(['subjectname','condname','trialname']).agg({'metabolics_all_avg':['mean']})
exp_df.columns = ['metabolics_all_avg_mean']
exp_df = exp_df.reset_index()
print('exp_df')
print(exp_df)


both_df = pd.merge(test_df, exp_df, how='right', on=['subjectname','condname','trialname'])
pd.set_option('display.max_rows',None) #,'display.max_columns',None)

bothtrim_df = both_df.dropna()
print(bothtrim_df)


### going to make a print out of the actual reductions for swing and stance, as well as percent 
# grab the raw differences between them
stance_means = stance_df.groupby(['condname']).agg({'metabolics_stance_avg_mean':['mean']})
swing_means = swing_df.groupby(['condname']).agg({'metabolics_swing_avg_mean':['mean']})

print(stance_means)
print(swing_means)

# get all the values in a workable format
# swings_exo = swing_df.loc[swing_df['condname'] == 'welkexo']
# swings_natural = swing_df.loc[swing_df['condname'] == 'welknatural']
# stances_exo = stance_df.loc[stance_df['condname'] == 'welkexo']
# stances_natural = stance_df.loc[stance_df['condname'] == 'welknatural']
### edits on 20241105 - pretty sure these are what we were going for. 
swings_exo = swing_means.loc[swing_means.index == 'welkexo']
swings_natural = swing_means.loc[swing_means.index == 'welknatural']
stances_exo = stance_means.loc[stance_means.index == 'welkexo']
stances_natural = stance_means.loc[stance_means.index == 'welknatural']




# get the average raw differences for both stance and swing - check with above
stance_change = np.mean(stances_exo) - np.mean(stances_natural)
swing_change = np.mean(swings_exo) - np.mean(swings_natural)
# now to get the percent changes for stance and swing
stance_perc_change = (np.mean(stances_exo) - np.mean(stances_natural)) / np.mean(stances_natural) * 100
swing_perc_change = (np.mean(swings_exo) - np.mean(swings_natural)) / np.mean(swings_natural) * 100

print("Raw stance difference: %f" % stance_change)
print('Raw swing  difference: %f' % swing_change)
print('Percent Difference stance: %f' % stance_perc_change)
print('Percent Difference swing: %f' % swing_perc_change)



# make a figure for the simulations vs exp
f, ax = plt.subplots()
tempx = bothtrim_df['metabolics_all_avg_mean_x']
tempy = bothtrim_df['metabolics_all_avg_mean_y']
mse = mean_squared_error(tempy, tempx)
ax.scatter(tempx, tempy, marker='o', alpha=0.6, c='blue', label='MSE = %f'%mse)
# line = mlines.Line2D([0,1], [0,1], color='red')
# transform = ax.transAxes
# line.set_transform(transform)
# ax.add_line(line)
unitx = np.linspace(8,12,100)
ax.plot(unitx,unitx,color='red',alpha=0.5, label='y = x')
# some nice 
ax.set_title('Simulation vs Experimental')
plt.xlim([8, 13])
plt.ylim([8, 13])
ax.set_aspect('equal', adjustable='box')
plt.grid()
plt.xlabel('Simulated [W/kg]')
plt.ylabel('Experimental [W/kg]')
plt.legend()
# plt.show()

# make a figure for exo vs natural stance and swing
# thinking box plots
# pdb.set_trace()

swingexo = np.array([])
swingnatural = np.array([])
stanceexo = np.array([])
stancenatural = np.array([])

# split the swings up
for i, row in swing_df.iterrows():
    # print(i)
    # print(row)
    tempcond = row['condname']
    if 'welkexo' in tempcond:
        swingexo = np.append(swingexo, row['metabolics_swing_avg_mean'])
    elif 'welknatural' in tempcond:
        swingnatural = np.append(swingnatural, row['metabolics_swing_avg_mean'])
for i, row in stance_df.iterrows():
    # print(i)
    # print(row)
    tempcond = row['condname']
    if 'welkexo' in tempcond:
        stanceexo = np.append(stanceexo, row['metabolics_stance_avg_mean'])
    elif 'welknatural' in tempcond:
        stancenatural = np.append(stancenatural, row['metabolics_stance_avg_mean'])

# now have 4 vectors, just need to get stats and put in a bar plot

# function for setting the colors of the box plots pairs
def setBoxColors(bp):

    setp(bp['boxes'][1], color='blue')
    setp(bp['caps'][2], color='blue')
    setp(bp['caps'][3], color='blue')
    setp(bp['whiskers'][2], color='blue')
    setp(bp['whiskers'][3], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][1], color='blue')

    setp(bp['boxes'][0], color='red')
    setp(bp['caps'][0], color='red')
    setp(bp['caps'][1], color='red')
    setp(bp['whiskers'][0], color='red')
    setp(bp['whiskers'][1], color='red')
    setp(bp['fliers'][0], color='red')
    setp(bp['fliers'][1], color='red')
    setp(bp['medians'][0], color='red')

swings = [swingnatural, swingexo]
stances = [stancenatural, stanceexo]


# pdb.set_trace()
fig = figure()
ax = axes()
#hold(True)

# first pair - stance
bp = boxplot(stances, positions=[1,2], widths=1)
setBoxColors(bp)
# second pair - swing
bp = boxplot(swings, positions=[4,5], widths=1)
setBoxColors(bp)
# set axes limits and labels
xlim(0,6)
ylim(2,9)


ax.set_xticks([1.5, 4.5])
# ax.set_xticklabels(['Stance\n13% Reduction Avg.', 'Swing\n1% Increase Avg.'],fontsize=16)
ax.set_xticklabels(['Stance', 'Swing'],fontsize=16)
# ax.yticks(fontsize=16)
# draw temporary red and blue lines and use them to create a legend
hB, = plot([1,1],'r-')
hR, = plot([1,1],'b-')
plt.legend((hB, hR),('Natural', 'Exo'),fontsize=16)
plt.ylabel('Metabolic Cost [W/kg]',fontsize=16)
plt.yticks(fontsize=16)
hB.set_visible(False)
hR.set_visible(False)
plt.grid()
colorsjit = ['red', 'blue']
for i in [1,2]: # ,4,5]:
    y = stances[i-1] # titanic.age[titanic.pclass==i].dropna()
    # Add some random "jitter" to the x-axis
    x = np.random.normal(i, 0.0, size=len(y))
    plot(x, y, '.', alpha=0.8, color='grey') # color=colorsjit[i-1])
for i in [4,5]:
    y = swings[i-4]
    x = np.random.normal(i, 0.0, size=len(y))
    plot(x, y, '.', alpha=0.8, color='grey') # color=colorsjit[i-4])

# figure out how to put in subject lines


# plt.text(0.5,8.5,'Mean: %.1f' % np.mean(stances_natural),fontsize=16)
# plt.text(1.75,8.5,'Mean: %.1f' % np.mean(stances_exo),fontsize=16)
# plt.text(3.5,5.5,'Mean: %.1f' % np.mean(swings_natural),fontsize=16)
# plt.text(4.75,5.5,'Mean: %.1f' % np.mean(swings_exo),fontsize=16)

show()
pdb.set_trace()
pdb.set_trace()



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