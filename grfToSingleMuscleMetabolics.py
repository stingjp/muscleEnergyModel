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
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
import seaborn as sns
import scipy


from keras import layers
from keras.layers import Input, Dense, BatchNormalization, Dropout #, Activation, ZeroPadding2D, , Flatten, Conv2D
# from keras.layers import AveragePooling2D, MaxPooling2D, Dropout, GlobalMaxPooling2D, GlobalAveragePooling2D
from keras.regularizers import l2

from keras.models import Model, Sequential
from keras.wrappers.scikit_learn import KerasRegressor
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline



# from keras.preprocessing import image
from keras.utils import layer_utils
from keras.utils.data_utils import get_file
# from keras.applications.imagenet_utils import preprocess_input
import pydot
from IPython.display import SVG
from keras.utils.vis_utils import model_to_dot
from keras.utils import plot_model
# from kt_utils import *

import keras.backend as K
# K.set_image_data_format('channels_last')
from matplotlib.pyplot import imshow

# %matplotlib inline
# pd.set_option('display.max_rows',None,'display.max_columns',None)

#################################################################
# functions
#################################################################

# create the base NN model
def baseline_model(input_shape):
    
    # # create the model
    # model = Sequential()
    # model.add(Dense(5, input_dim=5, kernel_initializer='normal', activation='relu'))
    # model.add(Dense(4, kernel_initializer='normal'))
    # model.add(Dense(2, kernel_initializer='normal'))
    # model.add(Dense(1, kernel_initializer='normal'))
    # # compile model
    # model.compile(loss='mean_squared_error', optimizer='adam')
    
    X_input = Input(input_shape)
    # X = Dropout(0.2, input_shape=input_shape)(X_input)
    X = Dense(100, input_dim=input_shape, kernel_initializer='normal', activation='relu', kernel_regularizer=l2(0.01))(X_input)
    X = Dropout(0.2)(X)
    X = Dense(100, activation='relu', kernel_initializer='normal', kernel_regularizer=l2(0.01))(X)
    X = Dropout(0.2)(X)
    # X = Dense(300, activation='relu', kernel_initializer='normal', kernel_regularizer=l2(0.01))(X)
    # X = Dropout(0.2)(X)
    X = Dense(1, kernel_initializer='normal', activation='linear', kernel_regularizer=l2(0.01))(X)
    
    model = Model(inputs=X_input, outputs=X, name='baselineModel')
    return model


# create a function that will return the X and Y objects for a given muscle and datatype
def importGRFs(grffilelist):
    # create a templist
    templist = []
    for file in grffilelist:
        # load the file in to a dataframe with its info
        tempsubj = file[0:7]
        tempcond = file[8:22]
        temptrial = file[23:-7]
        tempexp = file[0:4]
        tempfile = os.path.join(grffilespath,file)
        tempdata = pd.read_csv(tempfile, skiprows=6, delimiter='\t')
        # select the data we want
        temptime = tempdata['time'].to_numpy().reshape(-1,1)
        grfx = tempdata['ground_force_r_vx'].to_numpy().reshape(-1,1)
        grfy = tempdata['ground_force_r_vy'].to_numpy().reshape(-1,1)
        grfz = tempdata['ground_force_r_vz'].to_numpy().reshape(-1,1)

        ### get all the different datatypes
        # resample curves for a set size
        tempgrfx = np.trim_zeros(grfx)
        tempgrfy = np.trim_zeros(grfy)
        tempgrfz = np.trim_zeros(grfz)

        testx = np.argwhere(tempgrfx == 0)
        testy = np.argwhere(tempgrfy == 0)
        testz = np.argwhere(tempgrfz == 0)
        if testx.shape[0] != 0:
            tempgrfx = tempgrfx[0:testx[0,0],:]
        if testy.shape[0] != 0:
            tempgrfy = tempgrfy[0:testy[0,0],:]
        if testz.shape[0] != 0:
            tempgrfz = tempgrfz[0:testz[0,0],:]

        # fig = plt.figure()
        # plt.plot(tempgrfy)

        # fig, axs = plt.subplots(3, sharex=True, sharey=False)
        # fig.suptitle('GRF Values', fontsize=16)
        # axs[0].plot(tempgrfx, color='red', lw=3)
        # # axs[0].set_xticks(fontsize=16)
        # # axs[0].set_yticks(fontsize=16)
        # axs[0].set(ylabel='Force [N]') #, fontsize=16)

        # axs[1].plot(tempgrfy, color='darkgreen', lw=3)
        # axs[1].set(ylabel='Force[N]') #, fontsize=16)
        # # axs[1].set_yticks(fontsize=16)
        # # axs[1].set_ylabel('Force [N]', fontsize=16)

        # axs[2].plot(tempgrfz, color='blue', lw=3)
        # # plt.tight_layout()
        # # axs[2].set_xticks(fontsize=16)
        # # axs[2].set_yticks(fontsize=16)
        # axs[2].set(ylabel='Force [N]') #, fontsize=16)

        # plt.show()
        # import time
        # time.sleep(30)
        newgrfx = scipy.signal.resample(x=tempgrfx, num=25) #, t=temptime) newtime
        newgrfy = scipy.signal.resample(tempgrfy, 25)
        newgrfz = scipy.signal.resample(tempgrfz, 25)
        newgrf = np.vstack((newgrfx, newgrfy, newgrfz))
        
        # binned samples
        newgrf = scipy.signal.resample(newgrf, (newgrf.shape[0]//30)*30)
        test = newgrf.reshape(30, (newgrf.shape[0] // 30))
        test = test.mean(axis=1)
        test = test.reshape(-1,1)
        
        # grab peak values
        peakx = max(grfx)
        peaky = max(grfy)
        peakz = max(grfz)
        # grab integral values
        tempintx = np.abs(grfx)
        tempinty = np.abs(grfy)
        tempintz = np.abs(grfz)
        integx = np.trapz(x=temptime, y=tempintx, axis=0)
        integy = np.trapz(x=temptime, y=tempinty, axis=0)
        integz = np.trapz(x=temptime, y=tempintz, axis=0)
        # do one with a combo of the integral and peaks
        # peakintegx = np.array([peakx,integx])
        peakinteg = np.array([peakx, peaky, peakz, integx, integy, integz])
        # peakintegy = np.array([peaky,integy])
        # peakintegz = np.array([peakz,integz])

        ## make a temp dataframe
        tempmuscle_df = pd.DataFrame({'subjectname':[tempsubj], 'condname':[tempcond],
                                    'trialname':[temptrial], 'experimentname':[tempexp],
                                    'grf x curve':[newgrfx], 'grf y curve':[newgrfy], 'grf z curve':[newgrfz],
                                    'grf xyz peaks':[np.array([peakx, peaky, peakz])],
                                    'grf xyz integral':[np.array([integx, integy, integz])],
                                    'grf xyz peak and integral':[peakinteg],
                                    'grf xyz curves':[newgrf],
                                    'grf xyz binned':[test]})
        templist.append(tempmuscle_df)

    
    # combine each of the dataframes into one
    grf_df = pd.concat(templist, ignore_index=True)
    grf_df.sort_values(by=['subjectname', 'condname','trialname'], inplace=True)
    grf_df.drop([40], inplace=True)
    print('\nMake sure that you actually want to drop this one entry from subject 14!\n')
    grf_df.reset_index(inplace=True, drop=True)

    return grf_df

################################################################



################################################################
# ## import the experimental values
# set paths
repobasedir = os.getcwd()
experimentalfile = os.path.join(repobasedir, 'experimentalMetabolics_all.csv')
# read in the file to dataframe
expmetcost_df = pd.read_csv(experimentalfile)
# print(expmetcost_df)

################################################################
## import the simulation muscle metabolic values
# set all the paths
simresultspath = os.path.join(repobasedir,'..\\metabolicsResults\\')
muscleInversepath = os.path.join(simresultspath,'muscleInverse\\')
muscleInverseWithEMGpath = os.path.join(simresultspath,'muscleInverseWithEMG\\')

## first handle the values in the regular muscle driven inverse problem
# get all the filenames
musclefiles = glob.glob(os.path.join(muscleInversepath,'*.csv'))
# load them all into a single dataframe
df_from_each_file = (pd.read_csv(f) for f in musclefiles)
# print(df_from_each_file)
muscle_df = pd.concat(df_from_each_file, ignore_index=True)
muscle_df.sort_values(by=['subjectname','condname','trialname'],inplace=True)
muscle_df.drop(['Row'],inplace=True,axis=1)
print('\ncheck if you want to drop this one value\n')
muscle_df.drop([40], inplace=True)
muscle_df.reset_index(inplace=True,drop=True)
# print(muscle_df)


################################################################
## import all the GRFs  -  X
# set all the paths
grffilespath = os.path.join(repobasedir,'..\\expfiles\\grffiles\\')
grffilelist = os.listdir(grffilespath)
# print(grffilelist)

muscles = ['metabolics_sol_avg'] # ['metabolics_bifemlh_avg', 'metabolics_recfem_avg','metabolics_sol_avg','metabolics_gas_avg'] 
grf_datatype = ['grf xyz binned'] #, 'grf x curve', 'grf y curve', 'grf z curve', 'grf xyz curves', 'grf xyz peaks'] #, 'grf xyz integral', 'grf xyz peak and integral', 'grf xyz binned']
modeltypes = ['Random Forest','Neural Network'] # 'Linear Regression'

# get the grf dataframe
grf_df = importGRFs(grffilelist)


# now have two dataframes, one with all the grf data in different forms
# and one with the muscle metabolics data with one value 
# for each muscle for each gait cycle
# grf_df -> X
# muscle_df -> Y


# now for each muscle 
# for each model type
# for each input type
# do the cross validation for each subject
full_analysis = {}

for muscle in muscles:
    print('\nStarting muscle analysis: %s' % muscle)
    for model in modeltypes:
        print('\n    Starting model analysis: %s' % model)
        for datatype in grf_datatype:
            print('\n        Starting input data type: %s' % datatype)
            # get the X
            X_df =  np.stack(grf_df[datatype])
            X_df = X_df[:,:,0]


            # get the model masses
            model_masses = muscle_df['model_mass'].values.reshape(-1,1)
            X_df = np.insert(X_df, [0], model_masses, axis=1)
            print('            X shape: ' + str(X_df.shape))
            # get the Y
            Y_df = muscle_df[muscle].values.reshape(-1,1)
            print('            Y shape: ' + str(Y_df.shape))


            # now to cycle through test subjects and create models and save the errors
            rmse_list_train = []
            mae_list_train = []
            rmse_list_test = []
            mae_list_test = []

            r2_list_train = []
            r2_list_test = []

            isubject = 0
            while isubject <= X_df.shape[0]:
                # print(isubject)


                if isubject == 0:
                    x_test = X_df[isubject:isubject+6,:]
                    y_test = Y_df[isubject:isubject+6]
                    x_train = X_df[isubject+6:-1,:]
                    y_train = Y_df[isubject+6:-1]
                
                elif X_df.shape[0] - isubject <= 6:
                    # we are at the end
                    x_test = X_df[isubject:-1,:]
                    y_test = Y_df[isubject:-1]
                    x_train = X_df[0:isubject,:]
                    y_train = Y_df[0:isubject,:]
                else:
                    x_test = X_df[isubject:isubject+6,:]
                    y_test = Y_df[isubject:isubject+6]
                    x_train = X_df[0:isubject, :]
                    tempxtrain = X_df[isubject+6:-1,:]
                    x_train = np.append(x_train, tempxtrain, axis=0)
                    y_train = Y_df[0:isubject,:]
                    tempytrain = Y_df[isubject+6:-1]
                    y_train = np.append(y_train, tempytrain, axis=0)

                # have the training and test sets now do stuff
                if model == 'Linear Regression':
                    lin_reg = LinearRegression()
                    lin_reg.fit(x_train, y_train)

                    # make predictions
                    y_pred_train = lin_reg.predict(x_train)
                    y_pred_test = lin_reg.predict(x_test)

                    # rmse + mae - train
                    mse_train = mean_squared_error(y_true=y_train, y_pred=y_pred_train)
                    rmse_train = np.sqrt(mse_train)
                    rmse_list_train.append(rmse_train)
                    mae_train = mean_absolute_error(y_true=y_train, y_pred=y_pred_train)
                    mae_list_train.append(mae_train)
                    # rmse + mae - train
                    mse_test = mean_squared_error(y_true=y_test, y_pred=y_pred_test)
                    rmse_test = np.sqrt(mse_test)
                    rmse_list_test.append(rmse_test)
                    mae_test = mean_absolute_error(y_true=y_test, y_pred=y_pred_test)
                    mae_list_test.append(mae_test)
                    # r2 score for linear reg
                    r2_list_train.append(r2_score(y_true=y_train, y_pred=y_pred_train))
                    r2_list_test.append(r2_score(y_true=y_test, y_pred=y_pred_test))

                if model == 'Neural Network':
                    basenn = baseline_model(x_train.shape[1:])
                    basenn.compile(optimizer='Adam', loss='mse', 
                        metrics=['mean_absolute_error','mean_squared_error', 'mape'])
                    basenn.fit(x=x_train, y=y_train, epochs=3000, verbose=0, 
                        shuffle=True, batch_size=None)
                    
                    preds = basenn.evaluate(x=x_test, y=y_test, verbose=0)
                    # print('#######################################')
                    # print('Loss = %f' % preds[0])
                    # print('RMSE = %f' % np.sqrt(preds[2]))
                    # print('MAE = %f' % preds[1])
                    # print('MAPE= %f ' % preds[3])

                    y_pred_train = basenn.predict(x=x_train)
                    y_pred_test = basenn.predict(x=x_test)

                    # rmse + mae - train
                    mse_train = mean_squared_error(y_true=y_train, y_pred=y_pred_train)
                    rmse_train = np.sqrt(mse_train)
                    rmse_list_train.append(rmse_train)
                    mae_train = mean_absolute_error(y_true=y_train, y_pred=y_pred_train)
                    mae_list_train.append(mae_train)
                    # rmse + mae - train
                    mse_test = mean_squared_error(y_true=y_test, y_pred=y_pred_test)
                    rmse_test = np.sqrt(mse_test)
                    rmse_list_test.append(rmse_test)
                    mae_test = mean_absolute_error(y_true=y_test, y_pred=y_pred_test)
                    mae_list_test.append(mae_test)
                    # r2 score for linear reg
                    r2_list_train.append(r2_score(y_true=y_train, y_pred=y_pred_train))
                    r2_list_test.append(r2_score(y_true=y_test, y_pred=y_pred_test))

                
                else:
                    ran_for = RandomForestRegressor(random_state=13, n_estimators=1000)
                    ran_for.fit(x_train, y_train[:,0])

                    # make predictions
                    y_pred_train = ran_for.predict(x_train)
                    y_pred_test = ran_for.predict(x_test)

                    # rmse + mae - train
                    mse_train = mean_squared_error(y_true=y_train, y_pred=y_pred_train)
                    rmse_train = np.sqrt(mse_train)
                    rmse_list_train.append(rmse_train)
                    mae_train = mean_absolute_error(y_true=y_train, y_pred=y_pred_train)
                    mae_list_train.append(mae_train)
                    # rmse + mae - train
                    mse_test = mean_squared_error(y_true=y_test, y_pred=y_pred_test)
                    rmse_test = np.sqrt(mse_test)
                    rmse_list_test.append(rmse_test)
                    mae_test = mean_absolute_error(y_true=y_test, y_pred=y_pred_test)
                    mae_list_test.append(mae_test)
                    # r2_score
                    r2_list_train.append(r2_score(y_true=y_train, y_pred=y_pred_train))
                    r2_list_test.append(r2_score(y_true=y_test, y_pred=y_pred_test))
                # cycle through subjects
                isubject+=6


            print('            Average RMSE training: %f W' % np.mean(rmse_list_train))
            print('            Average RMSE testing: %f W' % np.mean(rmse_list_test))
            print('            Average MAE training: %f W' % np.mean(mae_list_train))
            print('            Average MAE testing: %f W' % np.mean(mae_list_test))
            print('            Average R^2 training: %f W' % np.mean(r2_list_train))
            print('            Average R^2 testing: %f W' % np.mean(r2_list_test))

            # have scores for each test set on this model 
            full_analysis[muscle + model + datatype + ' rmse_train'] = rmse_list_train
            full_analysis[muscle + model + datatype + ' rmse_test'] = rmse_list_test
            full_analysis[muscle + model + datatype + ' mae_train'] = mae_list_train
            full_analysis[muscle + model + datatype + ' mae_test'] = mae_list_test
            full_analysis[muscle + model + datatype + ' r2_train'] = r2_list_train
            full_analysis[muscle + model + datatype + ' r2_test'] = r2_list_test


print('\nDone?')
