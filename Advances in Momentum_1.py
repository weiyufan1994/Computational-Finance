
#%%
import numpy as np
import pandas as pd
import datetime as dt
from scipy.stats import skew
import matplotlib.pyplot as plt


#%% [markdown]
## Load Features Data

#%%
features = pd.read_csv("CRSP_features.csv")

# ## Replace not valid with NaN
# features.replace(["B", "C"], np.nan, inplace = True)

# #%%
# # convert percent string to float
# features["RET"] = features["RET"].str.rstrip('%').astype('float') / 100.0 
# features["ivol"] = features["ivol"].str.rstrip('%').astype('float') / 100.0 
# features["tvol"] = features["tvol"].str.rstrip('%').astype('float') / 100.0 
# features["R2"] = features["R2"].str.rstrip('%').astype('float') / 100.0 
# features["exret"] = features["exret"].str.rstrip('%').astype('float') / 100.0 

# # Convert "date" to datetime
# features["DATE"] = pd.to_datetime(features["DATE"], format = "%Y%m%d")
# features["DATE"] = features["DATE"] + pd.offsets.MonthEnd(0)
# # Add explicit year and month as new columns
# features["Year"] = features["DATE"].dt.year
# features["Month"] = features["DATE"].dt.month




######### Load Volatility Data ###########
#features_tech = pd.read_csv("CRSP_features_tech.csv")


#%%

# Create training set for winning decile
# Load Momentum file 
def PS3_Q1(dat):
    #### STEP One: define the universe of monthly returns and ranking returns ####
    ##############################################################################

    # Look at the data types
    dat.dtypes

    # Fill na with -99, then transform the datatype into int
    dat[["SHRCD", "EXCHCD", "SHROUT"]] = dat[["SHRCD", "EXCHCD", "SHROUT"]].fillna(-99)

    dat[["PERMNO", "SHRCD", "EXCHCD", "SHROUT"]] = dat[["PERMNO", "SHRCD", "EXCHCD", "SHROUT"]].astype(int)

    ## Replace -99 with NaN
    dat.replace([-99], np.nan, inplace = True)

    # Convert "date" to datetime
    dat["date"] = pd.to_datetime(dat["date"], format = "%Y%m%d")

    # Push the "last trading date of a month 7-29" to the "last calender date of a month 7-31" 
    dat["date"] = dat["date"] + pd.offsets.MonthEnd(0)

    # Add explicit year and month as new columns
    dat.loc[:,"Year"] = dat["date"].dt.year
    dat.loc[:,"Month"] = dat["date"].dt.month


    ################# Slicing the data and Do computation #########
    # subset the data by Share Code and Exchange Code, time
    dat = dat.query("SHRCD == [10, 11] and EXCHCD == [1, 2, 3]")

    # calculate market cap t
    dat["mktval"] = np.abs(dat["PRC"])*dat["SHROUT"]

    # Groupby each stock, take lag on its market value
    dat["lag_Mkt_Cap"] = dat.groupby(["PERMNO"])["mktval"].shift(1)

    # Convert to million
    dat["lag_Mkt_Cap"] = dat["lag_Mkt_Cap"]/1000

    ## Replace characters with NaN
    dat.replace(["A", "S", "C", "T", "P"], np.nan, inplace = True)

    # convert datatype to float
    dat["RET"] = dat["RET"].astype(float)

    ##### orignal method
    # ## calculate real return: Aggregate delisting return and HPR
    # Ret = np.where(dat["DLRET"].notna(), dat["DLRET"].astype(float), 
    #                                          dat["RET"])

    ##### ignore delisting return
    Ret = dat["RET"]

    Ret = np.where(dat["RET"].notna() & dat["DLRET"].notna(), 
                (dat["DLRET"].astype(float) + 1)*(dat["RET"] + 1) - 1, Ret)

    # create new real return
    dat["Ret"] = Ret


    ################# Create Ranking Return #########
    # create a subset to work on with
    dat_temp = dat[['PERMNO','date','RET']].sort_values(['PERMNO','date'])\
        .set_index('date')

    # create log return at period s
    dat_temp["logret"] = np.log(1 + dat_temp.RET)

    # shift return by 2 
    dat_temp["shift_log_ret"] = dat_temp.groupby(["PERMNO"])["logret"].shift(2)

    # cumulative log return
    cum_ret = dat_temp.groupby(["PERMNO"])["shift_log_ret"].rolling(11).sum()

    # reset index and rename column, merge ranking return
    cum_ret = cum_ret.reset_index().rename(columns = {"shift_log_ret" : "Ranking_Ret"})
    dat = pd.merge(dat, cum_ret, on = ["PERMNO", "date"])



    ################# Apply Data Restriction #################
    def data_restriction(df):

        tmp_ret = df.Ranking_Ret

        # initiate an indicator vector, length equal to the rows of input df
        valid_indicator = np.ones(df.shape[0])

        # check price t-13 is not missing 
        valid_indicator[np.isnan(df.PRC.shift(13))] = np.nan

        # check ret(t-2) is not missing 
        valid_indicator[np.isnan(df.RET.shift(2))] = np.nan    

        # check me(t-1) is not missing
        valid_indicator[np.isnan(df.lag_Mkt_Cap)] = np.nan

        # check a minimum of eight monthly return over past 11 months
        valid_indicator[np.isnan(df.Ret.shift(1)).rolling(11).sum() > 3] = np.nan

        # update valid return
        valid_return = valid_indicator*tmp_ret

        return valid_return

    out2 = dat.groupby(["PERMNO"]).apply(data_restriction)

    # updata the Ranking_Ret with restricted ranking return
    dat.Ranking_Ret = out2.reset_index(0, drop=True)


    # subset the output dataframe
    out_df = dat[["Year", "Month", "PERMNO", "EXCHCD", "lag_Mkt_Cap", "Ret", "Ranking_Ret"]]
    out_df = out_df.sort_values(['Year', 'Month'])
    out_df = out_df[out_df.Year >= 1927]
    
    return out_df

def PS3_Q2(dat):

    #### STEP Two: define the monthly momentum portfolio decile ####
    ################################################################

    #### DM CUT - based on all stocks ####
    # cut into decile
    dat["DM_decile"] = dat.groupby(["Year", "Month"])["Ranking_Ret"].transform(lambda x: pd.qcut(x, 10, labels=False) + 1)

    #### NYSE CUT - based on nyse stocks ####
    # subset stocks in NYSE
    nyse = dat[dat['EXCHCD']==1]

    # cut nyse by percentiles
    nyse_break = nyse.groupby(["Year", "Month"])["Ranking_Ret"]\
                    .describe(percentiles=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])\
                    .reset_index()

    # get nyse break points
    nyse_break = nyse_break[["Year", "Month", "10%", "20%", "30%", "40%", \
                                   "50%", "60%", "70%", "80%", "90%"]]

    def nyse_cut(df, nyse_break):

        year = df.Year.min()
        month = df.Month.min()

        # EACH POINT OF TIME get nyse break points
        bk = nyse_break[(nyse_break.Year == year) & (nyse_break.Month == month)]
        bk = np.array(bk).flatten()[2:]
        bk = np.insert(np.append(bk, np.inf), 0, -np.inf)

        # cut all firms into nyse bins
        rank_labels = pd.cut(df.Ranking_Ret, bk, labels = False) + 1

        return rank_labels

    out = dat.groupby(["Year", "Month"]).apply(nyse_cut, nyse_break = (nyse_break))

    # record KRF_decile to new column 
    dat["KRF_decile"] = out.reset_index([0,1]).Ranking_Ret

    return dat

################ run functions ##################################
pd.set_option('mode.chained_assignment', None)

# Read Stock Data
dat = pd.read_csv("CRSP_Stocks.csv")

CRSP_Stocks_Momentum = PS3_Q1(dat)

CRSP_Stocks_Momentum_decile = PS3_Q2(CRSP_Stocks_Momentum.copy())


#%%

# Create labels for the winner decile
CRSP_Stocks_Momentum_decile["lag_ret"] = CRSP_Stocks_Momentum_decile.groupby(["PERMNO"])["Ret"].shift(1)
winning_decile = CRSP_Stocks_Momentum_decile[CRSP_Stocks_Momentum_decile.KRF_decile == 10 ]
winning_decile["Performance"] = np.where(winning_decile.Ret > 0, 1, 0)
winning_decile = winning_decile[['Performance', 'Year', 'Month', 'PERMNO', 'Ret', 'lag_ret', 'Ranking_Ret', 'lag_Mkt_Cap']]

# Create labels for the loser decile
Loser_decile = CRSP_Stocks_Momentum_decile[CRSP_Stocks_Momentum_decile.KRF_decile == 1 ]
Loser_decile["Performance"] = np.where(Loser_decile.Ret < 0, 1, 0)
Loser_decile = Loser_decile[['Performance', 'Year', 'Month', 'PERMNO', 'Ret', 'lag_ret', 'Ranking_Ret', 'lag_Mkt_Cap']]

# Join Labels with Features
winning_decile_data = pd.merge(winning_decile, features, how="left", on=["PERMNO", "Year", "Month"])
Loser_decile_data = pd.merge(Loser_decile, features, how="left", on=["PERMNO", "Year", "Month"])

# Add more stuff here
winning_decile_data["vol"] = np.where(winning_decile_data.vol == 0, np.nan, winning_decile_data.vol)
winning_decile_data["vol_spx"] = np.where(winning_decile_data.vol_spx == 0, np.nan, winning_decile_data.vol_spx)
winning_decile_data["norm_ret"] = winning_decile_data["Ret"]/ (winning_decile_data["vol"])
winning_decile_data["vol2vol_spx"] = winning_decile_data["vol"]/  (winning_decile_data["vol_spx"])
winning_decile_data["zscore_vol_spread"] = winning_decile_data["zscore_vol"] - winning_decile_data["zscore_vol_spx"]


Loser_decile_data["vol"] = np.where(Loser_decile_data.vol == 0, np.nan, Loser_decile_data.vol)
Loser_decile_data["vol_spx"] = np.where(Loser_decile_data.vol_spx == 0, np.nan, Loser_decile_data.vol_spx)
Loser_decile_data["norm_ret"] = Loser_decile_data["Ret"]/ Loser_decile_data["vol"]
Loser_decile_data["vol2vol_spx"] = Loser_decile_data["vol"]/ Loser_decile_data["vol_spx"]
Loser_decile_data["zscore_vol_spread"] = Loser_decile_data["zscore_vol"] - Loser_decile_data["zscore_vol_spx"]

# Take lag on these columns
lag_columns = winning_decile_data.drop(columns = ['Performance', 'Year', 'Month', 'PERMNO', 'Ret', 'Ranking_Ret', 'lag_Mkt_Cap', 'lag_ret']).columns
winning_decile_data[lag_columns] = winning_decile_data.groupby("PERMNO")[lag_columns].shift(1)
Loser_decile_data[lag_columns] = Loser_decile_data.groupby("PERMNO")[lag_columns].shift(1)


# Sort Value Here
winning_decile_data = winning_decile_data.sort_values(["Year", "Month"])
Loser_decile_data = Loser_decile_data.sort_values(["Year", "Month"])


# Split Data into Training and Testing
training_w = winning_decile_data[(winning_decile_data.Year >= 1973) & (winning_decile_data.Year < 2010)]
testing_w = winning_decile_data[(winning_decile_data.Year >= 2010)]

training_l = Loser_decile_data[(Loser_decile_data.Year >= 1973) & (Loser_decile_data.Year < 2010)]
testing_l = Loser_decile_data[(Loser_decile_data.Year >= 2010)]

feature_labels = ["lag_ret", "Ranking_Ret"] + list(lag_columns)


#%%
print(training_w.shape)
print(testing_w.shape)

# Check for NA values
training_w.dropna(inplace = True)
testing_w.dropna(inplace = True)
training_l.dropna(inplace = True)
testing_l.dropna(inplace = True)

# keep year and month in testing set 
testing_info_w = testing_w[["Year", "Month", "lag_Mkt_Cap", 'Ret']]
testing_info_l = testing_l[["Year", "Month", "lag_Mkt_Cap", 'Ret']]

# drop irrelatvent information
training_w.drop(columns = ["PERMNO", "Year", "Month", "lag_Mkt_Cap", 'Ret'], inplace = True)
testing_w.drop(columns = ["PERMNO", "Year", "Month", "lag_Mkt_Cap", 'Ret'], inplace = True)
training_l.drop(columns = ["PERMNO", "Year", "Month", "lag_Mkt_Cap", 'Ret'], inplace = True)
testing_l.drop(columns = ["PERMNO", "Year", "Month", "lag_Mkt_Cap", 'Ret'], inplace = True)

print(training_w.shape)
print(testing_w.shape)

#%%

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import datasets, metrics
from sklearn.model_selection import TimeSeriesSplit
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import GridSearchCV
import time
#%%

x_train_w = np.array(training_w.iloc[:,1:])
y_train_w = np.array(training_w.iloc[:,0])

x_test_w = np.array(testing_w.iloc[:,1:])
y_test_w = np.array(testing_w.iloc[:,0])

x_train_l = np.array(training_l.iloc[:,1:])
y_train_l = np.array(training_l.iloc[:,0])

x_test_l = np.array(testing_l.iloc[:,1:])
y_test_l = np.array(testing_l.iloc[:,0])

#%%

def Run_Models(x_train, y_train, x_test, y_test):
    
    Logistic = Pipeline([('scl',StandardScaler()),
                     ('clf',LogisticRegression(solver='newton-cg'))])

    RandomForest = RandomForestClassifier(n_estimators=100, n_jobs = -1)

    AdaBoost = Pipeline([('scl',StandardScaler()),
                    ('clf', AdaBoostClassifier(n_estimators=100))])


    summary = pd.DataFrame(columns = ['Name','OOS Accuracy','Time'])

    clf_labels = ['Logistic','RandomForest','AdaBoost']

    clf_list = [Logistic, RandomForest, AdaBoost]


    
    for i,clf,clf_label in zip(range(len(clf_list)), clf_list, clf_labels):
        
        start_time = time.time()

        clf.fit(x_train, y_train)
        y_predict = clf.predict(x_test)
        accuracy = metrics.accuracy_score(y_test, y_predict)
        clf_time = time.time() - start_time

        summary.loc[i] = [clf_label, accuracy, clf_time]

        
        
    return summary, clf_list

summary_w , clf_list_w = Run_Models(x_train_w, y_train_w, x_test_w, y_test_w)
summary_l , clf_list_l = Run_Models(x_train_l, y_train_l, x_test_l, y_test_l)

pred_probs_w = [clf.predict_proba(x_test_w)[:,1] for clf in clf_list_w]
pred_probs_l = [clf.predict_proba(x_test_l)[:,1] for clf in clf_list_l]



#%%

clf_label = ["LG", "RF", "ADA"]
Prob_label = ["Prob" + label for label in clf_label]
Scaled_Prob_label = ["ScaledProb" + label for label in clf_label]

for label, prob in zip(Prob_label, pred_probs_w):

    testing_w[label] = prob

testing_w["Year"] = testing_info_w["Year"]
testing_w["Month"] = testing_info_w["Month"]
testing_w["lag_Mkt_Cap"] = testing_info_w["lag_Mkt_Cap"]
testing_w["Ret"] = testing_info_w["Ret"]


for label, prob in zip(Prob_label, pred_probs_l):
    
    testing_l[label] = prob

testing_l["Year"] = testing_info_l["Year"]
testing_l["Month"] = testing_info_l["Month"]
testing_l["lag_Mkt_Cap"] = testing_info_l["lag_Mkt_Cap"]
testing_l["Ret"] = testing_info_l["Ret"]

#%%
def min_max_scale(df):
    return (df - df.min(axis = 0))/(df.max(axis = 0) - df.min(axis = 0))

testing_w[Scaled_Prob_label] = testing_w.groupby(["Year", "Month"])[Prob_label].apply(min_max_scale)
testing_l[Scaled_Prob_label] = testing_l.groupby(["Year", "Month"])[Prob_label].apply(min_max_scale)

def get_prob_weights(df):
    return df/np.sum(df, axis = 0)

testing_w[Scaled_Prob_label] = testing_w.groupby(["Year", "Month"])[Scaled_Prob_label].apply(get_prob_weights)
testing_l[Scaled_Prob_label] = testing_l.groupby(["Year", "Month"])[Scaled_Prob_label].apply(get_prob_weights)


# Ensemble Prob
testing_w["EnsembleProb"] = testing_w["ScaledProbRF"] + testing_w["ScaledProbADA"]
testing_l["EnsembleProb"] = testing_l["ScaledProbRF"] + testing_l["ScaledProbADA"]

testing_w["EnsembleProb"] = testing_w.groupby(["Year", "Month"])["EnsembleProb"].apply(get_prob_weights)
testing_l["EnsembleProb"] = testing_l.groupby(["Year", "Month"])["EnsembleProb"].apply(get_prob_weights)


#%%

def decile_weighted_return(df):

    lag_MV = df.lag_Mkt_Cap.sum()

    weights = df.lag_Mkt_Cap/lag_MV

    weighted_ret = weights*df.Ret

    return weighted_ret.sum()


def prob_weighted_return(df):
    weights = df.EnsembleProb

    weighted_ret = weights*df.Ret

    return weighted_ret.sum()

# calculate return for winner decile 
out = testing_w.groupby(["Year", "Month"]).apply(decile_weighted_return)
out = out.reset_index([0,1]).rename(columns = {0 : "Winner_Ret_vwret"})

out2 = testing_w.groupby(["Year", "Month"]).apply(prob_weighted_return)
out2 = out2.reset_index([0,1]).rename(columns = {0 : "Winner_Ret_pwret"})

out3 = testing_l.groupby(["Year", "Month"]).apply(decile_weighted_return)
out3 = out3.reset_index([0,1]).rename(columns = {0 : "Loser_Ret_vwret"})

out4 = testing_l.groupby(["Year", "Month"]).apply(prob_weighted_return)
out4 = out4.reset_index([0,1]).rename(columns = {0 : "Loser_Ret_pwret"})


out = pd.merge(out, out2, how = "inner", on=["Year", "Month"])
out = pd.merge(out, out3, how = "inner", on=["Year", "Month"])
out = pd.merge(out, out4, how = "inner", on=["Year", "Month"])

out["WML_vwret"] = out["Winner_Ret_vwret"] - out["Loser_Ret_vwret"]
out["WML_pwret"] = out["Winner_Ret_pwret"] - out["Loser_Ret_pwret"]

#%%

import statsmodels.formula.api as smf
result = smf.ols('WML_pwret ~ WML_vwret', data= out).fit()
result.summary()


#%%
out["date"] = out["Year"].astype(str) + out["Month"].astype(str)

out["date"] = pd.to_datetime(out["date"], format = "%Y%m")

plt.plot(out["date"], np.cumsum(np.log(1 + out.WML_vwret)), color = "red")
plt.plot(out["date"], np.cumsum(np.log(1 + out.WML_pwret)), color = "blue")

#%%
clf_RandomForest = clf_list_w[1]
importances = clf_RandomForest.feature_importances_

std = np.std([tree.feature_importances_ for tree in clf_RandomForest.estimators_],axis=0)

indices = np.argsort(-importances)[::-1]

features_text = [feature_labels[i] for i in indices]

# Print the feature ranking
print("Feature ranking:")

for f in range(x_train_w.shape[1]):
    print("%d. feature %s (%f)" % (f + 1, features_text[f], importances[indices[f]]))

# Plot the feature importances of the forest
plt.figure(figsize=(15,12))
plt.title("Feature importances")
ax = plt.barh(y = np.arange(x_train_w.shape[1]), width = importances[indices], alpha=.8,
       color="#6495ed", xerr=std[indices], error_kw={'ecolor':'r'})



for i,j in zip(ax.patches, features_text):
    plt.text(i.get_width()/2, i.get_y() + i.get_height()/2, j, fontsize=15, ha='center',va='center', color='black')


#%%
