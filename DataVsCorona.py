# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 13:51:40 2020

@author: Utkarsh Panara
"""
import pandas as pd
from torch.utils.data import Dataset
import torch

path = "C:/Users/Utkarsh Panara/Desktop/DataVsCorona/Local_Covid19_master/data/"

time_series_covid19_confirmed_US=pd.read_csv(path+"time_series_covid19_confirmed_US.csv", encoding='latin1')
print('time_series_covid19_confirmed_US', time_series_covid19_confirmed_US.shape)

time_series_covid19_deaths_US=pd.read_csv(path+"time_series_covid19_deaths_US.csv", encoding='latin1')
print('time_series_covid19_deaths_US', time_series_covid19_deaths_US.shape)

time_series_covid19_confirmed_US=time_series_covid19_confirmed_US.drop(['UID','iso2','iso3',
                                                                        'code3','Admin2',
                                                                        'Province_State',
                                                                        'Country_Region','Lat',
                                                                        'Long_','Combined_Key'],axis=1)

confirmed_cases = []
for i in range(time_series_covid19_confirmed_US.shape[0]):
    trial=time_series_covid19_confirmed_US.iloc[i:i+1, 0:]
    trial = trial.transpose()
    trial.columns = ['confirmed_cases']
    fips=trial.iloc[0]['confirmed_cases']
    trial = trial.drop('FIPS', axis=0)
    trial['FIPS'] = pd.Series([fips for x in range(len(trial.index))], index=trial.index)
    confirmed_cases.append(trial)
    
final_confirmed_cases=pd.concat(confirmed_cases, sort=False, ignore_index=False)


time_series_covid19_deaths_US=time_series_covid19_deaths_US.drop(['UID','iso2','iso3',
                                                                        'code3','Admin2',
                                                                        'Province_State',
                                                                        'Country_Region','Lat',
                                                                        'Long_','Combined_Key','Population'],axis=1)

deaths = []
for i in range(time_series_covid19_deaths_US.shape[0]):
    trial=time_series_covid19_deaths_US.iloc[i:i+1, 0:]
    trial = trial.transpose()
    trial.columns = ['deaths']
    fips=trial.iloc[0]['deaths']
    trial = trial.drop('FIPS', axis=0)
    trial['FIPS'] = pd.Series([fips for x in range(len(trial.index))], index=trial.index)
    deaths.append(trial)
    
final_deaths=pd.concat(deaths, sort=False, ignore_index=False)

combined_df = pd.concat([final_confirmed_cases, final_deaths], axis=1)

#Split Data
training_set = combined_df.iloc[0:192623, 0:]
print('Size of training data is :' , training_set.shape )

test_set = combined_df.iloc[192623: , 0:]
print('Size of test data is :' , test_set.shape )

window_size=30
#Dataloader
    
class train_dataset(Dataset):
    def __init__(self, training_set):
        self.input=training_set.iloc[0:, 0:2]
        self.target=training_set.iloc[0:,2:]
        
        # self.input_scaled = scaler_train.fit_transform(self.input)
        
        self.input_groupped = self.input.groupby('FIPS', sort = False)
        
        self.input_train = []
        for FIPS, confirmed_cases in self.input_groupped:
            # print(FIPS,confirmed_cases.shape)
            confirmed_cases = confirmed_cases.T
            for p in range(window_size, confirmed_cases.shape[1]+1):  
                xtrain = confirmed_cases.iloc[0:1, p-window_size:p]
                xtrain = xtrain.to_numpy()
                self.input_train.append(xtrain)
          
        self.input_train = torch.FloatTensor(self.input_train)
        print('The shape of self.input_train' , self.input_train.shape) 
        
        self.target_groupped = self.target.groupby('FIPS', sort = False)
        
        self.target_train = []
        for FIPS, deaths in self.target_groupped:
            # print(FIPS,deaths.shape)
            for p in range(window_size, deaths.shape[0]+1): 
                ytrain = deaths.iloc[p-1:p, 0:1]
                ytrain = ytrain.to_numpy()
                self.target_train.append(ytrain)
          
        self.target_train = torch.FloatTensor(self.target_train)
        self.target_train = self.target_train.view(-1,1)
        print('The shape of self.target_train' , self.target_train.shape) 

    def __getitem__(self, index):
        return self.input_train[index],  self.target_train[index]
    
    def __len__(self):
        return self.input_train.shape[0]

tr_dataset = train_dataset(training_set)
print('Length of training data' , len(tr_dataset))
        


























