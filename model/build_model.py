import pandas as pd
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.pipeline import FeatureUnion
from sklearn.base import BaseEstimator,TransformerMixin
from sklearn import preprocessing
import pickle
import xgboost as xgb
import time
from sklearn.model_selection import train_test_split as TTS
from sklearn.model_selection import GridSearchCV as GSCV
from pathlib import Path
from sklearn import metrics
from sklearn.metrics import roc_auc_score,recall_score,precision_score
from warnings import simplefilter
import tqdm
import os
import glob
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import StratifiedKFold, KFold
from sklearn.preprocessing import Normalizer
#import shap
import matplotlib.pyplot as plt
import time
import warnings

warnings.filterwarnings('ignore')
simplefilter(action='ignore', category=FutureWarning)
datatime = time.strftime("%Y-%m-%d",time.localtime())

ohencoder = preprocessing.OneHotEncoder(handle_unknown='ignore', sparse=False)
normalscl = preprocessing.Normalizer()

otu = pd.read_excel('../data/OUT1.xlsx',index_col=0)
Q300 = pd.read_excel('../data/fecal_metabolitie_1.xlsx',index_col=3).iloc[:,4:]
Q600 = pd.read_excel('../data/serum_metabolitie_1.xlsx',index_col=3).iloc[:,4:]
mre = pd.read_csv('../data/MRE_1.csv',index_col=0).iloc[:,:].transpose()
cov = pd.read_csv('../data/COV_1.csv',index_col=0).iloc[:,:].transpose()

otu2 = pd.read_excel('../data/OTU_2.xlsx',index_col=0)
Q3002 = pd.read_excel('../data/fecal_metabolitie_2.xlsx',index_col=3)
Q3002 = Q3002.iloc[:, 4:]
Q6002 = pd.read_excel('../data/serum_metabolitie_2.xlsx',index_col=3)
Q6002 = Q6002.iloc[:, 4:]
mre2 = pd.read_excel('../data/MRE_2.xlsx',index_col=0).iloc[:,2:].transpose()
cov2 = pd.read_excel('../data/COV_2.xlsx',index_col=0).iloc[:,2:].transpose()

otu.index = [x.split('g__')[-1] for x in otu.index]
otu2.index = [x.split('g__')[-1] for x in otu2.index]
Q300.index = Q300.index.map('Fecal_{}'.format)
Q600.index = Q600.index.map('Serum_{}'.format)
Q3002.index = Q3002.index.map('Fecal_{}'.format)
Q6002.index = Q6002.index.map('Serum_{}'.format)
Q3002.rename(index={'Fecal_L-Alanine': 'Fecal_Alanine'}, inplace=True)

tmpscl = normalscl
treat = pd.read_csv('../data/treat.csv',index_col=0)
group_df = treat.loc[treat['Group'] == 'CD',:]
cov = cov.apply(lambda row: row.fillna(row.mean()), axis=1)
cov2 = cov2.apply(lambda row: row.fillna(row.mean()), axis=1)

cov_select = ['Gender','Age','Location','BMI']
Q300_select = ['Fecal_Alanine',
               'Fecal_Pimelic acid',
               'Fecal_Suberic acid',
               'Fecal_Arachidonic acid',
               'Fecal_Glucose 6-phosphate',
               'Fecal_1-Methylhistidine']
Q600_select = ['Serum_CE(19:3)', 'Serum_ePE(36:4)', 'Serum_ePS(38:2)']
mre_select = ['ADC', 'Thickness', 'Perianal diseases']
otu_select = ['[Ruminococcus]_gnavus_group', 'Erysipelatoclostridium', 'Saccharimonadaceae']


# 合并DataFrame，使用keys区分不同来源的数据
combined_otu = pd.concat([otu, otu2], axis=1)

# 初始化归一化器
normalizer = Normalizer()

# 对合并后的数据进行归一化处理
combined_normalized_otu = normalizer.fit_transform(combined_otu)

# 将归一化后的numpy数组转换回DataFrame
combined_normalized_otu_df = pd.DataFrame(combined_normalized_otu, index=combined_otu.index, columns=combined_otu.columns)

# 分割归一化后的DataFrame
otu_normalized = combined_normalized_otu_df.iloc[:, :otu.shape[1]]
otu2_normalized = combined_normalized_otu_df.iloc[:, otu.shape[1]:]
otu = otu_normalized

# 找到所有以"QC"开头的列
qc_Q600 = Q600.loc[Q600_select,].filter(regex='^QC')
# 计算这些列的均值
qc_Q600_means = qc_Q600.mean(axis=1)
# 打印质控均值
print("Quality Control Means:")
print(qc_Q600_means)
Q600 = Q600.loc[Q600_select,].div(qc_Q600_means, axis=0)

# 找到所有以"QC"开头的列
qc_Q300 = Q300.loc[Q300_select,].filter(regex='^QC')
# 计算这些列的均值
qc_Q300_means = qc_Q300.mean(axis=1)
# 打印质控均值
print("Quality Control Means:")
print(qc_Q300_means)
Q300 = Q300.loc[Q300_select,].div(qc_Q300_means, axis=0)


save_file_path = './result_multi_algorithm'
os.makedirs(save_file_path, exist_ok=True)
os.makedirs(save_file_path+'/总结', exist_ok=True)

idx_raw1 = (otu.columns & group_df.index) & (Q600.columns & Q300.columns) #获得152个共有的CD
y_ = group_df.loc[idx_raw1,'Group_BD']

y = np.array(y_=='BD2',dtype='int')
x_idx = group_df.loc[idx_raw1,:].index

skfold = StratifiedKFold(n_splits=5,random_state=47,shuffle=True)

k_index = 0
for train_index, test_index in skfold.split(x_idx,y):
    print(y[test_index])
    X_train_cov = cov.loc[cov_select,x_idx[train_index]].transpose()
    X_test_cov = cov.loc[cov_select,x_idx[test_index]].transpose()
    print(X_train_cov.shape)
    X_train_otu = otu.loc[otu_select,x_idx[train_index]].transpose()
    X_test_otu = otu.loc[otu_select,x_idx[test_index]].transpose()
    print(X_train_otu.shape)
    X_train_Q300 = Q300.loc[Q300_select,x_idx[train_index]].transpose()
    X_test_Q300 = Q300.loc[Q300_select,x_idx[test_index]].transpose()
    print(X_train_Q300.shape)
    X_train_Q600 = Q600.loc[Q600_select,x_idx[train_index]].transpose()
    X_test_Q600 = Q600.loc[Q600_select,x_idx[test_index]].transpose()
    print(X_train_Q600.shape)
    X_train_mre = mre.loc[mre_select,x_idx[train_index]].transpose()
    X_test_mre = mre.loc[mre_select,x_idx[test_index]].transpose()
    print(X_train_mre.shape)

    X_train_concat_cov = pd.concat([X_train_otu,X_train_Q300,X_train_Q600,X_train_mre,X_train_cov],axis=1)
    X_test_concat_cov = pd.concat([X_test_otu,X_test_Q300,X_test_Q600,X_test_mre,X_test_cov],axis=1)

    X_train_otu_Q300_Q600_cov = pd.concat([X_train_otu,X_train_Q300,X_train_Q600,X_train_cov],axis=1)
    X_test_otu_Q300_Q600_cov = pd.concat([X_test_otu,X_test_Q300,X_test_Q600,X_test_cov],axis=1)

    X_train_otu_cov = pd.concat([X_train_otu,X_train_cov],axis=1)
    X_test_otu_cov = pd.concat([X_test_otu,X_test_cov],axis=1)

    X_train_Q300_cov = pd.concat([X_train_Q300,X_train_cov],axis=1)
    X_test_Q300_cov = pd.concat([X_test_Q300,X_test_cov],axis=1)

    X_train_Q600_cov = pd.concat([X_train_Q600,X_train_cov],axis=1)
    X_test_Q600_cov = pd.concat([X_test_Q600,X_test_cov],axis=1)

    X_train_mre_cov = pd.concat([X_train_mre,X_train_cov],axis=1)
    X_test_mre_cov = pd.concat([X_test_mre,X_test_cov],axis=1)

    cv = StratifiedKFold(5)
    max_features = None
    alpha = 0.01

    rfclf = RandomForestClassifier(random_state=0)
    dtclf = DecisionTreeClassifier(random_state=0)
    xgbclf = xgb.XGBClassifier(random_state=0)
    lr = LogisticRegression()
    svmclf_linear = SVC(kernel='linear',probability=True)
    svmclf_rbf = SVC(kernel='rbf',probability=True)
    knnclf = KNeighborsClassifier()

    rawmodels = [dtclf,rfclf,xgbclf,lr,svmclf_linear,svmclf_rbf,knnclf]

    dtpara = {'max_depth':range(3,15)}
    rfpara = {'n_estimators':range(10,150,10),'max_depth':range(3,15),'random_state':[0]}
    xgbpara = {'n_estimators':range(10,150,10),'max_depth':range(3,15),'random_state':[0]}
    lrpara = {'C' :np.linspace(0.01,1,100),'penalty' :["l1", "l2"]}
    svmclf_linearpara = { 'C':np.linspace(0.1,10,100)}
    svmclf_rbfpara = { 'C':np.linspace(0.1,10,100)}
    knnpara = {"n_neighbors": np.arange(1, 5, 1),"weights": ["uniform", "distance"],"p": [1, 2]}
    parameters = [dtpara,rfpara,xgbpara,lrpara,svmclf_linearpara,svmclf_rbfpara,knnpara]

    clfnames = ['DT','RFC','XGBC','LR','Linear_SVM','Rbf_SVM','KNN']

    for filename,(X_train,y_train,X_val,y_val) in zip(['cov','otu','Q300','Q600','MRE','otu_Q300_Q600','combined'],
                                                      [(X_train_cov,y[train_index],X_test_cov,y[test_index]),
                                                       (X_train_otu_cov,y[train_index],X_test_otu_cov,y[test_index]),
                                                       (X_train_Q300_cov,y[train_index],X_test_Q300_cov,y[test_index]),
                                                       (X_train_Q600_cov,y[train_index],X_test_Q600_cov,y[test_index]),
                                                       (X_train_mre_cov,y[train_index],X_test_mre_cov,y[test_index]),
                                                       (X_train_otu_Q300_Q600_cov,y[train_index],X_test_otu_Q300_Q600_cov,y[test_index]),
                                                       (X_train_concat_cov,y[train_index],X_test_concat_cov,y[test_index])]):
        os.makedirs(save_file_path+'/'+filename+str(k_index), exist_ok=True)
        X_train.columns = X_train.columns.str.replace('\[', '', regex=True)
        X_val.columns = X_val.columns.str.replace('\[', '', regex=True)
        X_train.columns = X_train.columns.str.replace('\]', '', regex=True)
        X_val.columns = X_val.columns.str.replace('\]', '', regex=True)
        record = pd.DataFrame(columns=['Classifiers'])
        if Path.exists(Path(save_file_path+'/'+filename+str(k_index)+'/'+datatime+'_Record(CV).csv')):
            record = pd.read_csv(save_file_path+'/'+filename+str(k_index)+'/'+datatime+'_Record(CV).csv')
        trainpred_report = pd.DataFrame({'id':X_train.index,'Group':y_train})
        trainprob_report = pd.DataFrame({'id':X_train.index,'Group':y_train})
        testpred_report = pd.DataFrame({'id':X_val.index,'Group':y_val})
        testprob_report = pd.DataFrame({'id':X_val.index,'Group':y_val})

        trainpred_report = pd.DataFrame({'id':X_train.index,'Group':y_train})
        trainprob_report = pd.DataFrame({'id':X_train.index,'Group':y_train})
        testpred_report = pd.DataFrame({'id':X_val.index,'Group':y_val})
        testprob_report = pd.DataFrame({'id':X_val.index,'Group':y_val})
        print('============================',filename,'============================')

        for i,(model,para) in enumerate(zip(rawmodels,parameters)):
            print(str(model))
            clf = GSCV(model,para,cv=cv).fit(X_train,y_train)
            trainpred,trainprob = clf.predict(X_train),clf.predict_proba(X_train)
            trainacc = metrics.accuracy_score(y_train,trainpred)
            trainauc = metrics.roc_auc_score(y_train,trainprob[:,1])
            print('Training,','ACC:',trainacc,'AUC:',trainauc)

            testpred,testprob = clf.predict(X_val),clf.predict_proba(X_val)
            testacc = metrics.accuracy_score(y_val,testpred)
            testauc = metrics.roc_auc_score(y_val,testprob[:,1])
            print('Test,','ACC:',testacc, 'AUC:',testauc)

            record.loc[i,'Classifiers'] = str(model)
            record.loc[i,'Training ACC'] = trainacc
            record.loc[i,'Test ACC'] = testacc
            record.loc[i,'Training AUC'] = trainauc
            record.loc[i,'Test AUC'] = testauc

            trainpred_report[clfnames[i]] = list(trainpred)
            trainprob_report[clfnames[i]] = list(trainprob[:,1])
            testpred_report[clfnames[i]] = list(testpred)
            testprob_report[clfnames[i]] = list(testprob[:,1])

            best_params = clf.best_params_
            print(best_params)

            trainpred_report.to_excel(save_file_path+'/'+filename+str(k_index)+'/'+'TrainPredict(CV).xlsx',index=False)
            trainprob_report.to_excel(save_file_path+'/'+filename+str(k_index)+'/'+'TrainProbability(CV).xlsx',index=False)
            testpred_report.to_excel(save_file_path+'/'+filename+str(k_index)+'/'+'ValPredict(CV).xlsx',index=False)
            testprob_report.to_excel(save_file_path+'/'+filename+str(k_index)+'/'+'ValProbability(CV).xlsx',index=False)

            record.to_csv(save_file_path+'/'+filename+str(k_index)+'/'+datatime+'_Record(CV).csv',index=False)
            pickle.dump(clf,open(save_file_path+'/'+filename+str(k_index)+'/'+datatime+'_'+str(cv).split(',')[0]+')_'+clfnames[i]+'(CV).dat','wb'))
            print('----------------------------------------------------------------')

        print('\n\n\n')
    k_index += 1




