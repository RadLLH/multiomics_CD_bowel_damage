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
otu = otu2_normalized

cov = cov2
mre = mre2


# 找到所有以"QC"开头的列
qc_Q6002 = Q6002.loc[Q600_select,].filter(regex='^QC')
# 计算这些列的均值
qc_Q6002_means = qc_Q6002.mean(axis=1)
# 打印质控均值
print("Quality Control Means:")
print(qc_Q6002_means)
Q600 = Q6002.loc[Q600_select,].div(qc_Q6002_means, axis=0)

# 找到所有以"QC"开头的列
qc_Q3002 = Q3002.loc[Q300_select,].filter(regex='^QC')
# 计算这些列的均值
qc_Q3002_means = qc_Q3002.mean(axis=1)
# 打印质控均值
print("Quality Control Means:")
print(qc_Q3002_means)
Q300 = Q3002.loc[Q300_select,].div(qc_Q3002_means, axis=0)

treat = pd.read_csv('../data/treat.csv',index_col=0)
group_df = treat.loc[treat['Group'] == 'CD',:]
idx_raw1 = (otu.columns & group_df.index) & (Q600.columns & Q300.columns)
y_ = group_df.loc[idx_raw1,'Group_BD']
y = np.array(y_=='BD2',dtype='int')
x_idx = group_df.loc[idx_raw1,:].index

X_test_cov = cov.loc[cov_select,x_idx].transpose()
print(X_test_cov.shape)
X_test_otu = otu.loc[otu_select,x_idx].transpose()
print(X_test_otu.shape)
X_test_Q300 = Q300.loc[Q300_select,x_idx].transpose()
print(X_test_Q300.shape)
X_test_Q600 = Q600.loc[Q600_select,x_idx].transpose()
print(X_test_Q600.shape)
X_test_mre = mre.loc[mre_select,x_idx].transpose()
print(X_test_mre.shape)

X_test_concat_cov = pd.concat([X_test_otu,X_test_Q300,X_test_Q600,X_test_mre,X_test_cov],axis=1)

X_test_otu_Q300_Q600_cov = pd.concat([X_test_otu,X_test_Q300,X_test_Q600,X_test_cov],axis=1)

X_test_otu_cov = pd.concat([X_test_otu,X_test_cov],axis=1)

X_test_Q300_cov = pd.concat([X_test_Q300,X_test_cov],axis=1)

X_test_Q600_cov = pd.concat([X_test_Q600,X_test_cov],axis=1)

X_test_mre_cov = pd.concat([X_test_mre,X_test_cov],axis=1)

save_file_path = './result_multi_algorithm'

for k_index in range(5):
    for filename,(X,Y) in zip(['cov','otu','Q300','Q600','MRE','otu_Q300_Q600','combined'],
                              [(X_test_cov,y),
                               (X_test_otu_cov,y),
                               (X_test_Q300_cov,y),
                               (X_test_Q600_cov,y),
                               (X_test_mre_cov,y),
                               (X_test_otu_Q300_Q600_cov,y),
                               (X_test_concat_cov,y)]):
        os.makedirs(save_file_path+'/'+filename+str(k_index), exist_ok=True)

        record = pd.DataFrame(columns=['Classifiers'])
        if Path.exists(Path(save_file_path+'/'+filename+str(k_index)+'/'+datatime+'_Record(CV).csv')):
            record = pd.read_csv(save_file_path+'/'+filename+str(k_index)+'/'+datatime+'_Record(CV).csv')

        testpred_report = pd.DataFrame({'id':X.index,'Group':Y})
        testprob_report = pd.DataFrame({'id':X.index,'Group':Y})
        print('============================',filename,'============================')
        clfnames = ['DT','RFC','XGBC','LR','Linear_SVM','Rbf_SVM','KNN']
        i = 0
        for classifier in clfnames:


            clf = pickle.load(open(save_file_path+'/'+filename+str(k_index)+'/2024-11-19_StratifiedKFold(n_splits=5)_'+classifier+'(CV).dat', 'rb'))

            testpred,testprob = clf.predict(X),clf.predict_proba(X)
            testacc = metrics.accuracy_score(Y,testpred)
            testauc = metrics.roc_auc_score(Y,testprob[:,1])
            print('Test,','ACC:',testacc, 'AUC:',testauc)

            record.loc[i,'Classifiers'] = 'RFC'
            record.loc[i,'Test ACC'] = testacc


            testpred_report[clfnames[i]] = list(testpred)
            testprob_report[clfnames[i]] = list(testprob[:,1])
            i += 1


        testpred_report.to_excel(save_file_path+'/'+filename+str(k_index)+'/'+'TestPredict(CV).xlsx',index=False)
        testprob_report.to_excel(save_file_path+'/'+filename+str(k_index)+'/'+'TestProbability(CV).xlsx',index=False)

        print('----------------------------------------------------------------')

        print('\n\n\n')

#ROC曲线
clfnames = ['DT','RFC','XGBC','LR','Linear_SVM','Rbf_SVM','KNN']
for classifier in clfnames:
    prob_otu_series = []
    prob_cov_series = []
    prob_Q300_series = []
    prob_Q600_series = []
    prob_mre_series = []
    prob_concat_2_series = []
    prob_concat_series = []

    y_true_cov_series = []
    y_true_otu_series = []
    y_true_Q300_series = []
    y_true_Q600_series = []
    y_true_mre_series = []
    y_true_concat_2_series = []
    y_true_concat_series = []


    prob_series = []
    y_true_series = []
    for k_index in range(5):
        prob_cov = pd.read_excel(save_file_path+'/cov'+str(k_index)+'/TestProbability(CV).xlsx')
        prob_otu = pd.read_excel(save_file_path+'/otu'+str(k_index)+'/TestProbability(CV).xlsx')
        prob_Q300 = pd.read_excel(save_file_path+'/Q300'+str(k_index)+'/TestProbability(CV).xlsx')
        prob_Q600 = pd.read_excel(save_file_path+'/Q600'+str(k_index)+'/TestProbability(CV).xlsx')
        prob_mre = pd.read_excel(save_file_path+'/MRE'+str(k_index)+'/TestProbability(CV).xlsx')
        prob_concat = pd.read_excel(save_file_path+'/combined'+str(k_index)+'/TestProbability(CV).xlsx')
        prob_concat_2 = pd.read_excel(save_file_path+'/otu_Q300_Q600'+str(k_index)+'/TestProbability(CV).xlsx')

        prob_cov_series.extend(prob_cov.loc[:,classifier])
        prob_otu_series.extend(prob_otu.loc[:,classifier])
        prob_Q300_series.extend(prob_Q300.loc[:,classifier])
        prob_Q600_series.extend(prob_Q600.loc[:,classifier])
        prob_mre_series.extend(prob_mre.loc[:,classifier])
        prob_concat_2_series.extend(prob_concat_2.loc[:,classifier])
        prob_concat_series.extend(prob_concat.loc[:,classifier])


        y_true_cov_series.extend(prob_cov.loc[:,'Group'])
        y_true_otu_series.extend(prob_otu.loc[:,'Group'])
        y_true_Q300_series.extend(prob_Q300.loc[:,'Group'])
        y_true_Q600_series.extend(prob_Q600.loc[:,'Group'])
        y_true_mre_series.extend(prob_mre.loc[:,'Group'])
        y_true_concat_2_series.extend(prob_concat_2.loc[:,'Group'])
        y_true_concat_series.extend(prob_concat.loc[:,'Group'])



    fpr_binary = dict()
    tpr_binary = dict()
    roc_auc_binary = dict()
    prob_series = [prob_cov_series,
                   prob_otu_series,
                   prob_Q300_series,
                   prob_Q600_series,
                   prob_mre_series,
                   prob_concat_2_series,
                   prob_concat_series]
    y_true_series = [y_true_cov_series,
                     y_true_otu_series,
                     y_true_Q300_series,
                     y_true_Q600_series,
                     y_true_mre_series,
                     y_true_concat_2_series,
                     y_true_concat_series]

    for i,(prob,y_label) in enumerate(zip(prob_series,y_true_series)):

        fpr_binary[i], tpr_binary[i], _ = metrics.roc_curve(y_label, prob)
        roc_auc_binary[i] = metrics.roc_auc_score(y_label, prob)

    plt.figure(dpi=300,figsize=(5,4))
    ax1=plt.gca()
    ax1.patch.set_facecolor("white")
    ax1.patch.set_alpha(0.95)
    plt.grid(color='w',
             linestyle='--',
             linewidth=1,
             alpha=0.3)
    lw = 1
    colors = ['#42d4f4','#3cb44b','#ffe119','#e6194B','#4363d8','#f58231','#911eb4']
    for i, color in zip(range(7), colors):

        plt.plot(fpr_binary[i], tpr_binary[i], color=color, linestyle='-', lw=lw,
                 label='{0} (AUC = {1:0.3f})'
                       ''.format(['COV','Microbiota','Fecal','Serum','MRE','Microbiota+Fecal+Serum\n','Microbiota+Fecal+Serum+MRE\n'][i], roc_auc_binary[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([-0.02, 1.02])
    plt.ylim([-0.02, 1.02])
    plt.xlabel('1 - Specificity',fontsize=10)
    plt.ylabel('Sensitivity',fontsize=10)
    plt.title('Receiver Operating Characteristic Curve',fontsize=10)
    plt.legend(loc="lower right",fontsize=7)
    plt.savefig(save_file_path+'/'+classifier+'_test_roc.pdf', format='pdf', dpi=None, bbox_inches='tight')
    #plt.show()

for classifier in clfnames:
    prob_otu_series = []
    prob_cov_series = []
    prob_Q300_series = []
    prob_Q600_series = []
    prob_mre_series = []
    prob_concat_2_series = []
    prob_concat_series = []

    y_true_cov_series = []
    y_true_otu_series = []
    y_true_Q300_series = []
    y_true_Q600_series = []
    y_true_mre_series = []
    y_true_concat_2_series = []
    y_true_concat_series = []


    prob_series = []
    y_true_series = []
    for k_index in range(5):
        prob_cov = pd.read_excel(save_file_path+'/cov'+str(k_index)+'/ValProbability(CV).xlsx')
        prob_otu = pd.read_excel(save_file_path+'/otu'+str(k_index)+'/ValProbability(CV).xlsx')
        prob_Q300 = pd.read_excel(save_file_path+'/Q300'+str(k_index)+'/ValProbability(CV).xlsx')
        prob_Q600 = pd.read_excel(save_file_path+'/Q600'+str(k_index)+'/ValProbability(CV).xlsx')
        prob_mre = pd.read_excel(save_file_path+'/MRE'+str(k_index)+'/ValProbability(CV).xlsx')
        prob_concat = pd.read_excel(save_file_path+'/combined'+str(k_index)+'/ValProbability(CV).xlsx')
        prob_concat_2 = pd.read_excel(save_file_path+'/otu_Q300_Q600'+str(k_index)+'/ValProbability(CV).xlsx')

        prob_cov_series.extend(prob_cov.loc[:,classifier])
        prob_otu_series.extend(prob_otu.loc[:,classifier])
        prob_Q300_series.extend(prob_Q300.loc[:,classifier])
        prob_Q600_series.extend(prob_Q600.loc[:,classifier])
        prob_mre_series.extend(prob_mre.loc[:,classifier])
        prob_concat_2_series.extend(prob_concat_2.loc[:,classifier])
        prob_concat_series.extend(prob_concat.loc[:,classifier])


        y_true_cov_series.extend(prob_cov.loc[:,'Group'])
        y_true_otu_series.extend(prob_otu.loc[:,'Group'])
        y_true_Q300_series.extend(prob_Q300.loc[:,'Group'])
        y_true_Q600_series.extend(prob_Q600.loc[:,'Group'])
        y_true_mre_series.extend(prob_mre.loc[:,'Group'])
        y_true_concat_2_series.extend(prob_concat_2.loc[:,'Group'])
        y_true_concat_series.extend(prob_concat.loc[:,'Group'])



    fpr_binary = dict()
    tpr_binary = dict()
    roc_auc_binary = dict()
    prob_series = [prob_cov_series,
                   prob_otu_series,
                   prob_Q300_series,
                   prob_Q600_series,
                   prob_mre_series,
                   prob_concat_2_series,
                   prob_concat_series]
    y_true_series = [y_true_cov_series,
                     y_true_otu_series,
                     y_true_Q300_series,
                     y_true_Q600_series,
                     y_true_mre_series,
                     y_true_concat_2_series,
                     y_true_concat_series]

    for i,(prob,y_label) in enumerate(zip(prob_series,y_true_series)):

        fpr_binary[i], tpr_binary[i], _ = metrics.roc_curve(y_label, prob)
        roc_auc_binary[i] = metrics.roc_auc_score(y_label, prob)

    plt.figure(dpi=300,figsize=(5,4))
    ax1=plt.gca()
    ax1.patch.set_facecolor("white")
    ax1.patch.set_alpha(0.95)
    plt.grid(color='w',
             linestyle='--',
             linewidth=1,
             alpha=0.3)
    lw = 1
    colors = ['#42d4f4','#3cb44b','#ffe119','#e6194B','#4363d8','#f58231','#911eb4']
    for i, color in zip(range(7), colors):

        plt.plot(fpr_binary[i], tpr_binary[i], color=color, linestyle='-', lw=lw,
                 label='{0} (AUC = {1:0.3f})'
                       ''.format(['COV','Microbiota','Fecal','Serum','MRE','Microbiota+Fecal+Serum\n','Microbiota+Fecal+Serum+MRE\n'][i], roc_auc_binary[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([-0.02, 1.02])
    plt.ylim([-0.02, 1.02])
    plt.xlabel('1 - Specificity',fontsize=10)
    plt.ylabel('Sensitivity',fontsize=10)
    plt.title('Receiver Operating Characteristic Curve',fontsize=10)
    plt.legend(loc="lower right",fontsize=7)
    plt.savefig(save_file_path+'/'+classifier+'_val_roc.pdf', format='pdf', dpi=None, bbox_inches='tight')

for classifier in clfnames:
    prob_otu_series = []
    prob_cov_series = []
    prob_Q300_series = []
    prob_Q600_series = []
    prob_mre_series = []
    prob_concat_2_series = []
    prob_concat_series = []

    y_true_cov_series = []
    y_true_otu_series = []
    y_true_Q300_series = []
    y_true_Q600_series = []
    y_true_mre_series = []
    y_true_concat_2_series = []
    y_true_concat_series = []


    prob_series = []
    y_true_series = []
    for k_index in range(5):
        prob_cov = pd.read_excel(save_file_path+'/cov'+str(k_index)+'/TrainProbability(CV).xlsx')
        prob_otu = pd.read_excel(save_file_path+'/otu'+str(k_index)+'/TrainProbability(CV).xlsx')
        prob_Q300 = pd.read_excel(save_file_path+'/Q300'+str(k_index)+'/TrainProbability(CV).xlsx')
        prob_Q600 = pd.read_excel(save_file_path+'/Q600'+str(k_index)+'/TrainProbability(CV).xlsx')
        prob_mre = pd.read_excel(save_file_path+'/MRE'+str(k_index)+'/TrainProbability(CV).xlsx')
        prob_concat = pd.read_excel(save_file_path+'/combined'+str(k_index)+'/TrainProbability(CV).xlsx')
        prob_concat_2 = pd.read_excel(save_file_path+'/otu_Q300_Q600'+str(k_index)+'/TrainProbability(CV).xlsx')

        prob_cov_series.extend(prob_cov.loc[:,classifier])
        prob_otu_series.extend(prob_otu.loc[:,classifier])
        prob_Q300_series.extend(prob_Q300.loc[:,classifier])
        prob_Q600_series.extend(prob_Q600.loc[:,classifier])
        prob_mre_series.extend(prob_mre.loc[:,classifier])
        prob_concat_2_series.extend(prob_concat_2.loc[:,classifier])
        prob_concat_series.extend(prob_concat.loc[:,classifier])


        y_true_cov_series.extend(prob_cov.loc[:,'Group'])
        y_true_otu_series.extend(prob_otu.loc[:,'Group'])
        y_true_Q300_series.extend(prob_Q300.loc[:,'Group'])
        y_true_Q600_series.extend(prob_Q600.loc[:,'Group'])
        y_true_mre_series.extend(prob_mre.loc[:,'Group'])
        y_true_concat_2_series.extend(prob_concat_2.loc[:,'Group'])
        y_true_concat_series.extend(prob_concat.loc[:,'Group'])



    fpr_binary = dict()
    tpr_binary = dict()
    roc_auc_binary = dict()
    prob_series = [prob_cov_series,
                   prob_otu_series,
                   prob_Q300_series,
                   prob_Q600_series,
                   prob_mre_series,
                   prob_concat_2_series,
                   prob_concat_series]
    y_true_series = [y_true_cov_series,
                     y_true_otu_series,
                     y_true_Q300_series,
                     y_true_Q600_series,
                     y_true_mre_series,
                     y_true_concat_2_series,
                     y_true_concat_series]

    for i,(prob,y_label) in enumerate(zip(prob_series,y_true_series)):

        fpr_binary[i], tpr_binary[i], _ = metrics.roc_curve(y_label, prob)
        roc_auc_binary[i] = metrics.roc_auc_score(y_label, prob)

    plt.figure(dpi=300,figsize=(5,4))
    ax1=plt.gca()
    ax1.patch.set_facecolor("white")
    ax1.patch.set_alpha(0.95)
    plt.grid(color='w',
             linestyle='--',
             linewidth=1,
             alpha=0.3)
    lw = 1
    colors = ['#42d4f4','#3cb44b','#ffe119','#e6194B','#4363d8','#f58231','#911eb4']
    for i, color in zip(range(7), colors):

        plt.plot(fpr_binary[i], tpr_binary[i], color=color, linestyle='-', lw=lw,
                 label='{0} (AUC = {1:0.3f})'
                       ''.format(['COV','Microbiota','Fecal','Serum','MRE','Microbiota+Fecal+Serum\n','Microbiota+Fecal+Serum+MRE\n'][i], roc_auc_binary[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([-0.02, 1.02])
    plt.ylim([-0.02, 1.02])
    plt.xlabel('1 - Specificity',fontsize=10)
    plt.ylabel('Sensitivity',fontsize=10)
    plt.title('Receiver Operating Characteristic Curve',fontsize=10)
    plt.legend(loc="lower right",fontsize=7)
    plt.savefig(save_file_path+'/'+classifier+'_train_roc.pdf', format='pdf', dpi=None, bbox_inches='tight')
