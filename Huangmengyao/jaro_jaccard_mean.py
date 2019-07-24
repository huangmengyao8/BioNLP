# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:00:25 2019

@author: huangmengyao
"""
import pandas as pd
import numpy as np
import jieba
import os
import Levenshtein

#设置路径
os.chdir('D:\\Bionlp\\data\\dicts')

#读入词典
taxdicts = pd.read_csv("TAX1_trim.txt",sep='|',encoding='utf-8')
taxdicts = taxdicts.iloc[:,0:2]
taxdicts.columns = ['pre_id','name'] #pre_id  match_name
taxdicts = taxdicts.iloc[:,0:2]
taxdicts= taxdicts.reset_index(drop = True)

#输入
NormInput  = pd.read_csv("NormInput1.csv",sep='\t')
NormInput = NormInput.loc[NormInput['dict_type']=='NCBI_Taxonomy']
NormInput = NormInput[['dict_id','entity']]
NormInput.columns = ['real_id','word'] #real_id word
NormInput = NormInput.drop_duplicates() 
NormInput = NormInput.dropna(axis=0)
NormInput[['real_id']]=NormInput[['real_id']].astype(int)
NormInput = NormInput.reset_index(drop = True)

#定义Jaccard距离函数
def Jaccrad(model,reference):#model为原句子，reference为候选句子
    terms_reference=jieba.cut(reference)
    terms_model=jieba.cut(model)
    grams_reference=set(terms_reference)
    grams_model=set(terms_model)
    temp=0
    for i in grams_reference:
        if i in grams_model:
            temp=temp+1 #交集
    fenmu=len(grams_model)+len(grams_reference)-temp #并集
    jaccard_coefficient=float(temp/fenmu)
    return jaccard_coefficient

#定义match函数
def match(df4,df1):
    #jaro
    matched_name1 = []           # 用来存储概率最大的匹配名称
    matched_id1 = []             # 用来存储概率最大的匹配id
    matched_probability1 = []    # 用来存储概率最大的匹配名称 对应的概率
    #result = [] 
    #jaccard
    matched_name2 = []           # 用来存储概率最大的匹配名称
    matched_id2 = []             # 用来存储概率最大的匹配id
    matched_probability2 = []    # 用来存储概率最大的匹配名称 对应的概率
    #mean
    matched_name3 = []           # 用来存储概率最大的匹配名称
    matched_id3 = []             # 用来存储概率最大的匹配id
    matched_probability3 = []    # 用来存储概率最大的匹配名称 对应的概率
    
    for word in df1['word'].tolist():
        probability1 = []
        probability2 = []
        probability3 = []
        #print(word)
        for name in df4['name'].tolist():
            jaro=Levenshtein.jaro(word,name)#word待匹配项。name词典的name
            jaccard=Jaccrad(word,name)
            mean=(jaro+jaccard)/2
            
            probability1.append(jaro) 
            probability2.append(jaccard) 
            probability3.append(mean) 
        #result.append(probability) # PD
        matched_name1.append(df4['name'].tolist()[np.argmax(probability1)])
        matched_name2.append(df4['name'].tolist()[np.argmax(probability2)])
        matched_name3.append(df4['name'].tolist()[np.argmax(probability3)])
        
        matched_id1.append(df4['pre_id'].tolist()[np.argmax(probability1)])
        matched_id2.append(df4['pre_id'].tolist()[np.argmax(probability2)])
        matched_id3.append(df4['pre_id'].tolist()[np.argmax(probability3)])
        
        matched_probability1.append(max(probability1))
        matched_probability2.append(max(probability2))
        matched_probability3.append(max(probability3))
    #result = pd.DataFrame(result, columns=df4['name'])
    #result = pd.concat([df1.iloc[:, :2], result], axis=1)
    #result.to_csv('matrix_jaro.csv', index=False, encoding='utf_8_sig')  
    
#    matched_type = []
#    real_id = df1['real_id']
#    for index in range(len(matched_id)):
#        if matched_id[index]== real_id[index]:
#            matched_type.append('T')
#        elif matched_id[index]!= real_id[index]:
#            matched_type.append('F')
#        else:
#            matched_type.append('-1')
    
    M4 = pd.DataFrame({'matched_name_jaro':matched_name1,
                       'matched_name_Jaccard':matched_name2,
                       'matched_name_mean':matched_name3,
                       'matched_id_jaro':matched_id1,
                       'matched_id_Jaccard':matched_id2,
                       'matched_id_mean':matched_id3,
                       'matched_probability_jaro':matched_probability1,
                       'matched_probability_Jaccard':matched_probability2,
                       'matched_probability_mean':matched_probability3})
    M4 = pd.concat([df1.reset_index(drop=True),M4],axis=1)
    M4.to_csv('jaro_jaccard_mean.csv', index=False, encoding='utf_8_sig',sep='\t') 
    #print(matched_type.count('T')/len(matched_type))
    
    return M4

match(taxdicts,NormInput)