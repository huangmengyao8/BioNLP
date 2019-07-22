# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 15:49:54 2019

@author: huangmengyao
"""

import pandas as pd
import numpy as np
import Levenshtein
import os

#设置路径
os.chdir('D:\\Bionlp\\data\\dicts')

#读入词典
taxdicts = pd.read_csv("TAX1_trim.txt",sep='|',encoding='utf-8')
taxdicts = taxdicts.iloc[:,0:2]
taxdicts.columns = ['pre_id','name'] #pre_id  match_name
#taxdicts = taxdicts.iloc[0:10000,0:2]
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

#定义jaro距离函数
def wordjaro_match(df4,df1):
    matched_name = []           # 用来存储概率最大的匹配名称
    matched_id = []             # 用来存储概率最大的匹配id
    matched_probability = []    # 用来存储概率最大的匹配名称 对应的概率
    result = [] 
    
    for word in df1['word'].tolist():
        probability = []
        #print(word)
        for name in df4['name'].tolist(): 
            probability.append(Levenshtein.jaro(word,name)) #word待匹配项。name词典的name
        result.append(probability) # PD
        matched_name.append(df4['name'].tolist()[np.argmax(probability)])
        matched_id.append(df4['pre_id'].tolist()[np.argmax(probability)])
        matched_probability.append(max(probability))
    result = pd.DataFrame(result, columns=df4['name'])
    result = pd.concat([df1.iloc[:, :2], result], axis=1)
    result.to_csv('matrix_jaro.csv', index=False, encoding='utf_8_sig')  
    
    matched_type = []
    real_id = df1['real_id']
    for index in range(len(matched_id)):
        if matched_id[index]== real_id[index]:
            matched_type.append('T')
        elif matched_id[index]!= real_id[index]:
            matched_type.append('F')
        else:
            matched_type.append('-1')
    
    M4 = pd.DataFrame({'matched_name':matched_name,'matched_id':matched_id,'type':matched_type})
    M4 = pd.concat([df1.reset_index(drop=True),M4],axis=1)
    M4.to_csv('jaro.csv', index=False, encoding='utf_8_sig') 
    print(matched_type.count('T')/len(matched_type))
    
    return M4


wordjaro_match(taxdicts,NormInput)