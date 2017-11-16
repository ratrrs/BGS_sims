import pandas as pd
import os
import sys


results_dir = sys.argv[1]#'../results'
model = sys.argv[2]#'tennessen'
def combine_results(input_dir,stats,selection,model):
    df_list = []
    for file in os.listdir(input_dir+ os.sep + model):
        if file.endswith("%s_%s.csv"%(selection,stats)):
            df_list.append(pd.read_csv(input_dir+os.sep+model+os.sep+file))
    combined = pd.concat(df_list)
    out_name = input_dir+os.sep+model+'_'+selection+'_'+ stats + '.csv'
    combined.to_csv(out_name,index=False)
    return(combined)


stats_list = ['pi','xi','tajD']
for i in stats_list:
    neutral = combine_results(results_dir,i,'neutral',model)
    bgs = combine_results(results_dir,i,'bgs',model)


