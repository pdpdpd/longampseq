# python2.7
# ref: hg19

import sys
import math
import os
import numpy as np
import csv
import pandas as pd
import pickle
import argparse
import re


def main():
    # use full path here
    bed_path = sys.argv[1]
    # bed_path = "/home/yp11/Documents/0219longamp/test_mergefiltered_2+.csv"

    df = pd.read_csv(bed_path, names=['chr', 'start', 'end', 'read_ID', 'score', 'strand'])

    # chr11	59136655	59136956	M04808:132:000000000-CTP75:1:2107:19955:2557	60	-
    # chr11	59138619	59138657	M04808:132:000000000-CTP75:1:2107:19955:2557	60	-
    # df.sort_values(by=['read_ID'])
    # ID_list = df.read_ID.unique()

    # Count the frequency
    id_counts = df['read_ID'].value_counts()
    df[df['read_ID'].isin(id_counts[id_counts == 2].index)]
    id_list = df.read_ID.unique()

    column_names = ['read_ID', 'start', 'length']
    df_output = pd.DataFrame(columns=column_names)

    for id in id_list:
        pos_list = []
        chr_list = []
        df_id = df.loc[df['read_ID'] == id]
        for row in df_id.itertuples():
            chr_list.append(row.chr)
            pos_list.append(row.start)
            pos_list.append(row.end)
        if len(pos_list) == 4 and chr_list[0] == chr_list[1] and pos_list[2] - pos_list[1] >= 200:
            pos_list.sort()
            df_output = df_output.append({'read_ID': id, 'start': pos_list[1], 'length': pos_list[2] - pos_list[1]},
                                         ignore_index=True)

    df_output.to_csv(sys.argv[2], index=False)
    # df_output.to_csv("/home/yp11/Documents/0219longamp/largedel_output.csv", index=False)

    total_large_deletion = len(df_output.index)
    #

    # Open a file with access mode 'a'
    file_object = open('log.txt', 'a')
    file_object.write(bed_path)
    file_object.write('\n')
    file_object.write(str(total_large_deletion))
    file_object.write('\n')
    file_object.close()

    # consolidation
    df_conso = df_output[['start', 'length']]
    df_conso_group = df_conso.groupby(df_conso.columns.tolist()).size().reset_index(). \
        rename(columns={0: 'repeat_num'})
    df_conso_group.to_csv(sys.argv[3], index=False)


if __name__ == '__main__':
    main()
