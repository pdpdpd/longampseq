#!/bin/bash

Filename=( BCL11A HBG1 OT18 R02 R66S HBB )
# BCL11A-1_501.csv   HBG1-2_502.csv         R02-1_502.csv
# BCL11A-1_502.csv   HBG1_all.csv           R02-2_501.csv
# BCL11A-2_501.csv   HBG1-UT_501.csv        R02-2_502.csv
# BCL11A-2_502.csv   HBG1-UT_502.csv        R66SHIFI-1_501.csv
# BCL11A_all.csv     HBG1-UT_all.csv        R66SHIFI-1_502.csv
# BCL11A-UT_501.csv  OT18_all.csv           R66SHIFI-2_501.csv
# BCL11A-UT_502.csv  OT18-R66SWT-1_501.csv  R66SHIFI-2_502.csv
# BCL11A-UT_all.csv  OT18-R66SWT-1_502.csv  R66SHIFI-2_large.csv
# HBB-UT_501.csv     OT18-R66SWT-2_501.csv  R66SHIFI_all.csv
# HBB-UT_502.csv     OT18-R66SWT-2_502.csv  R66SWT-1_501.csv
# HBB-UT_all.csv     OT18-UT_501.csv        R66SWT-1_502.csv
# HBG1-1_501.csv     OT18-UT_502.csv
for x in "${Filename[@]}"
do
    # echo $x
    if [ "$x" == BCL11A ]; then
        cutsite=60722401
        length=4267
    elif [ "$x" == HBG1 ]; then
        cutsite=2759
        length=6718
    elif [ "$x" == OT18 ]; then
        cutsite=59137387
        length=5045
    elif [ "$x" == R02 ]; then
        cutsite=5248214
        length=5465
    elif [ "$x" == R66S ]; then
        cutsite=5248229
        length=5465
    fi

    # echo $cutsite
    # echo $length
    for file in ${x}*.csv
    do
        echo ${file}
        echo ${cutsite}
        echo ${length}
        # python ~/Documents/scripts/0303_bedfile.py ${file} ${file/.csv/largedel_output.csv} ${file/.csv/largedel_group.csv} ${length}
        #can add more restrictions if excluding large dels not spanning cut site
        #draw figure
        python ~/Documents/scripts/0303_longampfigures.py ./output/${file/.csv/largedel_output.csv} ${cutsite}
        # mv ${file/.csv/largedel_output.csv} ./output
        # mv ${file/.csv/largedel_group.csv} ./output
        # mv ${file/.csv/largedel_output.eps} ./output
    done

    if [ "$x" == HBB ]; then
        for file in ${x}*.csv
        do
            python ~/Documents/scripts/0303_bedfile.py ${file} ${file/.csv/R02_largedel_output.csv} ${file/.csv/R02_largedel_group.csv} 5465
            python ~/Documents/scripts/0303_bedfile.py ${file} ${file/.csv/R66S_largedel_output.csv} ${file/.csv/R66S_largedel_group.csv} 5465
            #can add more restrictions if excluding large dels not spanning cut site
            #draw figure
            python ~/Documents/scripts/0303_longampfigures.py ${file/.csv/R02_largedel_output.csv} 5248214
            python ~/Documents/scripts/0303_longampfigures.py ${file/.csv/R66S_largedel_output.csv} 5248229
        done
    fi

done
