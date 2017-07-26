import os

def pix():
    files = []
    temp1 = []
    temp2 = []
    final = []

    # imports files lists with structure cat_hpx_#####.fits
    files = os.listdir("/data/des51.b/data/kadrlica/projects/y3q2/v7/cat/")

    # Removes cat_hpx_ from file name
    for i in range(0,len(files),1):
        temp1.append(files[i].replace("cat_hpx_", ""))

    # Removes .fits from file name
    for j in range(0,len(temp1),1):
        temp2.append(temp1[j].replace(".fits", ""))

    # Removes any remaining files that do not have numbers. 
    for k in range(0,len(temp2),1):
        if temp2[k].isdigit() == True:
            final.append(temp2[k])
    
    return final
