import os
#根据葵花8号文件20,30,40,50，查找对应的风云4a数据19,30,38,53
# h8 = /data/h08/202001/01/00/NC_H08_20200101_0030_L2CLP010_FLDK.02401_02401.nc
# fy4a = /data/fy4a/20200101/FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_20200101003000_20200101003417_4000M_V0001.HDF
def h8_find_fy4a(h08file, fy4a_root):
    h08date = h08file.split('_')[2]
    h08time = h08file.split('_')[3]
    hour = h08time[:2]
    minu = int(h08time[2:])
    if minu ==0 or minu ==10:
        return None
    elif minu ==20:
        fy4atime = h08date+hour+'1918'
        end = h08date+hour+'2335'
    elif minu == 30:
        fy4atime = h08date+hour+'3000'
        end = h08date+hour+'3417'
    elif minu == 40:
        fy4atime = h08date+hour+'3836'
        end = h08date+hour+'4253'
    elif minu == 50:
        fy4atime = h08date+hour+'4918'   
        end = h08date+hour+'5335' 
    fy4afilepath = os.path.join(fy4a_root,h08date)   
    fy4afilename = os.path.join(fy4afilepath, 'FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_%s_%s_4000M_V0001.HDF'%(fy4atime,end))   
    return fy4afilename

