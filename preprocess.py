import numpy as np
# from fy4a import FY4A_AGRI_L1
# import dboxio
import netCDF4 as nc
from h8_find_fy4a import h8_find_fy4a
import os, sys
import pandas as pd
import gdal
import osr
from fy4a import AGRI_L1
import time

BLOCK_SIZE = 256
OVERLAP_SIZE = 0
# dboxio.SetThreadLocalSymbol("DBOX_DEFAULT_PIXELS", str(BLOCK_SIZE))
# 各分辨率文件包含的通道号
CONTENTS = {'0500M': ('Channel02',),
            '1000M': ('Channel01', 'Channel02', 'Channel03'),
            '2000M': tuple(['Channel' + "%02d" % (x) for x in range(1, 8)]),
            '4000M': tuple(['Channel' + "%02d" % (x) for x in range(1, 15)])}

# geo_range = '5, 55, 70, 140, 0.05'
# # geo_range to geo_transform 
# minx = float(geo_range.split(',')[2])
# maxy = float(geo_range.split(',')[1])
# res = float(geo_range.split(',')[-1])
geo_range = [5, 54.95, 70, 139.95, 0.05] 
minx = geo_range[2]
maxy = geo_range[1]
res = geo_range[4]

fy4a_gt = (minx, res, 0.0, maxy, 0.0, (-1) * res)

NP2GDAL_CONVERSION = {
  "uint8": 1,
  "int8": 1,
  "uint16": 2,
  "int16": 3,
  "uint32": 4,
  "int32": 5,
  "float32": 6,
  "float64": 7,
  "complex64": 10,
  "complex128": 11,
}

# 数据根目录
fy4a_root = '/data/fy4a'
h08_root = '/data/h08'
# fy4a_tif = '/mnt/win/data/cloud/fy4a_tif'
# h08_tif = '/mnt/win/data/cloud/h08_tif'
# fy4a_tiles = '/mnt/win/data/cloud/fy4a_tiles'
# h08_tiles = '/mnt/win/data/cloud/h08_tiles'

fy4a_tif = '/tmp/fy4a_tif'
h08_tif = '/tmp/h08_tif'
fy4a_tiles = '/data/fy4a_tiles'
h08_tiles = '/data/h08_tiles'


def h08array_to_tiff(data, proj, gt, dst_nbands, dst_file, dtype):
    xsize = data.shape[0]
    ysize = data.shape[1]
    dst_format = 'GTiff'
#     dst_nbands = 1
    dst_datatype = gdal.GDT_Byte

    driver = gdal.GetDriverByName(dst_format)
    dst_ds = driver.Create(dst_file, ysize, xsize, dst_nbands, dst_datatype)
    dst_ds.SetGeoTransform(gt)
    dst_ds.SetProjection(proj)
#     dst_ds.GetRasterBand(1).WriteArray(data)
    
    if dst_nbands == 1:
        dst_ds.GetRasterBand(1).WriteArray(data)
    else:
        for i in range(dst_nbands):
            dst_ds.GetRasterBand(i + 1).WriteArray(data[i, :, :])
    del dst_ds


def read_fy4a_arr(h5name):
    
# file_path = '/data/fy4a/20200101/FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_20200101231918_20200101232335_4000M_V0001.HDF'
# geo_desc = [5, 54.95, 70, 139.95, 0.05]  # 顺序为南、北、西、东、分辨率，即[lat_s, lat_n, lon_w, lon_e, resolution]

    file = AGRI_L1(h5name, geo_range)
#     x =file.dataset['NOMChannel01']
#     type = x.dtype
    channels = CONTENTS['4000M']
    bands_list = []
    for channel in channels:
        channel_num = int(channel[-2:])
        if channel_num <= 6:
            bands_list.append(file.extract(channel))
        if channel_num > 6:
            bands_list.append(file.extract(channel, calibration='brightness_temperature'))
        
    imgarr = np.stack(band for band in bands_list)
#     print(imgarr.shape)    
    
#     fy4a_agri_l1 = FY4A_AGRI_L1(h5name)
# 
#     channels = CONTENTS['4000M']
#     for channel in channels:
#         fy4a_agri_l1.extract(channel, geo_range)
        
#     imgarr = np.stack(fy4a_agri_l1.channels[channel] for channel in channels)
#     print(imgarr.shape)
    
    return imgarr


def read_h08_label(h08filename, dst_dboxfile):
    ds = nc.Dataset(h08filename, 'r')
    CLTYPEarr = ds.variables['CLTYPE'][:].data
    
    Lon = ds.variables['longitude'][:]
    Lat = ds.variables['latitude'][:]
    LonMin, LatMax, LonMax, LatMin = [Lon.min(), Lat.max(), Lon.max(), Lat.min()] 
    # 分辨率计算
    N_Lat = len(Lat) 
    N_Lon = len(Lon)
    Lon_Res = (LonMax - LonMin) / (float(N_Lon) - 1)
    Lat_Res = (LatMax - LatMin) / (float(N_Lat) - 1)
    # 设置影像的显示范围
    # -Lat_Res一定要是-的
    geotransform = (LonMin, Lon_Res, 0, LatMax, 0, -Lat_Res)
    xsize = N_Lon
    ysize = N_Lat
    nbands = 1
    dtype = NP2GDAL_CONVERSION[str(CLTYPEarr.dtype)]
    
    dboxio.CreateDBoxDataset(dst_dboxfile, dboxio.EPSG_4326, geotransform, xsize, ysize, bands=nbands, datatype=dtype)
    return CLTYPEarr


def read_h08_by_tile(CLTYPEarr, geotransform):
    xsize = CLTYPEarr.shape[1]
    ysize = CLTYPEarr.shape[0]
    rnum_tile = int((ysize - tile_height) / (tile_height - OVERLAP_SIZE)) + 1
    cnum_tile = int((xsize - tile_width) / (tile_width - OVERLAP_SIZE)) + 1
    
    tile_list = []
    tile_list = []
    for i in range(rnum_tile + 1):
        for j in range(cnum_tile + 1):
            xoff = 0 + (tile_width - OVERLAP_SIZE) * j
            yoff = 0 + (tile_height - OVERLAP_SIZE) * i
            # the last column                 
            if j == cnum_tile:
                xoff = xsize - tile_width
            # the last row
            if i == rnum_tile:
                yoff = ysize - tile_height
            print("the row and column of tile is :", xoff, yoff)
            
            data = band.ReadAsArray(xoff, yoff, tile_width, tile_height)
            if np.all(data == noDataValue):
                continue
            feature.append(data)
            minx = geotransform[0] + xoff * geotransform[1]
            maxy = geotransform[3] + yoff * geotransform[5]
            
            tile_list.append([minx, maxy])
    
    
def gen_tile_bbox(h08file, BLOCK_SIZE, OVERLAP_SIZE):
    print('the image is :', h08file)
    dataset = gdal.Open(h08file)
    if dataset is None:
        print("Failed to open file: " + h08file)
        sys.exit(1)
    band = dataset.GetRasterBand(1)
    xsize = dataset.RasterXSize
    ysize = dataset.RasterYSize
    proj = dataset.GetProjection()
    geotrans = dataset.GetGeoTransform()
    tile_gt = list(geotrans)
    noDataValue = band.GetNoDataValue()
    
    minx_wgs, maxy_wgs = GeomTrans(proj, 'EPSG:4326').transform_point([geotrans[0], geotrans[3]])
    maxx_wgs, miny_wgs = GeomTrans(proj, 'EPSG:4326').transform_point([geotrans[0] + xsize * geotrans[1], geotrans[3] + ysize * geotrans[5]])
    region_bbox = [minx_wgs, maxy_wgs, maxx_wgs, miny_wgs]
    
    rnum_tile = int((ysize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    cnum_tile = int((xsize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    print('the number of tile is :', rnum_tile * cnum_tile)
    
    wgs_bbox_list = []
    
    for i in range(rnum_tile + 1):
        print(i)
        for j in range(cnum_tile + 1):
            xoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * j
            yoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * i
            # the last column                 
            if j == cnum_tile:
                xoff = xsize - BLOCK_SIZE
            # the last row
            if i == rnum_tile:
                yoff = ysize - BLOCK_SIZE
            print("the row and column of tile is :", xoff, yoff)
            
            data = band.ReadAsArray(xoff, yoff, BLOCK_SIZE, BLOCK_SIZE)
            if np.all(data == noDataValue):
                continue
               
            tile_gt[0] = geotrans[0] + xoff * geotrans[1]
            tile_gt[3] = geotrans[3] + yoff * geotrans[5]
            
            minx = tile_gt[0]
            maxy = tile_gt[3]
            maxx = tile_gt[0] + BLOCK_SIZE * geotrans[1]
            miny = tile_gt[3] + BLOCK_SIZE * geotrans[5]
            
            minx_wgs, maxy_wgs = GeomTrans(proj, 'EPSG:4326').transform_point([minx, maxy])
            maxx_wgs, miny_wgs = GeomTrans(proj, 'EPSG:4326').transform_point([maxx, miny])
            
            wgs_bbox_list.append([minx_wgs, maxy_wgs, maxx_wgs, miny_wgs, i, j])
              
    return wgs_bbox_list, rnum_tile, cnum_tile, region_bbox


def tif2tiles(fy4afile, h08file, BLOCK_SIZE, OVERLAP_SIZE):
    h08ds = gdal.Open(h08file)
    fy4ads = gdal.Open(fy4afile)
    if h08ds is None:
        print("Failed to open file: " + h08file)
        sys.exit(1)
    if fy4ads is None:
        print("Failed to open file: " + fy4afile)
        sys.exit(1)
    band = h08ds.GetRasterBand(1)
    xsize = h08ds.RasterXSize
    ysize = h08ds.RasterYSize
    proj = h08ds.GetProjection()
    noDataValue = band.GetNoDataValue()
    dataType = band.DataType
    geotrans = h08ds.GetGeoTransform()
    tile_gt = list(geotrans)
    noDataValue = band.GetNoDataValue()
    
    fyband = fy4ads.GetRasterBand(1)
    fydataType = fyband.DataType
    
    rnum_tile = int((ysize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    cnum_tile = int((xsize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    print('the number of tile is :', rnum_tile * cnum_tile)
    
    for i in range(rnum_tile + 1):
        print(i)
        for j in range(cnum_tile + 1):
            h08_tile_name = h08file.split('/')[-1].split('.')[0] + '_' + str(i) + '_' + str(j) + '.tif'
            h08_tile_file = os.path.join(h08_tiles, h08_tile_name)
            
            fy4a_tile_name = fy4afile.split('/')[-1].split('.')[0] + '_' + str(i) + '_' + str(j) + '.tif'
            fy4a_tile_file = os.path.join(fy4a_tiles, fy4a_tile_name)
            
            if os.path.exists(h08_tile_file) and os.path.exists(fy4a_tile_file):
                continue
            xoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * j
            yoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * i
            # the last column                 
            if j == cnum_tile:
                xoff = xsize - BLOCK_SIZE
            # the last row
            if i == rnum_tile:
                yoff = ysize - BLOCK_SIZE
            print("the row and column of tile is :", xoff, yoff)
            
            h08_tile_data = h08ds.ReadAsArray(xoff, yoff, BLOCK_SIZE, BLOCK_SIZE)
            if np.all(h08_tile_data == noDataValue) or np.all(h08_tile_data == 255):
                continue
            fy4a_tile_data = fy4ads.ReadAsArray(xoff, yoff, BLOCK_SIZE, BLOCK_SIZE)
            
            tile_gt[0] = geotrans[0] + xoff * geotrans[1]
            tile_gt[3] = geotrans[3] + yoff * geotrans[5]
            
            dst_format = 'GTiff'
            dst_datatype = dataType
        
            driver = gdal.GetDriverByName(dst_format)
            dst_ds = driver.Create(h08_tile_file, BLOCK_SIZE, BLOCK_SIZE, 1, dst_datatype)
            dst_ds.SetGeoTransform(tile_gt)
            dst_ds.SetProjection(proj)
            dst_ds.GetRasterBand(1).WriteArray(h08_tile_data)
            
            del dst_ds
            
            dst_nbands = len(fy4a_tile_data)
            driver = gdal.GetDriverByName(dst_format)
            dst_ds1 = driver.Create(fy4a_tile_file, BLOCK_SIZE, BLOCK_SIZE, dst_nbands, fydataType)
            dst_ds1.SetGeoTransform(tile_gt)
            dst_ds1.SetProjection(proj)
            if dst_nbands == 1:
    #             tile_data[tile_data == noDataValue] = np.NaN
                dst_ds1.GetRasterBand(1).WriteArray(fy4a_tile_data)
            else:
                for b in range(dst_nbands):
    #                 tile_data[i][tile_data[i] == noDataValue] = np.NaN
                    dst_ds1.GetRasterBand(b + 1).WriteArray(fy4a_tile_data[b])
            del dst_ds1
            
    h08_rm_cmd = 'rm -rf %s' % (h08file)
    os.system(h08_rm_cmd)
    fy4a_rm_cmd = 'rm -rf %s' % (fy4afile)
    os.system(fy4a_rm_cmd)

    
def fy4a2tiles(fy4afile, BLOCK_SIZE, OVERLAP_SIZE):
    print('the image is :', fy4afile)
    dataset = gdal.Open(fy4afile)
    if dataset is None:
        print("Failed to open file: " + fy4afile)
        sys.exit(1)
    band = dataset.GetRasterBand(1)
    xsize = dataset.RasterXSize
    ysize = dataset.RasterYSize
    proj = dataset.GetProjection()
    noDataValue = band.GetNoDataValue()
    dataType = band.DataType
    geotrans = dataset.GetGeoTransform()
    tile_gt = list(geotrans)
    noDataValue = band.GetNoDataValue()
    
    rnum_tile = int((ysize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    cnum_tile = int((xsize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    print('the number of tile is :', rnum_tile * cnum_tile)
    
    for i in range(rnum_tile + 1):
        print(i)
        for j in range(cnum_tile + 1):
            xoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * j
            yoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * i
            # the last column                 
            if j == cnum_tile:
                xoff = xsize - BLOCK_SIZE
            # the last row
            if i == rnum_tile:
                yoff = ysize - BLOCK_SIZE
            print("the row and column of tile is :", xoff, yoff)
            
            tile_data = dataset.ReadAsArray(xoff, yoff, BLOCK_SIZE, BLOCK_SIZE)
            if np.all(tile_data == noDataValue) or np.all(tile_data == 255):
                continue
            
            tile_gt[0] = geotrans[0] + xoff * geotrans[1]
            tile_gt[3] = geotrans[3] + yoff * geotrans[5]
            tile_name = fy4afile.split('/')[-1].split('.')[0] + '_' + str(i) + '_' + str(j) + '.tif'
            tile_file = os.path.join(fy4a_tiles, tile_name)
            dst_format = 'GTiff'
            dst_datatype = dataType
            dst_nbands = len(tile_data)
            
            driver = gdal.GetDriverByName(dst_format)
            dst_ds = driver.Create(tile_file, BLOCK_SIZE, BLOCK_SIZE, dst_nbands, dst_datatype)
            dst_ds.SetGeoTransform(tile_gt)
            dst_ds.SetProjection(proj)
            
            if dst_nbands == 1:
    #             tile_data[tile_data == noDataValue] = np.NaN
                dst_ds.GetRasterBand(1).WriteArray(tile_data)
            else:
                for b in range(dst_nbands):
    #                 tile_data[i][tile_data[i] == noDataValue] = np.NaN
                    dst_ds.GetRasterBand(b + 1).WriteArray(tile_data[b])
            del dst_ds
    rm_cmd = 'rm -rf %s' % (fy4afile)
    os.system(rm_cmd)
            
            
def h082tiles(h08file, BLOCK_SIZE, OVERLAP_SIZE):
    print('the image is :', h08file)
    dataset = gdal.Open(h08file)
    if dataset is None:
        print("Failed to open file: " + h08file)
        sys.exit(1)
    band = dataset.GetRasterBand(1)
    xsize = dataset.RasterXSize
    ysize = dataset.RasterYSize
    proj = dataset.GetProjection()
    noDataValue = band.GetNoDataValue()
    dataType = band.DataType
    geotrans = dataset.GetGeoTransform()
    tile_gt = list(geotrans)
    noDataValue = band.GetNoDataValue()
    
    rnum_tile = int((ysize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    cnum_tile = int((xsize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    print('the number of tile is :', rnum_tile * cnum_tile)
    
    for i in range(rnum_tile + 1):
        print(i)
        for j in range(cnum_tile + 1):
            xoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * j
            yoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * i
            # the last column                 
            if j == cnum_tile:
                xoff = xsize - BLOCK_SIZE
            # the last row
            if i == rnum_tile:
                yoff = ysize - BLOCK_SIZE
            print("the row and column of tile is :", xoff, yoff)
            
            tile_data = band.ReadAsArray(xoff, yoff, BLOCK_SIZE, BLOCK_SIZE)
            if np.all(tile_data == noDataValue) or np.all(tile_data == 255):
                continue
            
            tile_gt[0] = geotrans[0] + xoff * geotrans[1]
            tile_gt[3] = geotrans[3] + yoff * geotrans[5]
            tile_name = h08file.split('/')[-1].split('.')[0] + '_' + str(i) + '_' + str(j) + '.tif'
            print(tile_name)
            tile_file = os.path.join(h08_tiles, tile_name)
            dst_format = 'GTiff'
            dst_datatype = dataType
        
            driver = gdal.GetDriverByName(dst_format)
            dst_ds = driver.Create(tile_file, BLOCK_SIZE, BLOCK_SIZE, 1, dst_datatype)
            dst_ds.SetGeoTransform(tile_gt)
            dst_ds.SetProjection(proj)
            dst_ds.GetRasterBand(1).WriteArray(tile_data)
            
            del dst_ds
    rm_cmd = 'rm -rf %s' % (h08file)
    os.system(rm_cmd)


def epst_to_wkt(in_proj):     
    inSpatialRef = osr.SpatialReference()
    inSpatialRef.SetFromUserInput(in_proj)
    wkt = inSpatialRef.ExportToWkt()
    return wkt


def preprocessh08(h08file, h08_dst):
#     6雨层云Ns
#     9层云St
#     7积云Cu
#     1卷云Ci
#     3深对流云Dc
#     255其他
    wgs84_wkt = epst_to_wkt('EPSG:4326')
    
    ds = nc.Dataset(h08file, 'r')
    CLTYPEarr = ds.variables['CLTYPE'][:].data
    CLTYPEarr[CLTYPEarr==2]=255
    CLTYPEarr[CLTYPEarr==4]=255
    CLTYPEarr[CLTYPEarr==5]=255
    CLTYPEarr[CLTYPEarr==8]=255
    CLTYPEarr[CLTYPEarr==0]=255
    CLTYPEarr[CLTYPEarr==10]=255
    Lon = ds.variables['longitude'][:]
    Lat = ds.variables['latitude'][:]
    LonMin, LatMax, LonMax, LatMin = [Lon.min(), Lat.max(), Lon.max(), Lat.min()] 
    # 分辨率计算
    N_Lat = len(Lat) 
    N_Lon = len(Lon)
    Lon_Res = (LonMax - LonMin) / (float(N_Lon) - 1)
    Lat_Res = (LatMax - LatMin) / (float(N_Lat) - 1)
    # 设置影像的显示范围
    # -Lat_Res一定要是-的
    geotransform = (LonMin, Lon_Res, 0, LatMax, 0, -Lat_Res)
    xsize = N_Lon
    ysize = N_Lat
    nbands = 1
    dtype = NP2GDAL_CONVERSION[str(CLTYPEarr.dtype)]
    h08array_to_tiff(CLTYPEarr, wgs84_wkt, geotransform, 1, '/tmp/test1.tif', dtype)
    clip_cmd = 'gdalwarp -te 79.975  5 140 54 -crop_to_cutline /tmp/test1.tif %s' % (h08_dst)
    rm_cmd = 'rm -rf /tmp/test1.tif'
    os.system(clip_cmd)
    os.system(rm_cmd)


def preprocessfy4a(fy4a_file, fy4a_dst):
    proj = epst_to_wkt('EPSG:4326')
    print('the current fy4afile is %s' % fy4a_file)
    
    data = read_fy4a_arr(fy4a_file)
    xsize = data.shape[1]
    ysize = data.shape[2]
    
    fydataType = NP2GDAL_CONVERSION[str(data.dtype)]
    gt = fy4a_gt 
    dst_nbands = len(CONTENTS['4000M'])
    dst_file = '/tmp/test2.tif'

    dst_format = 'GTiff'
    driver = gdal.GetDriverByName(dst_format)
    dst_ds = driver.Create(dst_file, ysize, xsize, dst_nbands, fydataType)
    dst_ds.SetGeoTransform(gt)
    dst_ds.SetProjection(proj)
    
    if dst_nbands == 1:
        dst_ds.GetRasterBand(1).WriteArray(data)
    else:
        for i in range(dst_nbands):
            dst_ds.GetRasterBand(i + 1).WriteArray(data[i, :, :])
    del dst_ds
    
    clip_cmd = 'gdalwarp -te 79.975  5 140 54 -crop_to_cutline %s %s' % (dst_file,fy4a_dst)
    rm_cmd = 'rm -rf /tmp/test2.tif'
    os.system(clip_cmd)
    os.system(rm_cmd)                       

                
if __name__ == '__main__':
#     fy4a_file = '/data/fy4a/20200101/FY4A-_AGRI--_N_REGC_1047E_L1-_FDI-_MULT_NOM_20200101044918_20200101045335_4000M_V0001.HDF'
#     fy4a_dst = '/tmp/2020fy.tif'
#     preprocessfy4a(fy4a_file, fy4a_dst)
    h08_files = []
    fy4a_files = []
    start = time.time()
    for root, dirs, files in os.walk(h08_root):
        for file in files:
            if file.endswith('.nc'):
                h08date = file.split('_')[2]
                h08time = file.split('_')[3]
                hour = h08time[:2]
                minu = int(h08time[2:])
                if int(hour) >= 1 and int(hour) <= 8:
            #                 daytime in China
                    if minu == 20 or minu == 30 or minu == 40 or minu == 50:  
            #                     have corresponding fy4a data
                        h08_file = os.path.join(root, file)
                        fy4a_file = h8_find_fy4a(h08_file, fy4a_root) 
                        if os.path.exists(fy4a_file):  
            #                             print(h08_file)
#                             h08_dst = os.path.join(h08_tif, 'H08_' + h08date + '_' + h08time + '.tif')
                            h08_dst = os.path.join(h08_tif, h08date + '_' + h08time + '.tif')
#                             fy4a_dst = os.path.join(fy4a_tif, 'FY4A_' + h08date + '_' + h08time + '.tif')
                            fy4a_dst = os.path.join(fy4a_tif, h08date + '_' + h08time + '.tif')
 
            # #                             os.path.join(fy4a_tif,(fy4a_file.split('.')[0]+'.tif').split('/')[-1])            
                            try:
                                preprocessh08(h08_file, h08_dst)
                                preprocessfy4a(fy4a_file, fy4a_dst)
                            except:
                                continue
                            else:
                                h08_files.append(h08_dst)
                                fy4a_files.append(fy4a_dst)
                             
                                tif2tiles(fy4a_dst, h08_dst, BLOCK_SIZE, OVERLAP_SIZE)

#                 h08date = file.split('_')[2]
#                 h08time = file.split('_')[3]
#                 hour = h08time[:2]
#                 minu = int(h08time[2:])
#                 if int(hour)>=1 and int(hour)<=8:
#     #                 daytime in China
#                     if minu ==20 or minu ==30 or minu ==40 or minu ==50:  
#     #                     have corresponding fy4a data
#                         h08_file= os.path.join(root,file)
#                         fy4a_file = h8_find_fy4a(h08_file, fy4a_root) 
#                         if os.path.exists(fy4a_file):  
# #                             print(h08_file)
#                             h08_dst = os.path.join(h08_tif,'H08_'+h08date+'_'+h08time+'.tif')
# # #                             os.path.join(h08_tif,(h08_file.split('.')[0]+'.tif').split('/')[-1])
#                             fy4a_dst = os.path.join(fy4a_tif,'FY4A_'+h08date+'_'+h08time+'.tif')
# #                             os.path.join(fy4a_tif,'FY4A_'+fy4a_file.split('/')[-1][44 : 58]+'.tif')
# # #                             os.path.join(fy4a_tif,(fy4a_file.split('.')[0]+'.tif').split('/')[-1])
#                             preprocessh08(h08_file,h08_dst)
#                             preprocessfy4a(fy4a_file,fy4a_dst)
#                                        
#                             h08_files.append(h08_dst)
#                             fy4a_files.append(fy4a_dst)
#                             
#                             tif2tiles(fy4a_dst, h08_dst, BLOCK_SIZE, OVERLAP_SIZE)
#                             fy4a2tiles(fy4a_dst, BLOCK_SIZE, OVERLAP_SIZE)
#                             h082tiles(h08_dst, BLOCK_SIZE, OVERLAP_SIZE)
#     path_name = {'h08_path': h08_files, 'fy4a_path': fy4a_files}
#     df = pd.DataFrame(path_name)
#     df.to_csv('/data/cloud_files_path.csv', index = False)
    end = time.time()
    period = end - start
    print('the total time is %s seconds' % period)
