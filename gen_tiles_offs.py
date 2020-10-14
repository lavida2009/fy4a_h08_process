def gen_tiles_offs(xsize, ysize, BLOCK_SIZE,OVERLAP_SIZE):
    xoff_list = []
    yoff_list = []
    
    cnum_tile = int((xsize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    for j in range(cnum_tile + 1):
        xoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * j      
            # the last column                 
        if j == cnum_tile:
            xoff = xsize - BLOCK_SIZE
            # the last row
        xoff_list.append(xoff)
        
    rnum_tile = int((ysize - BLOCK_SIZE) / (BLOCK_SIZE - OVERLAP_SIZE)) + 1
    for i in range(rnum_tile + 1):
        yoff = 0 + (BLOCK_SIZE - OVERLAP_SIZE) * i
        if i == rnum_tile:
            yoff = ysize - BLOCK_SIZE
        yoff_list.append(yoff)
    
    if xoff_list[-1]==xoff_list[-2]:
        xoff_list.pop()
    if yoff_list[-1]==yoff_list[-2]:# the last tile overlap with the last second tile
        yoff_list.pop()

    return xoff_list,yoff_list

if __name__ == '__main__':
    
    xsize = 770
    ysize = 1024
    x,y = gen_tiles_offs(xsize, ysize, 256,150)
    print(x)
    print(y)
    
    