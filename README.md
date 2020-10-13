# fy4a_h08_process
根据葵花8号数据，查找对应时间的风云卫星数据，中国区域，选取1-8点，即中国时间上午9点到下午4点的数据
葵花8号是每10分钟一景，使用其20，30，40，50分钟拍摄图像，对应风云获取时间为对应小时的19，30，38，49分钟
葵花8号数据层包括：
Cloud optical thickness, 
Cloud effective radius, 
Cloud top temperature,
Cloud top height,
Cloud types (ISCCP Definition)
其中cloud types为云类别标签，具体如下：
description: 0=Clear, 1=Ci, 2=Cs, 3=Deep convection, 4=Ac, 5=As, 6=Ns, 7=Cu, 8=Sc, 9=St, 10=Unknown, 255=Fill
雨层云Ns 6
层云St 9
积云Cu 7
卷云Ci 1
深对流云Dc 3
卷层云Cs 2
高积云Ac 4
高层云As 5
层积云Sc 8

将葵花8号和风云4A数据在空间上对齐，然后进行切片划分。
