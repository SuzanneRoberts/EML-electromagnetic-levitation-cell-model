# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 10:10:33 2015

@author: Suzanne
"""

import numpy as np, pylab as pl

exp = np.array([
[0.137255, 110.769],
[0.392157, 286.154],
[0.604847, 370.513],
[0.693355, 376.923],
[0.795479, 447.436],
[0.999728, 402.564],
[1.3061,   344.872],
[1.59205,  223.077],
[1.878,    165.385],
[2.27968,  146.154],
[2.50436,  82.0513] ])

exp12mm = np.array([
[0.0923359, 149.926],
[0.322923,  444.674],
[0.512085,  607.386],
[0.610656,  653.749],
[0.807726,  723.155],
[0.889869,  761.791],
[1.20135,   698.721],
[1.29951,   612.933],
[1.49616,   550.188],
[1.79895,   331.669],
[2.11063,   330.788],
[2.39712,   143.409],
[2.68407,   103.73]  ])

exp15mm = np.array([
[0.224425, 421.631],
[0.522,    1151.51],
[0.710965, 1252.04],
[1.01503,  1437.75],
[1.32617,  1265.85],
[1.6044,   1062.94],
[1.91545,  859.949],
[1.96441,  782.074],
[2.39802,  431.033],
[2.783,    266.697],  ])

exp20mm = np.array([
[0.214064, 1445.75],
[0.394114, 2134.76],
[0.516799, 2738.59],
[0.623306, 2947.64],
[0.787241, 3133.5],
[0.967585, 3311.62],
[1.3042,   2738.99],
[1.62453,  1965.09],
[1.87092,  1392.41],
[2.22395,  804.306],
[2.4948,   487.079],  ])


zvecRoyer = np.array([
   -5.0000,
   -4.5455,
   -4.0909,
   -3.6364,
   -3.1818,
   -2.7273,
   -2.2727,
   -1.8182,
   -1.3636,
   -0.9091,
   -0.4545,
         0,
    0.4545,
    0.9091,
    1.3636,
    1.8182,
    2.2727,
    2.7273,
    3.1818,
    3.6364,
    4.0909,
    4.5455,
    5.0000,
    5.4545,
    5.9091,
    6.3636,
    6.8182,
    7.2727,
    7.7273,
    8.1818,
    8.6364,
    9.0909,
    9.5455,
   10.0000,
   10.4545,
   10.9091,
   11.3636,
   11.8182,
   12.2727,
   12.7273,
   13.1818,
   13.6364,
   14.0909,
   14.5455,
   15.0000,
   15.4545,
   15.9091,
   16.3636,
   16.8182,
   17.2727,
   17.7273,
   18.1818,
   18.6364,
   19.0909,
   19.5455,
   20.0000,
   20.4545,
   20.9091,
   21.3636,
   21.8182,
   22.2727,
   22.7273,
   23.1818,
   23.6364,
   24.0909,
   24.5455,
   25.0000,
   25.4545,
   25.9091,
   26.3636,
   26.8182,
   27.2727,
   27.7273,
   28.1818,
   28.6364,
   29.0909,
   29.5455,
   30.0000,
   30.4545,
   30.9091,
   31.3636,
   31.8182,
   32.2727,
   32.7273,
   33.1818,
   33.6364,
   34.0909,
   34.5455,
   35.0000,
   35.4545,
   35.9091,
   36.3636,
   36.8182,
   37.2727,
   37.7273,
   38.1818,
   38.6364,
   39.0909,
   39.5455,
   40.0000 ])
   
Fz_curveRoyer = np.array([
   -0.3743,
   -0.3533,
   -0.3290,
   -0.3016,
   -0.2712,
   -0.2381,
   -0.2025,
   -0.1647,
   -0.1252,
   -0.0842,
   -0.0424,
         0,
    0.0424,
    0.0842,
    0.1252,
    0.1647,
    0.2025,
    0.2381,
    0.2712,
    0.3016,
    0.3290,
    0.3533,
    0.3743,
    0.3921,
    0.4067,
    0.4180,
    0.4262,
    0.4315,
    0.4340,
    0.4340,
    0.4315,
    0.4269,
    0.4204,
    0.4122,
    0.4025,
    0.3916,
    0.3797,
    0.3670,
    0.3537,
    0.3399,
    0.3258,
    0.3116,
    0.2974,
    0.2832,
    0.2692,
    0.2554,
    0.2420,
    0.2289,
    0.2163,
    0.2041,
    0.1924,
    0.1811,
    0.1704,
    0.1601,
    0.1504,
    0.1412,
    0.1324,
    0.1241,
    0.1163,
    0.1089,
    0.1020,
    0.0955,
    0.0893,
    0.0836,
    0.0782,
    0.0731,
    0.0684,
    0.0639,
    0.0598,
    0.0559,
    0.0523,
    0.0489,
    0.0457,
    0.0427,
    0.0400,
    0.0374,
    0.0350,
    0.0328,
    0.0307,
    0.0287,
    0.0269,
    0.0252,
    0.0236,
    0.0221,
    0.0207,
    0.0194,
    0.0182,
    0.0171,
    0.0160,
    0.0151,
    0.0141,
    0.0133,
    0.0125,
    0.0117,
    0.0110,
    0.0104,
    0.0098,
    0.0092,
    0.0086,
    0.0081 ])
    
Fz_curveRoyer12mm = np.array([
   -0.6511,
   -0.6145,
   -0.5722,
   -0.5246,
   -0.4717,
   -0.4141,
   -0.3522,
   -0.2865,
   -0.2177,
   -0.1465,
   -0.0737,
         0,
    0.0737,
    0.1465,
    0.2177,
    0.2865,
    0.3522,
    0.4141,
    0.4717,
    0.5246,
    0.5722,
    0.6145,
    0.6511,
    0.6820,
    0.7073,
    0.7270,
    0.7414,
    0.7506,
    0.7549,
    0.7548,
    0.7505,
    0.7425,
    0.7311,
    0.7168,
    0.7000,
    0.6811,
    0.6604,
    0.6383,
    0.6152,
    0.5912,
    0.5667,
    0.5420,
    0.5172,
    0.4925,
    0.4682,
    0.4443,
    0.4209,
    0.3982,
    0.3762,
    0.3550,
    0.3346,
    0.3150,
    0.2963,
    0.2785,
    0.2616,
    0.2455,
    0.2303,
    0.2159,
    0.2023,
    0.1894,
    0.1774,
    0.1660,
    0.1553,
    0.1453,
    0.1359,
    0.1271,
    0.1189,
    0.1112,
    0.1039,
    0.0972,
    0.0909,
    0.0850,
    0.0795,
    0.0743,
    0.0695,
    0.0651,
    0.0609,
    0.0570,
    0.0533,
    0.0499,
    0.0467,
    0.0438,
    0.0410,
    0.0384,
    0.0360,
    0.0338,
    0.0317,
    0.0297,
    0.0279,
    0.0262,
    0.0246,
    0.0231,
    0.0217,
    0.0204,
    0.0192,
    0.0180,
    0.0170,
    0.0160,
    0.0150,
    0.0142 ])    
    
Fz_curveRoyer15mm = np.array([
   -1.2799,
   -1.2079,
   -1.1249,
   -1.0312,
   -0.9273,
   -0.8141,
   -0.6923,
   -0.5632,
   -0.4279,
   -0.2880,
   -0.1448,
         0,
    0.1448,
    0.2880,
    0.4279,
    0.5632,
    0.6923,
    0.8141,
    0.9273,
    1.0312,
    1.1249,
    1.2079,
    1.2799,
    1.3407,
    1.3904,
    1.4292,
    1.4573,
    1.4754,
    1.4840,
    1.4837,
    1.4753,
    1.4595,
    1.4372,
    1.4092,
    1.3761,
    1.3389,
    1.2983,
    1.2548,
    1.2093,
    1.1621,
    1.1140,
    1.0654,
    1.0167,
    0.9682,
    0.9203,
    0.8733,
    0.8274,
    0.7827,
    0.7395,
    0.6978,
    0.6577,
    0.6192,
    0.5825,
    0.5475,
    0.5142,
    0.4826,
    0.4527,
    0.4244,
    0.3976,
    0.3724,
    0.3487,
    0.3264,
    0.3054,
    0.2857,
    0.2672,
    0.2499,
    0.2337,
    0.2185,
    0.2043,
    0.1911,
    0.1787,
    0.1671,
    0.1562,
    0.1461,
    0.1367,
    0.1279,
    0.1197,
    0.1120,
    0.1048,
    0.0981,
    0.0919,
    0.0861,
    0.0806,
    0.0756,
    0.0708,
    0.0664,
    0.0623,
    0.0584,
    0.0548,
    0.0515,
    0.0483,
    0.0454,
    0.0427,
    0.0401,
    0.0377,
    0.0355,
    0.0334,
    0.0314,
    0.0295,
    0.0278  ])
    
Fz_curveRoyer20mm = np.array([
   -3.0533,
   -2.8816,
   -2.6835,
   -2.4600,
   -2.2123,
   -1.9421,
   -1.6516,
   -1.3436,
   -1.0209,
   -0.6870,
   -0.3454,
         0,
    0.3454,
    0.6870,
    1.0209,
    1.3436,
    1.6516,
    1.9421,
    2.2123,
    2.4600,
    2.6835,
    2.8816,
    3.0533,
    3.1984,
    3.3170,
    3.4094,
    3.4767,
    3.5198,
    3.5402,
    3.5395,
    3.5194,
    3.4818,
    3.4286,
    3.3617,
    3.2829,
    3.1942,
    3.0972,
    2.9935,
    2.8848,
    2.7724,
    2.6577,
    2.5416,
    2.4254,
    2.3098,
    2.1955,
    2.0834,
    1.9738,
    1.8673,
    1.7641,
    1.6646,
    1.5690,
    1.4773,
    1.3897,
    1.3061,
    1.2267,
    1.1513,
    1.0799,
    1.0123,
    0.9486,
    0.8884,
    0.8318,
    0.7786,
    0.7285,
    0.6815,
    0.6375,
    0.5962,
    0.5575,
    0.5213,
    0.4875,
    0.4558,
    0.4262,
    0.3986,
    0.3727,
    0.3486,
    0.3261,
    0.3051,
    0.2855,
    0.2671,
    0.2501,
    0.2341,
    0.2192,
    0.2053,
    0.1924,
    0.1803,
    0.1690,
    0.1585,
    0.1486,
    0.1394,
    0.1308,
    0.1228,
    0.1153,
    0.1083,
    0.1018,
    0.0957,
    0.0899,
    0.0846,
    0.0796,
    0.0749,
    0.0705,
    0.0664])
    
posVecSz = np.array([-0.05      , -0.04473684, -0.03947368, -0.03421053, -0.02894737,
       -0.02368421, -0.01842105, -0.01315789, -0.00789474, -0.00263158,
        0.00263158,  0.00789474,  0.01315789,  0.01842105,  0.02368421,
        0.02894737,  0.03421053,  0.03947368,  0.04473684,  0.05      ])
        
forceVecSz = np.array([[ -1.44469936e-05],
       [ -2.27309575e-05],
       [ -3.69335564e-05],
       [ -6.19909188e-05],
       [ -1.07101181e-04],
       [ -1.88120295e-04],
       [ -3.25174792e-04],
       [ -5.11007688e-04],
       [ -6.04291985e-04],
       [ -3.03236504e-04],
       [  3.03236504e-04],
       [  6.04291985e-04],
       [  5.11007688e-04],
       [  3.25174792e-04],
       [  1.88120295e-04],
       [  1.07101181e-04],
       [  6.19909188e-05],
       [  3.69335564e-05],
       [  2.27309575e-05],
       [  1.44469936e-05]])
       
posVecSz25 = np.array([-0.05      , -0.04473684, -0.03947368, -0.03421053, -0.02894737,
       -0.02368421, -0.01842105, -0.01315789, -0.00789474, -0.00263158,
        0.00263158,  0.00789474,  0.01315789,  0.01842105,  0.02368421,
        0.02894737,  0.03421053,  0.03947368,  0.04473684,  0.05      ])
        
forceVecSz25 = np.array([[  1.94582524e-06],
       [  3.29543446e-06],
       [  5.84081112e-06],
       [  1.08470869e-05],
       [  2.10008874e-05],
       [  4.15713493e-05],
       [  8.00326934e-05],
       [  1.34427299e-04],
       [  1.59355745e-04],
       [  7.73046342e-05],
       [ -7.73046342e-05],
       [ -1.59355745e-04],
       [ -1.34427299e-04],
       [ -8.00326934e-05],
       [ -4.15713493e-05],
       [ -2.10008874e-05],
       [ -1.08470869e-05],
       [ -5.84081112e-06],
       [ -3.29543446e-06],
       [ -1.94582524e-06]])

pl.figure()
pl.plot(exp[:,0]*10.0, exp[:,1]/1000.0, 'ro', markersize=10)
pl.plot(zvecRoyer, Fz_curveRoyer, linewidth=2)
#pl.plot(exp12mm[:,0]*10.0, exp12mm[:,1]/1000.0, 'go')
#pl.plot(zvecRoyer, Fz_curveRoyer12mm)
#pl.plot(exp15mm[:,0]*10.0, exp15mm[:,1]/1000.0, 'ro')
#pl.plot(zvecRoyer, Fz_curveRoyer15mm)
#pl.plot(exp20mm[:,0]*10.0, exp20mm[:,1]/1000.0, 'co')
#pl.plot(zvecRoyer, Fz_curveRoyer20mm)
#pl.plot(posVecSz*1000, forceVecSz/9.81*1000.0) # times 7 fits
#pl.plot(posVecSz25*1000, forceVecSz25/9.81*1000.0)
pl.xlabel('Height [mm]', fontsize=18)
pl.ylabel('Force equivalent [g]', fontsize=18)
#pl.title('Different model predictions for the Fromm & Jehn experiment')
pl.title('Comparison of the model with experimental data', fontsize=18)
pl.legend(['Fromm & Jehn (1965) experimental data',
           'Implementation of Royer et al. (2013) model',
           'Fromm & Jehn (1965) experimental data, dia=12mm',
           'Implementation of Royer et al. (2013) model, dia=12mm',
           'Fromm & Jehn (1965) experimental data, dia=15mm',
           'Implementation of Royer et al. (2013) model, dia=15mm',
           'Fromm & Jehn (1965) experimental data, dia=20mm',
           'Implementation of Royer et al. (2013) model, dia=20mm',
           'El-Kaddah & Szekely model, dia=10mm [20,20,20]',
           'El-Kaddah & Szekely model, dia=10mm [25,25,25]'], loc='lower right')

pl.show()