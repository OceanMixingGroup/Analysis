% Second order polynom fits better that first order polynom in all cases,
% however calibration coefficients for different CTD casts differ more for 
% second order polynom
% MHC should be calibrated with first order polynom
% MHT & T should be calibrated with second order polynom
% First one might want to modify the header coefficients
    head.coef.S1(1)=1e-4*head.coef.S1(1);
    head.coef.S2(1)=1e-4*head.coef.S2(1);
if q.script.num>545 && q.script.num<912 && ~isempty(strfind(head.instrument,'CHAM04-03'))
    head.sensor_id(7,:)='114         ';
    head.sensor_id(8,:)='114         ';
end
if q.script.num>=1631 && ~isempty(strfind(head.instrument,'CHAM04-03'))
    head.sensor_id(7,:)='04_062      ';
    head.sensor_id(8,:)='04_062      ';
end
if ~isempty(strfind(head.instrument,'CHAM04-01'))
    head.coef.TP(2)=0.094741;
elseif ~isempty(strfind(head.instrument,'CHAM04-02'))
    head.coef.TP(2)=0.086116;
elseif ~isempty(strfind(head.instrument,'CHAM04-03'))
    head.coef.TP(2)=0.0995492;
elseif ~isempty(strfind(head.instrument,'CHAM04-04'))
    head.coef.TP(2)=0.083883;
end
%% T calibration coefficients
if ~isempty(strfind(head.sensor_id(5,:),'91_01')) && ~isempty(strfind(head.instrument,'CHAM04-01'))% CHAM 04-01
    % CTD #2 & cham#6 - 2 profiles
    head.coef.T=[16.3928 3.3474 0.0280 0 1];
elseif ~isempty(strfind(head.sensor_id(5,:),'90_03')) && ~isempty(strfind(head.instrument,'CHAM04-04'))% CHAM 04-04
    % casts 16-66,358-541
    if q.script.num<=66
        % CTD #9 & cham#17
        head.coef.T=[14.8539 3.4715 0.0229 0 1];
    elseif q.script.num>=358 && q.script.num<=446
        % CTD #12 & cham#358
        head.coef.T=[15.0758    3.0871    0.1393 0 1]; 
    else
        % CTD #13 & cham#541
        head.coef.T=[15.0434 2.9989 0.1761 0 1]; % fits well with a second order polynom
    end
elseif ~isempty(strfind(head.sensor_id(5,:),'91_60')) && ~isempty(strfind(head.instrument,'CHAM04-04'))% CHAM 04-04
%     in water casts 1851 - 1866
    % CTD #23 & cham#1853
    head.coef.T=[15.4768 3.7544 0.0630 0 1]; % fits well with a second order polynom
elseif ~isempty(strfind(head.sensor_id(5,:),'02_21')) && ~isempty(strfind(head.instrument,'CHAM04-02'))% CHAM 04-02 
    % CTD #3 & cham#9 (CTD is 30 min later than Cham)
    head.coef.T=[10.5179 5.7503 -0.3640 0 1];
elseif ~isempty(strfind(head.sensor_id(5,:),'91_12')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num>=67 && q.script.num<=357
        head.coef.T=[15.8715    3.5692    0.0703 0 1];% CTD 14 & cham#548
    elseif q.script.num>=542 && q.script.num<=545
        head.coef.T=[15.6681    3.9531   -0.0534 0 1];% CTD 13 & cham#543
    elseif q.script.num>=546 && q.script.num<=550
        head.coef.T=[16.0938    2.9980    0.2606 0 1];% CTD 14 & cham#548
    elseif q.script.num==551
        head.coef.T=[15.6282    4.1564   -0.1286 0 1];% CTD 14 & cham#551
    elseif q.script.num==552
        head.coef.T=[15.9665    3.4285    0.1014 0 1];% CTD 14 & cham#552
    elseif q.script.num==553
        head.coef.T=[15.9330    3.5453    0.0638 0 1];% CTD 14 & cham#553
    elseif q.script.num>=554 && q.script.num<=796
        head.coef.T=[15.8486    3.6907    0.0222 0 1];% CTD 14 & cham#554
    elseif q.script.num>=797 && q.script.num<=1629
        head.coef.T=[15.8832 3.6521 0.0262 0 1];% CTD 14 & cham#554
    elseif q.script.num>=1631 && q.script.num<=1654
        head.coef.T=[15.9436 3.6941 -0.0116 0 1];% CTD #21 & cham#1624
    elseif q.script.num>=1655 && q.script.num<=1866
        head.coef.T=[15.8648    3.6527    0.0294 0 1];% CTD #22 & cham#1806
    elseif q.script.num>=1867 && q.script.num<=2668
        head.coef.T=[15.8577 3.5722 0.0656 0 1];% CTD 19 & cham#1416
    end
end
%% MHT calibration coefficients 
if ~isempty(strfind(head.sensor_id(8,:),'04_059')) && ~isempty(strfind(head.instrument,'CHAM04-01'))% CHAM 04-01
    % CTD #1 & cham#4. Not very good calibration
    head.coef.MHT=[-6.9594 -24.2379 -7.6707 -1.0119 1];
elseif ~isempty(strfind(head.sensor_id(8,:),'114')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==67
        head.coef.MHT=[9.6953   -3.4632    0.3580 0 1];
    elseif q.script.num==68
        head.coef.MHT=[9.3950   -3.8141    0.2830 0 1];
    elseif q.script.num==69
        head.coef.MHT=[9.3442   -3.8771    0.2743 0 1];
    elseif q.script.num==70
        head.coef.MHT=[9.4490   -3.6371    0.3375 0 1];
    elseif q.script.num==71
        head.coef.MHT=[8.8125   -4.2825    0.2271 0 1];
    elseif q.script.num==72
        head.coef.MHT=[8.4406   -4.7053    0.1577 0 1];
    elseif q.script.num==73
        head.coef.MHT=[8.6810   -4.2970    0.2668 0 1];
    elseif q.script.num==74
        head.coef.MHT=[8.6810   -4.2970    0.2668 0 1];
    elseif q.script.num==75
        head.coef.MHT=[9.0981   -3.3246    0.8689 0 1];
    elseif q.script.num==76
        head.coef.MHT=[9.0981   -3.3246    0.8689 0 1];
    elseif q.script.num>=546 && q.script.num<=547%!!! board has been changed before drop 546
         head.coef.MHT=[11.5756   -3.1038    0.0794 0 1];% BAD!!!
    elseif q.script.num>=548  && q.script.num<=552
         head.coef.MHT=[11.8766   -3.2092    0.0278 0 1];% BAD!!!
    elseif q.script.num>=548  && q.script.num<=551
         head.coef.MHT=[11.8766   -3.2092    0.0278 0 1];% BAD!!!
    elseif q.script.num==552
         head.coef.MHT=[12.3277   -2.6563    0.1453 0 1];
    elseif q.script.num==553
         head.coef.MHT=[12.1646   -2.7822    0.1150 0 1];
    elseif q.script.num>=554  && q.script.num<=558
         head.coef.MHT=[11.6900   -3.1609    0.0511 0 1];
    elseif q.script.num>=559  && q.script.num<=563
         head.coef.MHT=[11.6681   -2.9069    0.1337 0 1];
    elseif q.script.num>=564  && q.script.num<=575
         head.coef.MHT=[11.6681   -2.9069    0.1337 0 1];
    elseif q.script.num>=576  && q.script.num<=585
         head.coef.MHT=[11.6870   -2.8730    0.1648 0 1];
    elseif q.script.num>=586  && q.script.num<=595
         head.coef.MHT=[11.4088   -3.1388    0.1056 0 1];
    elseif q.script.num>=596  && q.script.num<=615
         head.coef.MHT=[11.5360   -2.9549    0.1344 0 1];
    elseif q.script.num>=616  && q.script.num<=635
         head.coef.MHT=[11.4640   -3.0186    0.0990 0 1];
    elseif q.script.num>=636  && q.script.num<=655
         head.coef.MHT=[11.5058   -2.8947    0.1142 0 1];
    elseif q.script.num>=656  && q.script.num<=675
         head.coef.MHT=[11.3966   -3.0531    0.0749 0 1];
    elseif q.script.num>=676  && q.script.num<=695
         head.coef.MHT=[11.4202   -2.9962    0.0845 0 1];
    elseif q.script.num>=696  && q.script.num<=715
         head.coef.MHT=[11.4717   -2.9107    0.1032 0 1];
    elseif q.script.num>=716  && q.script.num<=739
         head.coef.MHT=[11.5160   -2.8491    0.1132 0 1];
    elseif q.script.num>=740  && q.script.num<=796
         head.coef.MHT=[11.2778 -3.0973 0.0658 0 1];% CTD #15 & cham#740
    elseif q.script.num==797
         head.coef.MHT=[10.9202   -3.0776    0.0725 0 1];% cham was out of water
    elseif q.script.num==798
         head.coef.MHT=[11.4169   -2.9170    0.1000 0 1];
    elseif q.script.num>=799  && q.script.num<=801
         head.coef.MHT=[11.7671   -2.3936    0.2124 0 1];
    elseif q.script.num>=802  && q.script.num<=806
         head.coef.MHT=[11.5929   -2.6316    0.1634 0 1];
    elseif q.script.num>=807  && q.script.num<=811
         head.coef.MHT=[11.4258   -2.9134    0.1014 0 1];
    elseif q.script.num>=812  && q.script.num<=840
         head.coef.MHT=[11.4549   -2.9275    0.0922 0 1];
    elseif q.script.num>=841  && q.script.num<=911
         head.coef.MHT=[11.5916   -2.8090    0.1178 0 1];
    elseif q.script.num==912 % cham was out of water
         head.coef.MHT=[11.9138   -2.4834    0.1911 0 1];
    elseif q.script.num==913
         head.coef.MHT=[11.6361   -3.0023    0.0664 0 1];
    elseif q.script.num==914
         head.coef.MHT=[11.6406   -2.9539    0.0781 0 1];
    elseif q.script.num>=915  && q.script.num<=917
         head.coef.MHT=[11.6677   -2.8744    0.0977 0 1];
    elseif q.script.num>=918  && q.script.num<=922
         head.coef.MHT=[11.6387   -2.8595    0.1043 0 1];
    elseif q.script.num>=923  && q.script.num<=1020
         head.coef.MHT=[11.6229   -2.9435    0.0803 0 1];
    elseif q.script.num>=1021  && q.script.num<=1040
         head.coef.MHT=[11.6409   -2.9539    0.0760 0 1];
    elseif q.script.num>=1041  && q.script.num<=1079
         head.coef.MHT=[11.6128   -2.9915    0.0681 0 1];
    elseif q.script.num>=1080  && q.script.num<=1120
         head.coef.MHT=[11.6516   -2.9419    0.0785 0 1];
    elseif q.script.num>=1121  && q.script.num<=1131
         head.coef.MHT=[11.6636   -2.8795    0.0931 0 1];
    elseif q.script.num>=1131  && q.script.num<=1141
         head.coef.MHT=[11.6271   -2.9771    0.0708 0 1];
    elseif q.script.num==1142 % cham was out of water
         head.coef.MHT=[11.6599   -2.9514    0.0750 0 1];
    elseif q.script.num>=1143  && q.script.num<=1145
         head.coef.MHT=[11.6599   -2.9514    0.0750 0 1];
    elseif q.script.num>=1146  && q.script.num<=1150
         head.coef.MHT=[11.6550   -2.9233    0.0814 0 1];
    elseif q.script.num>=1151  && q.script.num<=1160
         head.coef.MHT=[11.6363   -2.9567    0.0750 0 1];
    elseif q.script.num>=1161  && q.script.num<=1180
         head.coef.MHT=[11.6456   -2.9361    0.0800 0 1];
    elseif q.script.num>=1181  && q.script.num<=1200
         head.coef.MHT=[11.5154   -2.9245    0.0941 0 1];
    elseif q.script.num>=1201  && q.script.num<=1230
         head.coef.MHT=[11.5027   -2.9709    0.0810 0 1];
    elseif q.script.num>=1231  && q.script.num<=1250
         head.coef.MHT=[11.5040   -2.9844    0.0784 0 1];
    elseif q.script.num>=1251  && q.script.num<=1270
         head.coef.MHT=[11.6215   -2.8695    0.0995 0 1];
    elseif q.script.num>=1271  && q.script.num<=1278
         head.coef.MHT=[11.6110   -2.3532    0.2395 0 1];
    end
elseif ~isempty(strfind(head.sensor_id(8,:),'113')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==77
        head.coef.MHT=[10.9709   -4.6259    0.1242 0 1];
    elseif q.script.num==78
        head.coef.MHT=[10.9378   -4.6833    0.1070 0 1];
    elseif q.script.num>=79 && q.script.num<=81
        head.coef.MHT=[10.9626   -4.6481    0.1171 0 1];
    elseif q.script.num>=82 && q.script.num<=87
        head.coef.MHT=[10.9554   -4.6637    0.1124 0 1];
    elseif q.script.num>=88 && q.script.num<=94
        head.coef.MHT=[10.9166   -4.7195    0.0973 0 1];
    elseif q.script.num>=95 && q.script.num<=104
        head.coef.MHT=[10.9527   -4.6768    0.1077 0 1];
    elseif q.script.num>=105 && q.script.num<=114
        head.coef.MHT=[10.9095   -4.7183    0.0992 0 1];
    elseif q.script.num>=115 && q.script.num<=134
        head.coef.MHT=[10.9202   -4.7314    0.0911 0 1];
    elseif q.script.num>=135 && q.script.num<=154
        head.coef.MHT=[10.9381   -4.6878    0.1063 0 1];
    elseif q.script.num>=155 && q.script.num<=164
        head.coef.MHT=[10.9389   -4.6881    0.1050 0 1];
    elseif q.script.num>=165 && q.script.num<=174
        head.coef.MHT=[10.9372   -4.6818    0.1088 0 1];
    elseif q.script.num>=175 && q.script.num<=196
        head.coef.MHT=[10.9341   -4.6774    0.1108 0 1];
    elseif q.script.num==197
        head.coef.MHT=[11.0049   -4.5868    0.1328 0 1];
    elseif q.script.num==198
        head.coef.MHT=[10.9725   -4.6146    0.1303 0 1];
    elseif q.script.num==199
        % CTD #11 & cham#199 - long time in the water
        head.coef.MHT=[10.9349   -4.7028    0.1021 0 1];% first order is a bit worse than second order
    elseif q.script.num==200
        % CTD #11 & cham#200 - first profile
        head.coef.MHT=[10.9349   -4.7028    0.1021 0 1];% first order is a bit worse than second order
    elseif q.script.num==201
        head.coef.MHT=[10.9574   -4.6481    0.1194 0 1];
    elseif q.script.num>=202 && q.script.num<=204
        head.coef.MHT=[10.9366   -4.6906    0.1051 0 1];
    elseif q.script.num>=205 && q.script.num<=207
        head.coef.MHT=[10.9574   -4.6481    0.1194 0 1];
    elseif q.script.num>=208 && q.script.num<=210
        head.coef.MHT=[10.9384   -4.6764    0.1120 0 1];
    elseif q.script.num>=211 && q.script.num<=215
        head.coef.MHT=[10.9429   -4.6805    0.1081 0 1];
    elseif q.script.num>=216 && q.script.num<=225
        head.coef.MHT=[10.9214   -4.7275    0.0933 0 1];
    elseif q.script.num>=226 && q.script.num<=235
        head.coef.MHT=[10.9387   -4.6746    0.1129 0 1];
    elseif q.script.num>=236 && q.script.num<=285
        head.coef.MHT=[10.9384   -4.6753    0.1112 0 1];
    elseif q.script.num>=286 && q.script.num<=305
        head.coef.MHT=[10.9500   -4.6745    0.1116 0 1];
    elseif q.script.num>=306 && q.script.num<=319
        head.coef.MHT=[10.9393   -4.6826    0.1077 0 1];
    elseif q.script.num>=320 && q.script.num<=339
        head.coef.MHT=[10.9893   -4.5816    0.1380 0 1];
    elseif q.script.num>=340 && q.script.num<=356
        head.coef.MHT=[11.0161   -4.7584    0.0874 0 1];
    elseif q.script.num==357
        % CTD #12 & cham#357 - long time in the water
        head.coef.MHT=[11.0652   -4.6890    0.1079 0 1]; % fits well with a second order polynom
    end
elseif ~isempty(strfind(head.sensor_id(8,:),'04_059')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==542
        head.coef.MHT=[11.8153   -3.5151   -0.0052 0 1];% CTD #13 & cham#542
    elseif q.script.num>=543 && q.script.num<=545
        head.coef.MHT=[11.8287   -3.4481    0.0148 0 1];% CTD #13 & cham#544
    end
elseif ~isempty(strfind(head.sensor_id(8,:),'04_062')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==1279 % cham was out of water
         head.coef.MHT=[12.1418   -2.9989    0.0598 0 1];
    elseif q.script.num==1280 
         head.coef.MHT=[12.1273   -3.0327    0.0531 0 1];
    elseif q.script.num>=1281 && q.script.num<=1287
         head.coef.MHT=[12.1331   -3.0118    0.0588 0 1];
    elseif q.script.num>=1288 && q.script.num<=1293
         head.coef.MHT=[12.1355   -3.0151    0.0579 0 1];
    elseif q.script.num>=1294 && q.script.num<=1310
         head.coef.MHT=[12.1387   -3.0052    0.0605 0 1];
    elseif q.script.num>=1311 && q.script.num<=1330
         head.coef.MHT=[12.2151   -2.7975    0.1098 0 1];
    elseif q.script.num>=1331 && q.script.num<=1360
         head.coef.MHT=[12.1558   -2.9628    0.0697 0 1];
    elseif q.script.num>=1361 && q.script.num<=1390
         head.coef.MHT=[12.1695   -2.9342    0.0778 0 1];
    elseif q.script.num>=1391 && q.script.num<=1416
         head.coef.MHT=[12.1382   -3.0099    0.0588 0 1];
    elseif q.script.num==1417
         head.coef.MHT=[12.1691   -2.9462    0.0739 0 1];
    elseif q.script.num==1418
         head.coef.MHT=[12.1410   -3.0143    0.0567 0 1];
    elseif q.script.num==1419
         head.coef.MHT=[12.1593   -2.9599    0.0692 0 1];
    elseif q.script.num>=1420 && q.script.num<=1424
         head.coef.MHT=[12.1479   -3.0159    0.0570 0 1];
    elseif q.script.num>=1425 && q.script.num<=1430
         head.coef.MHT=[12.1517   -2.9712    0.0669 0 1];
    elseif q.script.num>=1431 && q.script.num<=1440
         head.coef.MHT=[12.1143   -3.0201    0.0566 0 1];
    elseif q.script.num>=1441 && q.script.num<=1450
         head.coef.MHT=[12.1119   -3.0102    0.0591 0 1];
    elseif q.script.num>=1451 && q.script.num<=1463
         head.coef.MHT=[12.1653   -2.9506    0.0705 0 1];
    elseif q.script.num>=1464 && q.script.num<=1475
         head.coef.MHT=[12.1294   -3.0293    0.0528 0 1];
    elseif q.script.num>=1476 && q.script.num<=1506
         head.coef.MHT=[12.1249   -3.0184    0.0565 0 1];
    elseif q.script.num==1507
         head.coef.MHT=[12.0890   -3.0327    0.0538 0 1];
    elseif q.script.num==1508
         head.coef.MHT=[12.1314   -2.9335    0.0785 0 1];
    elseif q.script.num==1509
         head.coef.MHT=[12.1250   -2.9429    0.0754 0 1];
    elseif q.script.num==1510
         head.coef.MHT=[12.1141   -3.0386    0.0505 0 1];
    elseif q.script.num>=1511 && q.script.num<=1513
         head.coef.MHT=[12.0902   -3.0193    0.0566 0 1];
    elseif q.script.num>=1514 && q.script.num<=1518
         head.coef.MHT=[12.0902   -3.0193    0.0566 0 1];
    elseif q.script.num>=1519 && q.script.num<=1530
         head.coef.MHT=[12.1032   -3.0182    0.0562 0 1];
    elseif q.script.num>=1531 && q.script.num<=1540
         head.coef.MHT=[12.1230   -2.9746    0.0673 0 1];
    elseif q.script.num>=1541 && q.script.num<=1560
         head.coef.MHT=[12.1032   -3.0202    0.0562 0 1];
    elseif q.script.num>=1561 && q.script.num<=1580
         head.coef.MHT=[12.0903   -3.0439    0.0514 0 1];
    elseif q.script.num>=1581 && q.script.num<=1600
         head.coef.MHT=[12.1014   -3.0147    0.0579 0 1];
    elseif q.script.num>=1601 && q.script.num<=1624
         head.coef.MHT=[12.1318   -3.0061    0.0573 0 1];
    elseif q.script.num>=1631 && q.script.num<=1632
         head.coef.MHT=[12.0923   -3.1048    0.0308 0 1];
    elseif q.script.num>=1633 && q.script.num<=1640
         head.coef.MHT=[12.0704   -3.1229    0.0327 0 1];
    elseif q.script.num>=1641 && q.script.num<=1654
         head.coef.MHT=[12.1073   -3.0863    0.0367 0 1];
    elseif q.script.num>=1655 && q.script.num<=1662
         head.coef.MHT=[12.1121   -3.0090    0.0597 0 1];
    elseif q.script.num>=1663 && q.script.num<=1850
         head.coef.MHT=[12.1298   -2.9862    0.0635 0 1];
    elseif q.script.num==1867
         head.coef.MHT=[12.1334   -2.9781    0.0654 0 1];
    elseif q.script.num==1868
         head.coef.MHT=[12.0761   -3.0200    0.0579 0 1];
    elseif q.script.num>=1869 && q.script.num<=1874
         head.coef.MHT=[12.0932   -3.0267    0.0561 0 1];
    elseif q.script.num>=1875 && q.script.num<=1880
         head.coef.MHT=[12.1104   -3.0238    0.0555 0 1];
    elseif q.script.num>=1881 && q.script.num<=1910
         head.coef.MHT=[12.1418   -2.9619    0.0698 0 1];
    elseif q.script.num>=1911 && q.script.num<=1950
         head.coef.MHT=[12.1343   -2.9670    0.0683 0 1];
    elseif q.script.num>=1951 && q.script.num<=2250
         head.coef.MHT=[12.1343   -2.9670    0.0683 0 1];
    elseif q.script.num>=2251 && q.script.num<=2464
         head.coef.MHT=[12.1925   -2.8890    0.0857 0 1];
    elseif q.script.num>=2465 && q.script.num<=2668
         head.coef.MHT=[12.1940   -2.8508    0.0950 0 1];
    end
end
%% MHC calibration coefficients
if ~isempty(strfind(head.sensor_id(7,:),'08-01')) && ~isempty(strfind(head.instrument,'CHAM04-01'))% CHAM 04-01
    % CTD #2 & cham#6 - 2 profiles
    head.coef.MHC=[2.5169 0.6719 0 0 1]; 
elseif ~isempty(strfind(head.sensor_id(7,:),'08-04'))% CHAM 04-02
%     % CTD #9 & cham#13
%     head.coef.MHC=[33.2497 -25.4256 7.3760 -0.6907 1];
elseif ~isempty(strfind(head.sensor_id(7,:),'08-14')) && ~isempty(strfind(head.instrument,'CHAM04-04'))% CHAM 04-04
    % casts 16-66,358-541
    if q.script.num==16
        head.coef.MHC=[2.5588    0.6777 0 0 1];% CTD #9 & cham#16 - only 1 profile
    elseif q.script.num==17
        head.coef.MHC=[2.5498    0.6709 0 0 1];% CTD #9 & cham#17 - only 2 profiles
    elseif q.script.num==18
        head.coef.MHC=[2.5499    0.6685 0 0 1];% CTD #9 & cham#18
    elseif q.script.num==19
        head.coef.MHC=[2.5473    0.6673 0 0 1];
    elseif q.script.num==20
        head.coef.MHC=[2.5473    0.6663 0 0 1];
    elseif q.script.num==21
        head.coef.MHC=[2.5429    0.6666 0 0 1];
    elseif q.script.num==22
        head.coef.MHC=[2.5331    0.6683 0 0 1];
    elseif q.script.num==23
        head.coef.MHC=[2.5302    0.6681 0 0 1];
    elseif q.script.num==24
        head.coef.MHC=[2.5280    0.6680 0 0 1];
    elseif q.script.num==25
        head.coef.MHC=[2.5250    0.6680 0 0 1];
    elseif q.script.num==26
        head.coef.MHC=[2.5220    0.6680 0 0 1];
    elseif q.script.num==27
        head.coef.MHC=[2.5200    0.6680 0 0 1];
    elseif q.script.num==28
        head.coef.MHC=[2.5180    0.6680 0 0 1];
    elseif q.script.num==29
        head.coef.MHC=[2.5160    0.6680 0 0 1];
    elseif q.script.num==30
        head.coef.MHC=[2.5140    0.6680 0 0 1];
    elseif q.script.num==31
        head.coef.MHC=[2.5120    0.6680 0 0 1];
    elseif q.script.num==32
        head.coef.MHC=[2.510    0.6680 0 0 1];
    elseif q.script.num==33
        head.coef.MHC=[2.508    0.6680 0 0 1];
    elseif q.script.num==34
        head.coef.MHC=[2.5059 0.6671 0 0 1];
    elseif q.script.num==35
        head.coef.MHC=[2.5117    0.6725 0 0 1];
    elseif q.script.num==36
        head.coef.MHC=[2.5105    0.6700 0 0 1];
    elseif q.script.num==37
        head.coef.MHC=[2.5191    0.6675 0 0 1];
    elseif q.script.num==38
        head.coef.MHC=[2.5231    0.6647 0 0 1];
    elseif q.script.num==39
        head.coef.MHC=[2.5203    0.6648 0 0 1];
    elseif q.script.num==40
        head.coef.MHC=[2.5136    0.6662 0 0 1];
    elseif q.script.num==41
        head.coef.MHC=[2.5077    0.6666 0 0 1];
    elseif q.script.num==42
        head.coef.MHC=[2.5090    0.6654 0 0 1];
    elseif q.script.num==43
        head.coef.MHC=[2.5060    0.6654 0 0 1];
    elseif q.script.num>=44 && q.script.num<=46
        head.coef.MHC=[2.5030    0.6653 0 0 1];
    elseif q.script.num>=47 && q.script.num<=49
        head.coef.MHC=[2.4942    0.6649 0 0 1];
    elseif q.script.num>=50 && q.script.num<=52
        head.coef.MHC=[2.4975    0.6624 0 0 1];
    elseif q.script.num>=53 && q.script.num<=55
        head.coef.MHC=[2.4885    0.6634 0 0 1];
    elseif q.script.num>=56 && q.script.num<=58
        head.coef.MHC=[2.4850    0.6631 0 0 1];
    elseif q.script.num>=59 && q.script.num<=61
        head.coef.MHC=[2.4883    0.6607 0 0 1];
    elseif q.script.num>=62 && q.script.num<=63
        head.coef.MHC=[2.4836    0.6616 0 0 1];
    elseif q.script.num>=64 && q.script.num<=65
        head.coef.MHC=[2.4858    0.6606 0 0 1];
    elseif q.script.num==66
        head.coef.MHC=[2.4826    0.6608 0 0 1];
    elseif q.script.num==358
        head.coef.MHC=[2.8118    0.4900    0.0224 0 1];% CTD #12 & cham#358 - 1 profile
    elseif q.script.num==359
        head.coef.MHC=[2.7451    0.5203    0.0186 0 1];% CTD #12 & cham#359 - 2 profiles
    elseif q.script.num>=360 && q.script.num<=363
        head.coef.MHC=[2.8415    0.4655    0.0256 0 1];
    elseif q.script.num>=364 && q.script.num<=366
        head.coef.MHC=[2.8593    0.4545    0.0267 0 1];
    elseif q.script.num>=367 && q.script.num<=371
        head.coef.MHC=[2.8199    0.4752    0.0238 0 1];
    elseif q.script.num>=372 && q.script.num<=376
        head.coef.MHC=[2.8382    0.4632    0.0252 0 1];
    elseif q.script.num>=377 && q.script.num<=381
        head.coef.MHC=[2.8085    0.4813    0.0224 0 1];
    elseif q.script.num>=382 && q.script.num<=386
        head.coef.MHC=[2.8743    0.4401    0.0280 0 1];
    elseif q.script.num>=387 && q.script.num<=406
        head.coef.MHC=[2.9116    0.4178    0.0309 0 1];
    elseif q.script.num>=407 && q.script.num<=416
        head.coef.MHC=[2.8618    0.4429    0.0275 0 1];
    elseif q.script.num>=417 && q.script.num<=426
        head.coef.MHC=[2.8783    0.4320    0.0289 0 1];
    elseif q.script.num>=427 && q.script.num<=446
        head.coef.MHC=[2.8323    0.4576    0.0251 0 1];
    elseif q.script.num>=447 && q.script.num<=458
        head.coef.MHC=[2.6948    0.5284    0.0164 0 1];
    elseif q.script.num>=459 && q.script.num<=478
        head.coef.MHC=[2.8641    0.4301    0.0296 0 1];
    elseif q.script.num>=479 && q.script.num<=498
        head.coef.MHC=[2.8462    0.4401    0.0279 0 1];
    elseif q.script.num>=499 && q.script.num<=510
        head.coef.MHC=[2.8941    0.4123    0.0314 0 1];
    elseif q.script.num>=510 && q.script.num<=518
        head.coef.MHC=[2.8658    0.4263    0.0296 0 1];
    elseif q.script.num>=519 && q.script.num<=528
        head.coef.MHC=[2.7783    0.4789    0.0219 0 1];
    elseif q.script.num>=529 && q.script.num<=531
        head.coef.MHC=[2.9033    0.4030    0.0325 0 1];
    elseif q.script.num>=532 && q.script.num<=538
        head.coef.MHC=[2.9279    0.3869    0.0348 0 1];
    elseif q.script.num>=539 && q.script.num<=541
        % CTD #13 & cham#539,540,541
        head.coef.MHC=[2.9196    0.3935    0.0335 0 1];
    end
elseif ~isempty(strfind(head.sensor_id(7,:),'08-03')) && ~isempty(strfind(head.instrument,'CHAM04-04'))% CHAM 04-04
    %  in water casts 1851 - 1866 
    if q.script.num==1851% CTD #23 & cham#1851
        head.coef.MHC=[-0.0752    2.7418 0 0 1];
    elseif q.script.num==1852% CTD #23 & cham#1852
        head.coef.MHC=[-0.0528    2.7171 0 0 1];
    elseif q.script.num==1853 % CTD #23 & cham#1853
        head.coef.MHC=[-0.0549    2.7015 0 0 1];
    elseif q.script.num==1854 % CTD #23 & cham#1854
        head.coef.MHC=[-0.1204    2.7384 0 0 1];
    elseif q.script.num==1855 % CTD #23 & cham#1855
        head.coef.MHC=[-0.0651    2.7044 0 0 1];
    elseif q.script.num==1856 % CTD #23 & cham#1856
        head.coef.MHC=[-0.0615    2.6987 0 0 1];
    elseif q.script.num==1857
        head.coef.MHC=[-0.0495    2.6904 0 0 1];
    elseif q.script.num==1858
        head.coef.MHC=[-0.0479    2.6873 0 0 1];
    elseif q.script.num==1859
        head.coef.MHC=[-0.0567    2.6902 0 0 1];
    elseif q.script.num>=1860 && q.script.num<=1862
        head.coef.MHC=[-0.0529    2.6855 0 0 1];
    elseif q.script.num>=1863 && q.script.num<=1866
        head.coef.MHC=[-0.0530    2.6832 0 0 1];
    end
elseif ~isempty(strfind(head.sensor_id(7,:),'114')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==67 
        head.coef.MHC=[3.0056    0.7700 0 0 1];
    elseif q.script.num==68 
        head.coef.MHC=[3.0090    0.7668 0 0 1];
    elseif q.script.num==69 
        head.coef.MHC=[3.0035    0.7672 0 0 1];
    elseif q.script.num==70 
        head.coef.MHC=[2.9981    0.7693 0 0 1];
    elseif q.script.num==71 
        head.coef.MHC=[2.9896    0.7714 0 0 1];
    elseif q.script.num==72 
        head.coef.MHC=[2.9896    0.7714 0 0 1];
    elseif q.script.num==73 
        head.coef.MHC=[2.9871    0.7699 0 0 1];
    elseif q.script.num==74 
        head.coef.MHC=[2.9871    0.7699 0 0 1];
    elseif q.script.num==75 
        head.coef.MHC=[2.9865    0.7684 0 0 1];
    elseif q.script.num==76 
        head.coef.MHC=[2.9865    0.7684 0 0 1];
    elseif q.script.num==546 %!!! board has been changed before drop 546
        head.coef.MHC=[-2.5579    5.0439   -0.6285 0 1];% CTD #14 & cham#546 - 1 profile
    elseif q.script.num>=547  && q.script.num<=548
        head.coef.MHC=[3.4030   -0.9509    0.8503 0 1];% CTD #14 & cham#547
    elseif q.script.num>=549  && q.script.num<=550 
        head.coef.MHC=[3.5402   -1.0983    0.8898 0 1];
    elseif q.script.num==551 
        head.coef.MHC=[-4.2265    6.6519   -1.0296 0 1];
    elseif q.script.num==552 
        head.coef.MHC=[0.0527    2.4382 0 0 1];
    elseif q.script.num==553 
        head.coef.MHC=[-0.0002    2.4654 0 0 1];
    elseif q.script.num==554 
        head.coef.MHC=[-0.0679    2.4997 0 0 1];
    elseif q.script.num>=555  && q.script.num<=558
        head.coef.MHC=[-0.0608    2.4963 0 0 1];
    elseif q.script.num>=559  && q.script.num<=563
        head.coef.MHC=[-0.0651    2.4972 0 0 1];
    elseif q.script.num>=564  && q.script.num<=575
        head.coef.MHC=[-0.0591    2.4928 0 0 1];
    elseif q.script.num>=576  && q.script.num<=585
        head.coef.MHC=[-0.0544    2.4860 0 0 1];
    elseif q.script.num>=586  && q.script.num<=615
        head.coef.MHC=[-0.0653    2.4894 0 0 1];
    elseif q.script.num>=596  && q.script.num<=615
        head.coef.MHC=[-0.0381    2.4720 0 0 1];
    elseif q.script.num>=616  && q.script.num<=635
        head.coef.MHC=[-0.0634    2.4842 0 0 1];
    elseif q.script.num>=636  && q.script.num<=655
        head.coef.MHC=[-0.0340    2.4679 0 0 1];
    elseif q.script.num>=656  && q.script.num<=675
        head.coef.MHC=[-0.0542    2.4791 0 0 1];
    elseif q.script.num>=676  && q.script.num<=695
        head.coef.MHC=[-0.0565    2.4792 0 0 1];
    elseif q.script.num>=696  && q.script.num<=715
        head.coef.MHC=[-0.0628    2.4826 0 0 1];
    elseif q.script.num>=716  && q.script.num<=739
        head.coef.MHC=[-0.0664    2.4831 0 0 1];
    elseif q.script.num>=740 && q.script.num<=796
         head.coef.MHC=[-0.0652    2.4821 0 0 1];% CTD #15 & cham#740
    elseif q.script.num==797
         head.coef.MHC=[-0.0447    2.4856 0 0 1];% cham was out of water
    elseif q.script.num==798
         head.coef.MHC=[-0.0447    2.4856 0 0 1];
    elseif q.script.num>=799 && q.script.num<=801
         head.coef.MHC=[-0.0318    2.4763 0 0 1];
    elseif q.script.num>=802 && q.script.num<=806
         head.coef.MHC=[-0.0380    2.4769 0 0 1];
    elseif q.script.num>=807 && q.script.num<=811
         head.coef.MHC=[-0.0211    2.4646 0 0 1];
    elseif q.script.num>=812 && q.script.num<=840
         head.coef.MHC=[-0.0310    2.4714 0 0 1];
    elseif q.script.num>=841 && q.script.num<=860
         head.coef.MHC=[-0.0239    2.4658 0 0 1];
    elseif q.script.num>=861 && q.script.num<=880
         head.coef.MHC=[-0.0499    2.4809 0 0 1];
    elseif q.script.num>=881 && q.script.num<=911% CTD #16 & cham#910
         head.coef.MHC=[-0.0211    2.4633 0 0 1];
    elseif q.script.num==912% cham was out of water
         head.coef.MHC=[-0.0814    2.5023 0 0 1];
    elseif q.script.num==913
         head.coef.MHC=[-0.0772    2.4976 0 0 1];
    elseif q.script.num==914
         head.coef.MHC=[-0.0655    2.4899 0 0 1];
    elseif q.script.num>=915 && q.script.num<=917
         head.coef.MHC=[-0.0601    2.4850 0 0 1];
    elseif q.script.num>=918 && q.script.num<=922
         head.coef.MHC=[-0.0730    2.4887 0 0 1];
    elseif q.script.num>=923 && q.script.num<=980
         head.coef.MHC=[-0.0515    2.4758 0 0 1];
    elseif q.script.num>=981 && q.script.num<=1098
         head.coef.MHC=[-0.0515    2.4758 0 0 1];
    elseif q.script.num==1099
         head.coef.MHC=[-0.0591    2.4774 0 0 1];
    elseif q.script.num>=1021 && q.script.num<=1040
         head.coef.MHC=[-0.0602    2.4790 0 0 1];
    elseif q.script.num>=1041 && q.script.num<=1079
         head.coef.MHC=[-0.0756    2.4867 0 0 1];
    elseif q.script.num>=1080 && q.script.num<=1120
         head.coef.MHC=[-0.0700    2.4839 0 0 1];
    elseif q.script.num>=1121 && q.script.num<=1141
         head.coef.MHC=[-0.0633    2.4797 0 0 1];
    elseif q.script.num==1142 % cham was out of water
         head.coef.MHC=[-0.0795    2.4889 0 0 1];
    elseif q.script.num>=1143 && q.script.num<=1160
         head.coef.MHC=[-0.0737    2.4858 0 0 1];
    elseif q.script.num>=1161 && q.script.num<=1170
         head.coef.MHC=[-0.0793    2.4885 0 0 1];
    elseif q.script.num>=1171 && q.script.num<=1180
         head.coef.MHC=[-0.06    2.4785 0 0 1];
    elseif q.script.num>=1181 && q.script.num<=1190
         head.coef.MHC=[-0.055    2.4725 0 0 1];
    elseif q.script.num>=1191 && q.script.num<=1200
         head.coef.MHC=[-0.05    2.4705 0 0 1];
    elseif q.script.num>=1201 && q.script.num<=1230
         head.coef.MHC=[-0.0470    2.4685 0 0 1];
    elseif q.script.num>=1231 && q.script.num<=1250
         head.coef.MHC=[-0.0470    2.4685 0 0 1];
    elseif q.script.num>=1251 && q.script.num<=1270
         head.coef.MHC=[-0.0608    2.4752 0 0 1];
    elseif q.script.num>=1271 && q.script.num<=1278
         head.coef.MHC=[-0.0613    2.4752 0 0 1];% CTD #18 & cham#1278
    end
elseif ~isempty(strfind(head.sensor_id(7,:),'113 4')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==77
        head.coef.MHC=[2.9715    0.7807 0 0 1];
    elseif q.script.num==78
        head.coef.MHC=[2.9772    0.7759 0 0 1];
    elseif q.script.num>=79 && q.script.num<=81
        head.coef.MHC=[2.9724    0.7741 0 0 1];
    elseif q.script.num>=82 && q.script.num<=87
        head.coef.MHC=[2.9710    0.7716 0 0 1];
    elseif q.script.num>=88 && q.script.num<=94
        head.coef.MHC=[2.9643    0.7716 0 0 1];
    elseif q.script.num>=95 && q.script.num<=104
        head.coef.MHC=[2.9649    0.7688 0 0 1];
    elseif q.script.num>=105 && q.script.num<=114
        head.coef.MHC=[2.9580    0.7695 0 0 1];
    elseif q.script.num>=115 && q.script.num<=134
        head.coef.MHC=[2.9583    0.7678 0 0 1];
    elseif q.script.num>=135 && q.script.num<=154
        head.coef.MHC=[2.9561    0.7684 0 0 1];
    elseif q.script.num>=155 && q.script.num<=164
        head.coef.MHC=[2.9572    0.7662 0 0 1];
    elseif q.script.num>=165 && q.script.num<=174
        head.coef.MHC=[2.9515    0.7682 0 0 1];
    elseif q.script.num>=175 && q.script.num<=196
        head.coef.MHC=[2.9523    0.7670 0 0 1];
    elseif q.script.num==197
        head.coef.MHC=[2.9581    0.7625 0 0 1];
    elseif q.script.num==198
        head.coef.MHC=[2.9572    0.7633 0 0 1];
    elseif q.script.num==199
        % CTD #11 & cham#199 - long time in the water
        head.coef.MHC=[2.9535 0.7663 0 0 1];% first order is a bit worse than second order
    elseif q.script.num==200
        % CTD #11 & cham#200 - instrument was out of water; 3 hours after the CTD
        head.coef.MHC=[2.9598    0.7704 0 0 1];% first order is worse than second order
    elseif q.script.num==201
        head.coef.MHC=[2.9598    0.7704 0 0 1];
    elseif q.script.num>=202 && q.script.num<=204
        head.coef.MHC=[2.9605    0.7692 0 0 1];
    elseif q.script.num>=205 && q.script.num<=207
        head.coef.MHC=[2.9594    0.7675 0 0 1];
    elseif q.script.num>=208 && q.script.num<=210
        head.coef.MHC=[2.9550    0.7676 0 0 1];
    elseif q.script.num>=211 && q.script.num<=215
        head.coef.MHC=[2.9574    0.7660 0 0 1];
    elseif q.script.num>=216 && q.script.num<=225
        head.coef.MHC=[2.9575    0.7649 0 0 1];
    elseif q.script.num>=226 && q.script.num<=235
        head.coef.MHC=[2.9550    0.7658 0 0 1];
    elseif q.script.num>=236 && q.script.num<=240
        head.coef.MHC=[2.9509    0.7671 0 0 1];
    elseif q.script.num>=241 && q.script.num<=270
        head.coef.MHC=[2.9471    0.7669 0 0 1];
    elseif q.script.num>=271 && q.script.num<=285
        head.coef.MHC=[2.9440    0.7685 0 0 1];
    elseif q.script.num>=286 && q.script.num<=305
        head.coef.MHC=[2.9476    0.7685 0 0 1];
    elseif q.script.num>=306 && q.script.num<=319
        head.coef.MHC=[2.9385    0.7710 0 0 1];
    elseif q.script.num>=320 && q.script.num<=339
        head.coef.MHC=[2.9394    0.7694 0 0 1];
    elseif q.script.num>=340 && q.script.num<=356
        head.coef.MHC=[2.9472    0.7672 0 0 1];
    elseif q.script.num==357
        % CTD #12 & cham#357 - long time in the water
        head.coef.MHC=[2.9475    0.7664 0 0 1];% first order is a bit worse than second order
    end
elseif ~isempty(strfind(head.sensor_id(7,:),'04_059')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==542 % 1 profile
        head.coef.MHC=[-2.2143    4.6726   -0.5267 0 1];% CTD #13 & cham#542
    elseif q.script.num>=543 && q.script.num<=545
        head.coef.MHC=[-2.5579    5.0439   -0.6285 0 1];
    end
elseif ~isempty(strfind(head.sensor_id(7,:),'04_062')) && ~isempty(strfind(head.instrument,'CHAM04-03'))% CHAM 04-03
    if q.script.num==1279 % cham was out of water
         head.coef.MHC=[-0.0584    2.4772 0 0 1];% CTD #18 & cham#1279
    elseif q.script.num==1280 
         head.coef.MHC=[-0.0586    2.4755 0 0 1];
    elseif q.script.num==1281 
         head.coef.MHC=[-0.0647    2.4783 0 0 1];
    elseif q.script.num==1282 
         head.coef.MHC=[-0.0483    2.4697 0 0 1];
    elseif q.script.num==1283 
         head.coef.MHC=[-0.0754    2.4826 0 0 1];
    elseif q.script.num>=1284 && q.script.num<=1287
         head.coef.MHC=[-0.0713    2.4805 0 0 1];
    elseif q.script.num>=1288 && q.script.num<=1293
         head.coef.MHC=[-0.0550    2.4716 0 0 1];
    elseif q.script.num>=1294 && q.script.num<=1310
         head.coef.MHC=[-0.0626    2.4749 0 0 1];
    elseif q.script.num>=1311 && q.script.num<=1330
         head.coef.MHC=[-0.0709    2.4782 0 0 1];
    elseif q.script.num>=1331 && q.script.num<=1340
         head.coef.MHC=[-0.07    2.4765 0 0 1];
    elseif q.script.num>=1341 && q.script.num<=1416 
         head.coef.MHC=[-0.0695    2.4754 0 0 1];% CTD #19 & cham#1416
    elseif q.script.num==1417 % cham was out of water
         head.coef.MHC=[-0.0440    2.4620 0 0 1];% CTD #19 & cham#1417
    elseif q.script.num==1418
         head.coef.MHC=[-0.0790    2.4818 0 0 1];
    elseif q.script.num==1419
         head.coef.MHC=[-0.0824    2.4830 0 0 1];
    elseif q.script.num>=1420 && q.script.num<=1421
         head.coef.MHC=[-0.0550    2.4678 0 0 1];
    elseif q.script.num==1422 && q.script.num<=1424
         head.coef.MHC=[-0.0690    2.4751 0 0 1];
    elseif q.script.num>=1423 && q.script.num<=1424
         head.coef.MHC=[-0.0532    2.4709 0 0 1];
    elseif q.script.num>=1425 && q.script.num<=1430
         head.coef.MHC=[-0.0582    2.4731 0 0 1];
    elseif q.script.num>=1431 && q.script.num<=1434
         head.coef.MHC=[-0.0471    2.4675 0 0 1];
    elseif q.script.num>=1435 && q.script.num<=1438
         head.coef.MHC=[-0.0431    2.4673 0 0 1];
    elseif q.script.num>=1439 && q.script.num<=1440
         head.coef.MHC=[-0.0488    2.4674 0 0 1];
    elseif q.script.num>=1441 && q.script.num<=1443
         head.coef.MHC=[-0.0579    2.4729 0 0 1];
    elseif q.script.num>=1444 && q.script.num<=1447
         head.coef.MHC=[-0.0559    2.4750 0 0 1];
    elseif q.script.num>=1448 && q.script.num<=1450
         head.coef.MHC=[-0.0557    2.4716 0 0 1];
    elseif q.script.num==1451
         head.coef.MHC=[-0.0567    2.4749 0 0 1];
    elseif q.script.num>=1452 && q.script.num<=1463
         head.coef.MHC=[-0.0420    2.4639 0 0 1];
    elseif q.script.num>=1464 && q.script.num<=1465
         head.coef.MHC=[-0.0456    2.4657 0 0 1];
    elseif q.script.num>=1466 && q.script.num<=1467
         head.coef.MHC=[-0.0673    2.4811 0 0 1];
    elseif q.script.num>=1468 && q.script.num<=1475
         head.coef.MHC=[-0.0624    2.4736 0 0 1];
    elseif q.script.num>=1476 && q.script.num<=1501
         head.coef.MHC=[-0.0538    2.4697 0 0 1];
    elseif q.script.num>=1502 && q.script.num<=1506
         head.coef.MHC=[-0.0572    2.4757 0 0 1];
    elseif q.script.num==1507 % cham was out of water
         head.coef.MHC=[-0.0512    2.4721 0 0 1];% CTD #20 & cham#1507
    elseif q.script.num==1508
         head.coef.MHC=[-0.0577    2.4749 0 0 1];
    elseif q.script.num==1509
         head.coef.MHC=[-0.0471    2.4694 0 0 1];
    elseif q.script.num==1510
         head.coef.MHC=[-0.0433    2.4683 0 0 1];
    elseif q.script.num>=1511 && q.script.num<=1513
         head.coef.MHC=[-0.0539    2.4734 0 0 1];
    elseif q.script.num>=1514 && q.script.num<=1518
         head.coef.MHC=[-0.0685    2.4805 0 0 1];
    elseif q.script.num>=1519 && q.script.num<=1530
         head.coef.MHC=[-0.0653    2.4800 0 0 1];
    elseif q.script.num>=1531 && q.script.num<=1540
         head.coef.MHC=[-0.0845    2.4894 0 0 1];
    elseif q.script.num>=1541 && q.script.num<=1560
         head.coef.MHC=[-0.0674    2.4814 0 0 1];
    elseif q.script.num>=1561 && q.script.num<=1580
         head.coef.MHC=[-0.0693    2.4816 0 0 1];
    elseif q.script.num>=1581 && q.script.num<=1600
         head.coef.MHC=[-0.0538    2.4722 0 0 1];
    elseif q.script.num>=1601 && q.script.num<=1608
         head.coef.MHC=[-0.0561    2.4723 0 0 1]; 
    elseif q.script.num>=1609 && q.script.num<=1611
         head.coef.MHC=[-0.0562    2.4671 0 0 1]; 
    elseif q.script.num==1612
         head.coef.MHC=[-0.0562    2.4671 0 0 1]; 
    elseif q.script.num>=1613 && q.script.num<=1624
         head.coef.MHC=[-0.0545    2.4730 0 0 1];% CTD #21 & cham#1624 
    elseif q.script.num>=1631 && q.script.num<=1632 % cham was out of water
         head.coef.MHC=[-0.0116    2.4516 0 0 1];
    elseif q.script.num>=1633 && q.script.num<=1640
         head.coef.MHC=[-0.0369    2.4666 0 0 1];
    elseif q.script.num>=1641 && q.script.num<=1654
         head.coef.MHC=[-0.0088    2.4481 0 0 1];
    elseif q.script.num>=1655 && q.script.num<=1662 % cham was out of water
         head.coef.MHC=[-0.0496    2.4697 0 0 1];
    elseif q.script.num>=1663 && q.script.num<=1850
         head.coef.MHC=[-0.0412    2.4640 0 0 1];
    elseif q.script.num==1867
         head.coef.MHC=[-0.0598    2.4804 0 0 1];
    elseif q.script.num==1868
         head.coef.MHC=[-0.0506    2.4740 0 0 1];
    elseif q.script.num>=1869 && q.script.num<=1874
         head.coef.MHC=[-0.0487    2.4711 0 0 1];
    elseif q.script.num>=1875 && q.script.num<=1880
         head.coef.MHC=[-0.0480    2.4686 0 0 1];
    elseif q.script.num>=1875 && q.script.num<=1880
         head.coef.MHC=[-0.0480    2.4686 0 0 1];
    elseif q.script.num>=1881 && q.script.num<=1910
         head.coef.MHC=[-0.0429    2.4650 0 0 1];
    elseif q.script.num>=1911 && q.script.num<=1950
         head.coef.MHC=[-0.0504    2.4695 0 0 1];
    elseif q.script.num>=1951 && q.script.num<=2260
         head.coef.MHC=[-0.0504    2.4695 0 0 1];
    elseif q.script.num>=2261 && q.script.num<=2456
         head.coef.MHC=[-0.0539    2.4688 0 0 1];
    elseif q.script.num==2457
         head.coef.MHC=[-0.0539    2.4688 0 0 1];
    elseif q.script.num>=2458 && q.script.num<=2464
         head.coef.MHC=[-0.0539    2.4688 0 0 1];
    elseif q.script.num>=2465 && q.script.num<=2668
         head.coef.MHC=[-0.0417    2.4609 0 0 1];
    end
end

