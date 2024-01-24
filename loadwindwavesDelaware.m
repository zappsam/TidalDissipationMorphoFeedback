clear;close all;clc
a=dlmread('ST63154_v03.onlns');

a=[a(:,5:6) a(:,10) a(:,12) a(:,16)];
windDELAWARE=a;
