clear all
clc
mr=0.1;
kr=2*0.5;
wr=kr/mr;
syms 
AF=@(lamda) abs((1-wr*lamda^2)/((1+kr-lamda^2)*(1-wr*lamda^2)-kr));

fplot(AF,[0 2])
ylim([0 5])