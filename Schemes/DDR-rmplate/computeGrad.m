clc
clear all 
close all

syms x y

%%%% polynomial test
fct = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)

Df = [diff(fct,x); diff(fct,y)];

disp('gradient')
ccode(Df(1))
ccode(Df(2))

