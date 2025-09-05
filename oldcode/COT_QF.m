function [E, G] = COT_QF(coeffs, t, W, uField, vField, wField, dt, p, B)

Ufactor = 0;

[uField, vField, wField] = ChangeU(uField,vField,wField, Ufactor);


U = calculateRMS(uField,vField,wField);


[E, G] = COT_function(coeffs,t,W,StartLoc,uField,vField,wField,dt, p,U, B);