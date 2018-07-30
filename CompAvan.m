clear;
clc;
%%Dados requeridos
%{
Kp=lim (s->0) Gc(s)*G(s)
ess(tipo 0) = 1/(1+Kp)

Kv=lim (s->0) s*Gc(s)*G(s) 
%}
Kc=19;
qsi=0.707;
wn=36;
%%Calculo Tz e Tp
s1=-qsi*wn+i*wn*sqrt(1-qsi*qsi);

s=tf('s');
FTMA=tf(500/((s+10)*(s+50)));

Gs=evalfr(FTMA,s1);

Mg=abs(Gs);
Og=angle(Gs);

Ms=abs(s1);
Os= angle(s1);

Tz=(sin(Os)-Kc*Mg*sin(Og-Os))/(Kc*Mg*Ms*sin(Og))
Tp=-((Kc*Mg*sin(Os)+sin(Og+Os))/(Ms*sin(Og)))
%{
Avanço de Fase Tz>Tp

%}
Gc=tf(Kc*((Tz*s+1)/(Tp*s+1)));

%% Resultados
H=tf(1);

A = feedback(FTMA, H);
B = feedback(Gc*FTMA, H);

figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2, 4, [1 2]);
step(A, B);
grid;
legend('FTMA','FTMA*Gc');

subplot(2, 4, 3);
rlocus(FTMA);
title('Root Locus Original');

subplot(2, 4, 4);
rlocus(FTMA*Gc);
title('Root Locus Compensado');

subplot(2, 4, [5 6]);
bode(FTMA, FTMA*Gc);
grid;
legend('FTMA','FTMA*Gc');

subplot(2, 4, 7);
nyquist(FTMA);
title('Nyquist Original');

subplot(2, 4, 8);
nyquist(FTMA*Gc);
title('Nyquist Compensado');
