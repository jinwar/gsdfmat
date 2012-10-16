clear;

[sachdr data]=load_sac('IRM.facor');

x=sachdr.b:sachdr.delta:sachdr.e;

para=load('IRM.facor.fit');

ao=para(1);
so=para(2);
wo=para(3);
tp=para(4);
tg=para(5);

y=ao*exp(-0.5*(so*(x-tg)).^2).*cos(wo*(x-tp));


figure(1)
clf
hold on
plot(x,data)
plot(x,y,'r--');
set(gca,'FontSize',12);
title('Cross Correlation Waveform');
xlabel('Time Lag /s');
ylabel('Amplitude');
print('-depsc','xcorfit');
system('ps2pdf xcorfit.eps');
