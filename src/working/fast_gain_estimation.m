clear all; close all;
sys_case=39; method=10;
compute_gain;
[num,den] = tfdata(G{4}(2,3));
sys_order=size(den{1},2);

b=num{1}; a=den{1};
[r,p,k] = residue(b,a);
[y_impulse,t_impulse]=impulse(tf(b,a));
t_impulse(3)-t_impulse(3);

threshold_triv=1e-2; threshold_T=1e-1;
idx_nontriv=find((abs(real(r)./real(p))>threshold_triv).*(abs(real(p))>1e-5));
idx_triv=find((abs(real(r)./real(p))<threshold_triv).*(abs(real(p))>1e-5));
t_end=max(abs(log(abs(real(p(idx_nontriv))./real(r(idx_nontriv)))*threshold_T)./real(p(idx_nontriv))));
time_step=min(2*pi/10/max(abs(imag(p(idx_nontriv)))),0.1*min(real(p(idx_nontriv))*log(0.1)));
t_impulse_anly=0:time_step:t_end;
y_impulse_anly=real(r(idx_nontriv).'*exp(p(idx_nontriv)*t_impulse_anly));

%subplot(2,1,1); 
hold all;
plot(t_impulse,y_impulse)
plot(t_impulse_anly,y_impulse_anly,'r--')
%plot(t_impulse,exp(real(p(1))*t_impulse).*(2*real(r(1))*cos(imag(p(1)*t_impulse))-2*imag(r(1))*sin(imag(p(1)*t_impulse))),'g:')
%subplot(2,1,2)
%plot(t_impulse,diag(r)*exp(real(p)*t_impulse'),'r--')
disp(['Gain from simulation: ' num2str(trapz(t_impulse,abs(y_impulse)))])
disp(['Gain from estimation: ' num2str(trapz(t_impulse_anly,abs(y_impulse_anly))+sum(abs(real(r(idx_triv))./real(p(idx_triv)))))])
disp(['Trivial/Nontrivial: ' num2str(size(idx_triv,1)) '/' num2str(size(idx_nontriv,1))])
disp(['Contribution from nontrivial: ' num2str(sum(abs(real(r(idx_triv))./real(p(idx_triv)))))])

