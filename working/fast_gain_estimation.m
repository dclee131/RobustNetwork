clear all; close all;
sys_case=9; compute_gain_method=10; loading_level=100;
compute_gain;
[num,den] = tfdata(prescale(G{4}(1,1)));
sys_order=size(den{1},2);

b=num{1}; a=den{1};
tic;
[y_impulse,t_impulse]=impulse(tf(b,a));
gain_sim=trapz(t_impulse,abs(y_impulse));
cp_sim=toc;

tic;
[r,p,k] = residue(b,a);
threshold_reduction=1e-2;
idx_nontriv=find((abs(real(r)./real(p))>threshold_reduction).*(abs(real(p))>1e-5));
idx_triv=find((abs(real(r)./real(p))<threshold_reduction).*(abs(real(p))>1e-5));
[b_reduced,a_reduced] = residue(r(idx_nontriv),p(idx_nontriv),k);
[y_reduced,t_reduced]=impulse(tf(real(b_reduced),real(a_reduced)));
gain_reduced=trapz(t_reduced,abs(y_reduced))+sum(abs(real(r(idx_triv))./real(p(idx_triv))));
cp_reduced=toc;

figure; hold all;
plot(t_impulse,y_impulse)
plot(t_reduced,y_reduced,'r--')

disp(['Gain from simulation: ' num2str(gain_sim)])
disp(['Gain from estimation: ' num2str(gain_reduced)])
disp(['Contribution from nontrivial: ' num2str(size(idx_nontriv,1)) '/' num2str(sys_order)])
disp(['Computation Time: ' num2str(cp_sim) '/' num2str(cp_reduced)])

