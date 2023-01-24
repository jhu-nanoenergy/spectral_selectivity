clear all; close all
load('circ_all.mat')
load('PbS-CQD-Model.mat')
ASTMG173=readtable('ASTMG173.txt','ReadVariableNames',true);
colNames=ASTMG173(1,:);
ASTMG173(1,:) = [];
data=table2array(ASTMG173); %Converts data into matrix


phc_type = sqr_hol;

vis_indices = [1, 76];
k_norm = fdtd_k./max(fdtd_k);
w = 0.5;
for i = 1:length(phc_type)
    A_norm = phc_type(i).A./max(phc_type(i).A);
    vis_integral(i) = trapz(A_norm(vis_indices(1):vis_indices(2)));
    nir_integral(i) = trapz(A_norm(vis_indices(2):221));
    spec_integral(i)= (w.*nir_integral(i)) + (1-w).* (-1.*vis_integral(i));

end
[M,Ispec] = max(spec_integral);

a = phc_type(Ispec).a
r_a = phc_type(Ispec).r_a
t_a = phc_type(Ispec).t_a


%%%Calculating Equivalent absorption
d = phc_type(Ispec).thk_eff.*1e-9;
n1 = 1; %Air
n2 = fdtd_n + 1i.*fdtd_k;
r12 = (n2 - n1)./(n1 + n2);
r21 = (n1 - n2)./(n1 + n2);
t12 = (2*n1)./(n1 + n2);
t21 = (2*n2)./(n1 + n2);
delta = 2 .* pi .* d.* (1./(wvl_nm.*1e-9)) .* sqrt((n2).^2);
r123 = (r12 + r21.*exp(2*1i*delta))./(1 + r12.*r21.*exp(2*1i*delta));
t123 = (t12 .* t21 .* exp(1i*delta))./(1 + r12.*r21.*exp(2*1i*delta));
TMM_transmitted = (abs(t123)).^2;
TMM_reflected = (abs(r123)).^2;
TMM_absorbed = 1 - TMM_reflected - TMM_transmitted;


%Plotting
hold on
grid on
plot(wvl_nm, phc_type(Ispec).A,'r')
plot(wvl_nm, phc_type(Ispec).T,'b')
plot(wvl_nm, TMM_absorbed, 'r-.');
plot(wvl_nm, TMM_transmitted, 'b-.');
patch([400 750 750 400], [0 0 1 1], [0.8 0.8 0.8], 'EdgeAlpha', 0, 'FaceAlpha', 0.3)
legend('Abs', 'Trn', 'Abs-Eqv', 'Trn-Eqv', 'Location', 'best')
xlabel('wvl [nm]')
ylabel('Power')
axis([400 1500 0 1])

absorption_gain = mean(TMM_absorbed(1:76)-phc_type(Ispec).A(1:76).')
transmission_gain = mean(phc_type(Ispec).T(1:76).'-TMM_transmitted(1:76))
