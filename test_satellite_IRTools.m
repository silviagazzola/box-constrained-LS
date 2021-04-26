clear, clc
close all
n = 256;

% deblurring
% optbl.trueImage = 'satellite';
% optbl.PSF = 'speckle';
optbl.PSF = 'shake';
optbl.level = 'severe';
[A, b, x, PbInfo] = PRblur(n, optbl);
nl = 1e-2;
bn = PRnoise(b, nl);

% % inpainting
% p = 0.6; % p = 0.4; % p = 0.3; % p=0.5;
% rng(0)
% A1 = A;
% % Mask
% masktemp = rand(n);
% mask = zeros(n); mask(masktemp<p) = 1;
% window = zeros(21); hsw = 2; window(11-hsw:11+hsw, 11-hsw:11+hsw) = 1/(hsw+1)^2;
% mask1 = conv2(mask,window); hmask = 1.7; mask1 = mask1(11:(size(mask1,1)-10),11:(size(mask1,2)-10)); mask1(mask1<=hmask) = 0; mask1(mask1>hmask) = 1;
% ind = find(mask1==1);
% M = numel(ind);
% %Mask matrix (sparse matrix in matlab)
% S = sparse(1:M, ind, ones(M, 1), M, n^2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % finally generate the blur-inpaint operator
% A = @(xx,tflag) OPblurinpaint(xx,A1,S,tflag);
% 
% A2 = @(xx,tflag) OPblurinpaintvar(xx,A1,S,ind,n,tflag);
% 
% % generate data
% b = A_times_vec(A, x);
% bdisp = zeros(n^2,1); bdisp(ind) = b;
% undersampling = nnz(bdisp)/n^2;
% nl = 1e-3; % nl = 0; % 
% bn = PRnoise(b, nl);
% b0 = bdisp;
b0 = bn;
b0(b0<0)=0;

m_in = 20; m_out = 40; % m_out = 10; %

disp('NN-FCGLS')
stopC_1.kind = 'discr';
stopC_1.noise = nl;
[Xiter_RS1, RelErr_RS1, RelRes_RS1, Alpha_RS1, tot_it1, cycle_it1, stopC_it1] = fcgls_rp_nn_ReSt (A, bn, m_in, m_out, x, b0, 1, stopC_1);

disp('BC-FCGLS')
stopC_1.kind = 'discr';
stopC_1.noise = nl;
b0(b0>1)=1;
% m_out = 80; 
[Xiter_RS2, RelErr_RS2, RelRes_RS2, Alpha_RS2, tot_it2, cycle_it2, stopC_it2] = fcgls_bc_ReSt (A, bn, 0, 1, m_in, m_out, x, b0, 1, stopC_1);
% 
optcg.NoStop = 'on';
optcg.x_true = x;
[Xcg, infocg] = IRcgls(A, bn, optcg);
% 
% % figure, semilogy(RelErr_RS1)
% % hold on
% % semilogy(RelErr_RS2)
% % semilogy(infofista.Enrm)
% % 
% % figure, imagesc(reshape(Xiter_RS1(:,end), n, n)), axis image, axis off
% % figure, imagesc(reshape(Xiter_RS2(:,end), n, n)), axis image, axis off
% % figure, imagesc(reshape(Xfista(:,end), n, n)), axis image, axis off
% 
optfista.NoStop = 'on';
optfista.xMax = 1;
optfista.x_true = x;
[Xfista, infofista] = IRfista(A, bn, optfista);
% 
figure, semilogy(infocg.Enrm, 'LineWidth', 2)
hold on
semilogy(RelErr_RS1, 'LineWidth', 2)
semilogy(RelErr_RS2, 'LineWidth', 2)
semilogy(infofista.Enrm, 'LineWidth', 2)
legend('CGLS', 'NN-FCGLS', 'box-FCGLS', 'FISTA')
set(gca,'FontSize',22)
% 
% figure, semilogy(infocg.Rnrm, 'LineWidth', 2)
% hold on
% semilogy(RelRes_RS1, 'LineWidth', 2)
% semilogy(RelRes_RS2, 'LineWidth', 2)
% semilogy(infofista.Rnrm, 'LineWidth', 2)
% legend('CGLS', 'NN-FCGLS', 'box-FCGLS', 'FISTA')
% set(gca,'FontSize',22)
% 
figure, imagesc(reshape(infocg.BestReg.X, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
[val1, ind1] = min(RelErr_RS1);
figure, imagesc(reshape(Xiter_RS1(:,ind1), n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
[val2, ind2] = min(RelErr_RS2);
figure, imagesc(reshape(Xiter_RS2(:,ind2), n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
figure, imagesc(reshape(infofista.BestReg.X, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
% 
figure, imagesc(reshape(x, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
figure, imagesc(reshape(bn, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
% figure, imagesc(reshape(PbInfo.PSF, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
% 
% % [val1, ind1] = min(RelErr_RS1);
% % figure, imagesc(reshape(Xiter_RS1(:,ind1), n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
% % [val2, ind2] = min(RelErr_RS2);
% % figure, imagesc(reshape(Xiter_RS2(:,ind2), n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
% % figure, imagesc(reshape(infofista.BestReg.X, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
% 
% % objFcnTrain = @(xx) leastsquares(A, xx, bn); 
% % options = optimoptions(@fmincon, 'Display', 'iter', 'SpecifyObjectiveGradient', true); % , 'Algorithm', 'sqp'
% % xm = fmincon(fun,zeros(size(x)),[],[],[],[],zeros(size(x)),ones(size(x)));
