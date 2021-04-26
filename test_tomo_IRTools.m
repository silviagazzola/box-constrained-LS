clear, close all

% n = 128;
n = 256;
angles = 0:2:120;
% angles = 0:1:90;
% angles = 0:1:70;
phantomImage = 'threephasessmooth'; % 'grains';
ProbOptions = PRtomo('defaults');
ProbOptions = PRset(ProbOptions, 'angles', angles, 'phantomImage', phantomImage);
[A, b_true, x_true, ProbInfo] = PRtomo(n, ProbOptions);
NoiseLevel = 1e-2;
b = PRnoise(b_true, NoiseLevel);
nl = NoiseLevel;
b0 = A'* b; bn = b; x = x_true;
b0(b0<0)=0;
m_in = 20; m_out = 20;
disp('NN-FCGLS')
stopC_1.kind = 'discr';
stopC_1.noise = nl;
[Xiter_RS1, RelErr_RS1, RelRes_RS1, Alpha_RS1, tot_it1, cycle_it1, stopC_it1] = fcgls_rp_nn_ReSt (A, bn, m_in, m_out, x, b0, 1, stopC_1);
disp('BC-FCGLS')
stopC_1.kind = 'discr';
stopC_1.noise = nl;
b0(b0>1)=1;
m_out = 200;
[Xiter_RS2, RelErr_RS2, RelRes_RS2, Alpha_RS2, tot_it2, cycle_it2, stopC_it2] = fcgls_bc_ReSt (A, bn, 0, 1, m_in, m_out, x, b0, 1, stopC_1);
%
optfista.NoStop = 'on';
optfista.xMax = 1;
optfista.x_true = x;
optfista.MaxIter = 200;
[Xfista, infofista] = IRfista(A, bn, optfista);

optcg.NoStop = 'on';
optcg.x_true = x;
[Xcg, infocg] = IRcgls(A, bn, optcg);

%
figure, semilogy(infocg.Enrm, 'LineWidth', 2)
hold on
semilogy(RelErr_RS1, 'LineWidth', 2)
semilogy(RelErr_RS2, 'LineWidth', 2)
semilogy(infofista.Enrm, 'LineWidth', 2)
legend('CGLS', 'NN-FCGLS', 'box-FCGLS', 'FISTA')
set(gca,'FontSize',22)

figure, semilogy(infocg.Rnrm)
hold on
semilogy(RelRes_RS1)
semilogy(RelRes_RS2)
semilogy(infofista.Rnrm)
legend('CGLS', 'NN-FCGLS', 'box-FCGLS', 'FISTA')
set(gca,'FontSize',20)

% % figure, imagesc(reshape(Xiter_RS1(:,end), n, n)), axis image, axis off
% figure, imagesc(reshape(Xiter_RS2(:,end), n, n)), axis image, axis off

% figure, imagesc(reshape(x, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
figure, imagesc(reshape(infocg.BestReg.X, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
[val1, ind1] = min(RelErr_RS1);
figure, imagesc(reshape(Xiter_RS1(:,ind1), n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
[val2, ind2] = min(RelErr_RS2);
figure, imagesc(reshape(Xiter_RS2(:,ind2), n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)
figure, imagesc(reshape(infofista.BestReg.X, n, n)), axis image, axis off, colorbar, set(gca,'FontSize',20)


% fun = @(xx)norm(A*xx - bn);
% options = optimoptions(@fmincon, 'Display', 'iter', 'SpecifyObjectiveGradient', true); % , 'Algorithm', 'sqp'
% xm = fmincon(fun,zeros(size(x)),[],[],[],[],0,1);

% % plotting a bit more
% figure, imagesc(linspace(0,120,61),1:362,reshape(bn1, ProbInfo1.bSize))
% set(gca,'ytick',[])
% set(gca,'xtick',0:10:120)
% set(gca,'FontSize',18)
% figure, imagesc(linspace(0,90,91),1:362,reshape(bn, ProbInfo.bSize))
% set(gca,'ytick',[])
% set(gca,'xtick',0:10:90)
% set(gca,'FontSize',18)
% figure, imagesc(reshape(x, n, n)), axis image, axis off, colorbar
% set(gca,'FontSize',18)

% % small-scale test for matlab
% n = 16;
% angles = 0:2:90;
% phantomImage = 'threephasessmooth'; % 'grains';
% ProbOptions = PRtomo('defaults');
% ProbOptions = PRset(ProbOptions, 'angles', angles, 'phantomImage', phantomImage);
% [A, b_true, x_true, ProbInfo] = PRtomo(n, ProbOptions);
% NoiseLevel = 1e-2;
% b = PRnoise(b_true, NoiseLevel);
% A = full(A);
% objFcnTrain = @(xx) leastsquares(A, xx, b);
% options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp'); % 'SpecifyObjectiveGradient', true %
% xm = fmincon(fun,zeros(size(x_true)),[],[],[],[],zeros(size(x_true)),ones(size(x_true)));
% objFcnTrain = @(xx) leastsquares(A, xx, b);
% options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp'); % 'SpecifyObjectiveGradient', true %
% xm = fmincon(objFcnTrain,zeros(size(x_true)),[],[],[],[],zeros(size(x_true)),ones(size(x_true)));
% objFcnTrain = @(xx) myleastsquares(A, xx, b);
% options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp'); % 'SpecifyObjectiveGradient', true %
% xm = fmincon(objFcnTrain,zeros(size(x_true)),[],[],[],[],zeros(size(x_true)),ones(size(x_true)));
% objFcnTrain = @(xx) myleastsquares(A, xx, b);
% options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true); %
% xm = fmincon(objFcnTrain,zeros(size(x_true)),[],[],[],[],zeros(size(x_true)),ones(size(x_true)));
% objFcnTrain = @(xx) myleastsquares(A, xx, b);
% options = optimoptions(@fmincon, 'Display', 'iter', 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true); %
% xm = fmincon(objFcnTrain,0.5*ones(size(x_true)),[],[],[],[],zeros(size(x_true)),ones(size(x_true)));
% objFcnTrain = @(xx) myleastsquares(A, xx, b);
% options = optimoptions(@fmincon, 'Display', 'iter'); %, 'Algorithm', 'sqp', 'SpecifyObjectiveGradient', true); %
% xm = fmincon(objFcnTrain,0.5*ones(size(x_true)),[],[],[],[],zeros(size(x_true)),ones(size(x_true)));





