%% Cross-Sectional Dependence Tests
%
% LM (Lagrange Multiplier) test       >>> Breusch, T. S., & Pagan, A. R. (1980). The Lagrange multiplier test and its applications to model specification in econometrics. The Review of Economic Studies, 47(1), 239-253.
% CD test                             >>> Pesaran, M. H. 2004. General diagnostic tests for cross section dependence in panels. CESifo Working Paper Series No. 0435.
% and Pesaran, M. H. (2015). Testing weak cross-sectional dependence in large panels. Econometric Reviews, 34(6-10), 1089-1117.
% CD_LM test                          >>> Pesaran, M. H. 2004. "General diagnostic tests for cross section dependence in panels. CESifo Working Paper Series No. 0435.
% LM_adj test                         >>> Pesaran, M. H., Ullah, A., & Yamagata, T. (2008). A bias?adjusted LM test of error cross?section independence. The Econometrics Journal, 11(1), 105-127.
% LM_BC test (Fixed Effects)          >>> Baltagi, B. H., Feng, Q., & Kao, C. (2012). A Lagrange Multiplier test for cross-sectional dependence in a fixed effects panel data model. Journal of Econometrics, 170(1), 164-177.
%--------------------------------------------------------------------------
% Note: Pesaran (2015) establishes that the CD test proposed in Pesaran (2004) is best viewed as a test of weak cross-sectional dependence, namely the null of 
% weak dependence is more appropriate than the null of independence in the case of large panel data models where only pervasive cross-dependence is of concern.
%--------------------------------------------------------------------------
% This code is written by Amin Haghnejad, March 2021 (am.haghnejad@gmail.com)
%**************************************************************************
%%
clear
clc
%% === Read Data ===

data = xlsread('data.xlsx');          % An m-by-n matrix of observations with m = N*T and n = the number of all variables (dependent and regressors)
T = 49;                                   % Time dimension
N = 24;                                   % Cross-sectional dimension
alpha = 0.05;                             % Significance level
%%
Y = data(:,1);                            % Y denotes dependent variable
X = data(:,2:end);                        % X denotes a vector of regressors
n = size(X,2);                            % n denotes the number of regressors
k = n+1;                                  % k denotes the number of parameters (intercept and slope coefficients)
%
Z = [ones(N*T,1) X];
YY = cell(1,N);
ZZ = cell(1,N);
for i = 1:N
    YY{i} = Y(((1-T)+(T*i)):(T*i),:);
    ZZ{i} = Z(((1-T)+(T*i)):(T*i),:);
end
beta_ols = cell(1,N);
for i = 1:N
    beta_ols{i} = regress(YY{i},ZZ{i});   % The ordinary least squares (OLS) estimate of the individual regressions
end
U = cell(1,N);
for i = 1:N
    U{i} = YY{i}-(ZZ{i}*beta_ols{i});
end
co = zeros(N);
for i = 1:N-1
    for j = i+1:N
        co(i,j) = corr(U{i},U{j});
    end
end
roh = zeros(((N^2)-N)/2,1);
for i = 1:N-1
    for j = i+1:N
        roh(((i-1)*N)-((i*(i-1))/2)+(j-i)) = co(i,j);
    end
end
L = length(roh);
roh2 = zeros(L,1);
for i = 1:L
    roh2(i) = roh(i)^2;
end
%% LM test ---- Breusch and Pagan (1980)
%
LM = T*sum(roh2);
df_LM = (N*(N-1))/2;
pvalue_LM = 1-(chi2cdf(LM,df_LM));
CV_LM = chi2inv((1-alpha),df_LM);
%% Scaled LM (CD_LM) test ---- Pesaran (2004)
%
CDLM = ((1/(N*(N-1)))^(1/2))*sum((T*roh2)-1);
if CDLM > 0
    pvalue_CDLM = 2*(1-(normcdf(CDLM,0,1)));
else
    pvalue_CDLM = 2*(normcdf(CDLM,0,1));
end
UCV_CDLM = norminv((1-(alpha/2)),0,1);
LCV_CDLM = norminv((alpha/2),0,1);
%% CD test ---- Pesaran (2004, 2015)
%
CD = (((2*T)/(N*(N-1)))^(1/2))*sum(roh);
if CD > 0
    pvalue_CD = 2*(1-(normcdf(CD,0,1)));
else
    pvalue_CD = 2*(normcdf(CD,0,1));
end
%% Bias-adjusted LM (LM_adj) test ---- Pesaran, Ullah, and Yamagata (2008)
%
coo = zeros(N);
for i = 1:N-1
    for j = i+1:N
        coo(i,j) = (T-k)*co(i,j)^2;
    end
end
I = eye(T);
H = cell(1,N);
for i = 1:N
    H{i} = ZZ{i}*(((ZZ{i})'*(ZZ{i}))^-1)*(ZZ{i})';
end
M = cell(1,N);
for i = 1:N
    M{i} = I-H{i};
end
TM = zeros(1,N);
for i = 1:N
    TM(i) = trace(M{i});
end
MM = cell(N);
MM2 = cell(N);
for i = 1:N-1
    for j = i+1:N
        MM{i,j} = M{i}*M{j};
        MM2{i,j} = MM{i,j}^2;
    end
end
TMM = zeros(N);
TMM2 = zeros(N);
for i = 1:N-1
    for j = i+1:N
        TMM(i,j) = trace(MM{i,j});
        TMM2(i,j) = trace(MM2{i,j});
    end
end
E = zeros(N);
for i = 1:N-1
    for j = i+1:N
        E(i,j) = TMM(i,j)/(T-k);
    end
end
a2 = 3*((((T-k-8)*(T-k+2))+24)/((T-k+2)*(T-k-2)*(T-k-4)))^2;
a1 = a2-(1/(T-k)^2);
V = zeros(N);
for i = 1:N-1
    for j = i+1:N
        V(i,j) = (a1*((TMM(i,j))^2))+(2*a2*TMM2(i,j));
    end
end
cooo = zeros(N);
for i = 1:N-1
    for j = i+1:N
        cooo(i,j) = co(i,j)*((coo(i,j))/((V(i,j))^(1/2)));
    end
end
roh2a = zeros(((N^2)-N)/2,1);
for i = 1:N-1
    for j = i+1:N
        roh2a(((i-1)*N)-((i*(i-1))/2)+(j-i)) = cooo(i,j);
    end
end
LMadj = (((2*T)/(N*(N-1)))^(1/2))*sum(roh2a);
if LMadj > 0
    pvalue_LMadj = 2*(1-(normcdf(LMadj,0,1)));
else
    pvalue_LMadj = 2*(normcdf(LMadj,0,1));
end
%% Bias-corrected scaled LM (LM_BC) test ---- Baltagi, Feng, and Kao (2012)
%
IN = eye(N);
INT = eye(N*T);
IIT = ones(T,1);
D = kron(IN,IIT);
PD = D*((D'*D)^-1)*D';
QD = INT-PD;
W = [D X];
beta_fe = regress(Y,W);   % The fixed effects (FE) estimate of the parameters
UU = Y-((D*beta_fe(1:N,1))+(X*beta_fe(N+1:n+N,1)));
U = cell(1,N);
for i = 1:N
    U{i} = UU(((1-T)+(T*i)):(T*i),:);
end
co = zeros(N);
for i = 1:N-1
    for j = i+1:N
        co(i,j) = corr(U{i},U{j});
    end
end
roh = zeros(((N^2)-N)/2,1);
for i = 1:N-1
    for j = i+1:N
        roh(((i-1)*N)-((i*(i-1))/2)+(j-i)) = co(i,j);
    end
end
L = length(roh);
roh2 = zeros(L,1);
for i = 1:L
    roh2(i) = roh(i)^2;
end
LMBC = (((1/(N*(N-1)))^(1/2))*sum((T*roh2)-1))-(N/(2*(T-1)));
if LMBC > 0
    pvalue_LMBC = 2*(1-(normcdf(LMBC,0,1)));
else
    pvalue_LMBC = 2*(normcdf(LMBC,0,1));
end
%%
disp('  ------------------------------------------------------')
disp('The results of tests for cross-section dependence in panels')
disp('  ------------------------------------------------------')
fprintf('\n The number of time periods: %2.0f \n',T);
fprintf('\n The number of cross-sections: %2.0f \n',N);
fprintf('\n Null hypothesis: No cross-sectional dependence (weak dependence in the Pesaran''s (2004, 2015) CD test for large panels) \n');
disp(' ')
%
disp('  ------------------------------------------------------')
disp('  LM test ---- Breusch and Pagan (1980)')
disp('  ------------------------------------------------------')
fprintf('\n  Test Statistic \t\t d.f.  \t p-value \n\n');
fprintf('\n  %4.6f        \t %1.0f  \t %4.6f \n',LM,df_LM,pvalue_LM);
%%
disp('  ------------------------------------------------------')
disp('  Scaled LM (CD_LM) test ---- Pesaran (2004)')
disp('  ------------------------------------------------------')
fprintf('\n  Test Statistic \t\t\t p-value \n\n');
fprintf('\n  %4.6f \t\t\t\t\t %4.6f \n\n',CDLM,pvalue_CDLM);
%%
disp('  ------------------------------------------------------')
disp('  CD test ---- Pesaran (2004, 2015)')
disp('  ------------------------------------------------------')
fprintf('\n  Test Statistic \t\t\t p-value \n\n');
fprintf('\n  %4.6f \t\t\t\t\t %4.6f \n\n',CD,pvalue_CD);
%%
disp('  ----------------------------------------------------------------------------')
disp('  Bias-adjusted LM (LM_adj) test ---- Pesaran, Ullah, and Yamagata (2008)')
disp('  ----------------------------------------------------------------------------')
fprintf('\n  Test Statistic \t\t\t p-value \n\n');
fprintf('\n  %4.6f  \t\t\t\t %4.6f \n\n',LMadj,pvalue_LMadj);
%%
disp('  ----------------------------------------------------------------------------')
disp('  Bias-corrected scaled LM (LM_BC) test ---- Baltagi, Feng, and Kao (2012)')
disp('  ----------------------------------------------------------------------------')
fprintf('\n  Test Statistic \t\t\t p-value \n\n');
fprintf('\n  %4.6f \t\t\t\t\t %4.6f \n\n',LMBC,pvalue_LMBC);