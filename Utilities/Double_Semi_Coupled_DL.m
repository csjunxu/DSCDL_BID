% Main Function of Coupled Dictionary Learning
% Input:
% Alphap,Alphas: Initial sparse coefficient of two domains
% Xp    ,Xs    : Image Data Pairs of two domains
% Dp    ,Ds    : Initial Dictionaries
% Wp    ,Ws    : Initial Projection Matrix
% par          : Parameters 
%
%
% Output
% Alphap,Alphas: Output sparse coefficient of two domains
% Dp    ,Ds    : Output Coupled Dictionaries
% Up    ,Us    : Output Projection Matrix for Alpha
% 

function [Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, Up, Us, f] = Double_Semi_Coupled_DL(Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, par)

%% parameter setting

[dimX, numX]        =       size(Xp);
dimY                =       size(Alphap, 1);
numD                =       size(Dp, 2);
rho                 =       par.rho;
lambda1             =       par.lambda1;
lambda2             =       par.lambda2;
mu                  =       par.mu;
sqrtmu              =       sqrt(mu);
nu                  =       par.nu;
nIter               =       par.nIter;
t0                  =       par.t0;
epsilon             =       par.epsilon;
param.lambda        = 	    lambda1; % not more than 20 non-zeros coefficients
param.lambda2       =       lambda2;
param.mode          = 	    2;       % penalized formulation
param.approx=0;
param.K = par.K;
param.L = par.L;
f = 0;

%% Initialize Us, Up as I

Us = Ws; 
Up = Wp; 

%% Iteratively solve D A U

for t = 1 : nIter

    %% Updating Alphas and Alphap
    f_prev = f;
    Alphas = mexLasso([Xs;sqrtmu * Up * full(Alphap)], [Ds; sqrtmu * Us],param);
    Alphap = mexLasso([Xp;sqrtmu * Us * full(Alphas)], [Dp; sqrtmu * Up],param);
    dictSize = par.K;

    %% Updating Ds and Dp 
    for i=1:dictSize
       ai        =    Alphas(i,:);
       Y         =    Xs-Ds*Alphas+Ds(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Ds(:,i)    =    di;
    end

    for i=1:dictSize
       ai        =    Alphap(i,:);
       Y         =    Xp-Dp*Alphap+Dp(:,i)*ai;
       di        =    Y*ai';
       di        =    di./(norm(di,2) + eps);
       Dp(:,i)    =    di;
    end

    %% Updating Ws and Wp => Updating Us and Up
    Us = (1 - rho) * Us  + rho * Up * Alphap * Alphas' / ( Alphas * Alphas' + par.nu * eye(size(Alphas, 1)));
    Up = (1 - rho) * Up  + rho * Us * Alphas * Alphap' / ( Alphap * Alphap' + par.nu * eye(size(Alphap, 1)));
    Ws = Up /Us;
    Wp = Us /Up;

    %% Find if converge (NEED MODIFICATION)

    P1 = Xp - Dp * Alphap;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphap, 1);    
    P3 = Us * Alphas - Up * Alphap; %  20160725: Alphas - Wp * Alphap -> Us * Alphas - Up * Alphap
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Up, 'fro');
    fp = 1 / 2 * P1 + P2 + mu * (P3 + P4);
    
    P1 = Xs - Ds * Alphas;
    P1 = P1(:)'*P1(:) / 2;
    P2 = lambda1 *  norm(Alphas, 1);    
    P3 = Us * Alphas - Up * Alphap; %  20160725: Alphap - Ws * Alphas -> Us * Alphas - Up * Alphap
    P3 = P3(:)'*P3(:) / 2;
    P4 = nu * norm(Us, 'fro');  %%
    fs = 1 / 2 * P1 + P2 + mu * (P3 + P4);
    
    f = fp + fs;
	
        %% if converge then break
    if (abs(f_prev - f) / f < epsilon)
        break;
    end
    fprintf('Energy: %d\n',f);
    save tempDSCDL_BID_Dict Ds Dp Us Up Ws Wp par param;
end
    
