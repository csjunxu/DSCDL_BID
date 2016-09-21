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

function [Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, Up, Us, Ps, f] = ADPU_Double_Semi_Coupled_DL(Alphap, Alphas, Xp, Xs, Dp, Ds, Wp, Ws, par)

%% parameter setting
param.lambda        = 	    par.lambda1; % not more than 20 non-zeros coefficients
param.lambda2       =       par.lambda2;
param.mode          = 	    2;       % penalized formulation
param.approx=0;
param.K = par.K;
param.L = par.L;
f = 0;

%% Initialize Us, Up as I

Us = Ws; 
Up = Wp; 

%% Initialize Ps as 0 matrix
Ps = zeros(size(Xs));

%% Iteratively solve D A U

for t = 1 : par.nIter
    
    %% Updating Alphas and Alphap
    f_prev = f;
    Alphas = mexLasso([Xs - Ps; par.sqrtmu * Up * full(Alphap)], [Ds; par.sqrtmu * Us],param);
    Alphap = mexLasso([Xp; par.sqrtmu * Us * full(Alphas)], [Dp; par.sqrtmu * Up],param);
    dictSize = par.K;

    %% Updating Ds and Dp 
    for i=1:dictSize
       ai        =    Alphas(i,:);
       Y         =    Xs - Ps - Ds * Alphas + Ds(:,i) * ai;
       di        =    Y * ai';
       di        =    di ./ (norm(di,2) + eps);
       Ds(:,i)   =    di;
    end

    for i=1:dictSize
       ai        =    Alphap(i,:);
       Y         =    Xp - Dp * Alphap + Dp(:,i) * ai;
       di        =    Y * ai';
       di        =    di ./ (norm(di,2) + eps);
       Dp(:,i)  =    di;
    end
    
    %% Updating Ps
    Ps = (Xs - Ds * Alphas) / (1 + par.nup);

    %% Updating Ws and Wp => Updating Us and Up
    Us = (1 - par.rho) * Us  + par.rho * Up * Alphap * Alphas' / ( Alphas * Alphas' + par.nu * eye(size(Alphas, 1)));
    Up = (1 - par.rho) * Up  + par.rho * Us * Alphas * Alphap' / ( Alphap * Alphap' + par.nu * eye(size(Alphap, 1)));
    Ws = Up /Us;
    Wp = Us /Up;

    %% Find if converge (NEED MODIFICATION)
    P1 = Xp - Dp * Alphap;
    P1 = P1(:)'*P1(:) / 2;
    P2 = par.lambda1 *  norm(Alphap, 1);    
    P3 = Us * Alphas - Up * Alphap; %  20160725: Alphas - Wp * Alphap -> Us * Alphas - Up * Alphap
    P3 = P3(:)'*P3(:) / 2;
    P4 = par.nu * norm(Up, 'fro');
    fp = 1 / 2 * P1 + P2 + par.mu * (P3 + P4);
    
    P1 = Xs - Ds * Alphas - Ps;
    P1 = P1(:)'*P1(:) / 2;
    P2 = par.lambda1 *  norm(Alphas, 1);    
    P3 = Us * Alphas - Up * Alphap; %  20160725: Alphap - Ws * Alphas -> Us * Alphas - Up * Alphap
    P3 = P3(:)'*P3(:) / 2;
    P4 = par.nu * norm(Us, 'fro');  
    P5 = par.nup * norm(Ps, 'fro'); % 20160815: added by Jun Xu 
    fs = 1 / 2 * P1 + P2 + par.mu * (P3 + P4 + P5);
    
    f = fp + fs;
	
    %% if converge then break
    if (abs(f_prev - f) / f < par.epsilon)
        break;
    end
    fprintf('Energy: %d\n',f);
    save tempDSCDL_BID_Dict Ds Dp Us Up Ws Wp par param;
end
    
