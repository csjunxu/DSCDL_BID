% May 20, 2006 (last revised on December, 2007)
%
% This is a Matlab implementation of a coordinate gradient descent method
% for solving the nonsmooth optimization problem
%       min   f(x) + c||x||_1 ,
% where f is a real-valued smooth function of n real variables.
% The method, for a given x, computes a direction d as the solution of the subproblem
%       min   g'*d + d'*H*d/2 + c||x+d||_1   subject to   d(j)=0 for j not in J,
% where g is the gradient of f at x, H is a positive definite diagonal matrix,
% and J is a nonempty subset of {1,...,n}.  Then x is updated by
%       x = x + step*d,
% and this is repeated until a termination criterion is met.
% The index subset J is chosen by a Gauss-southwell-q rule (in the M-file dirq.m),
% and step is chosen by an Armijo stepsize rule.
%
% This code comes with three M-files:  dirq.m,  nz.m,  signx.m
% The user must supply the following M-files:
%      fnc.m (f-value)
%      grad.m (f-gradient value)
%      hessian.m (f-Hessian diagonals or approximation of, to be used in constructing H)
% The user should specify the initial x.
% The user may also change the termination tolerance tol from its default value of 1e-4.
% If the f-Hessian diagonals are unknown or expensive to compute, then you can set the diagonals
% to any positive constants.   The CGD method converges faster if the diagonals are close to the
% diagonals of the Hessian of f.
%
% This Matlab code is written by Paul Tseng and Sangwoon Yun.
% Reference: A Coordinate Gradient Descent Method for Nonsmooth Separable Minimization,
% Department of Mathematics, University of Washington, Seattle, June 2006; to appear in
% Mathematical Programming, Series B.

function a=cgdsq(x,D,a,param)

% clear;
% c = input('Input regularization parameter c > 0 ');    % input regularization parameter
c = param.lambda;
% %load x.dat;              % initial point
% x=ones(1000,1);           % initial point

tol=1e-4;		        % Termination tolerance (for diag. scaled residual maxRr).

n=size(a,1);              % number of variables
g=grad(a);                % f-gradient at x
objf=fnc(x-D*a);
a1 = a(1:size(a,1)/2,:);
a2 = a(size(a,1)/2+1:end,:);
temp_a = [a1 a2];
objl12 = norm_l12(temp_a);          % l1-norm
obj=objf+c*objl12;         % f-value at x
objval(1)=obj;		  % Store objective values for graph plot

k=0;                      % iteration count
cgdk=0;                   % iteration count for cgd
lbfgsk=0;                 % iteration count for L-BFGS
rank1k=0;                 % iteration count for rank-1 approx.

step=1;  	              % initial stepsize for coord. gradient iterations.
ups=.5;                   % initial threshold for choosing J.

ml=5;			        % L-BFGS memory size.
S=zeros(n,ml);		  % S stores the last ml changes in x
Y=zeros(n,ml);		  % Y stores the last ml changes in gradient of f.
rho=zeros(1,ml);          % rho stores the last ml   s'*y.
kl=0;			        % number of nonempty columns in S & Y.

% Parameters for Armijo stepsize rule
sigma = .1;
beta = .5;
gamma = 0;		        % parameter for dirderiv. in cgd iterations.
gamma2 = 0;		        % parameter for dirderiv. in rank-1 Hessian accel.
sigmal = .1;		  % Armijo parameter for rank-1 Hessian accel.

r0=norm(median([ a' ; (g'+c) ; (g'-c) ])); % initial residual (diagnostic)

stop = 0;		        % flag for termination

t1=cputime;		        % CPU time for the algorithm.

while (stop==0) && k< 5000 % original 50000
    
    % run cgd for first 10 iterations
    if k<10
        flagl=1;
    end
    
    %t2=clock;
    
    % if lower bd is small (e.g., 1e-8), roundoff error occurs on BT function.
    h=min(max(hessian(a),1e-2),1e9);	% positive approx of Hessian diagonal
    
    if flagl==1
        
        [maxRr,d]=dirq(c,a,g,h,ups);	% compute cgd direction
        
        %normd(1)=maxRr;			% Store residual norm for graph plot
        
        % Check for termination
        if (maxRr <= tol)
            stop=1;  break;
        end
        
        % Armijo stepsize rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        step = min(step/beta,1);
        newobjf = fnc(x-D*(a+step*d));
        nsx = sign(a+step*d);
        psx = sign(a);
        ad = a+d;
        a1 = ad(1:size(a,1)/2,:);
        a2 = ad(size(a,1)/2+1:end,:);
        temp_a = [a1 a2];
        dirderiv = g'*d+gamma*h*(d.^2)+c*(norm_l12(temp_a)-objl12);
        nf = 1;
        while (newobjf-objf)/step+c*(ones(1,n)*((nsx-psx).*a+step*nsx.*d))/step > sigma*dirderiv
            step = step*beta;
            newobjf = fnc(x-D*(a+step*d));
            nsx = sign(a+step*d);
            nf = nf+1;
            
            if  (step<1e-30)
%                 fprintf('CGD stepsize small! Check roundoff error in newobj-obj.\n');
                stop=1; break;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        spit(k+1)=nf;
        s = step*d;		%save change in x for L-BFGS.
        a = a + s;
        objf = newobjf;
        a1 = a(1:size(a,1)/2,:);
        a2 = a(size(a,1)/2+1:end,:);
        temp_a = [a1 a2];
        objl12 = norm_l12(temp_a);          % l1,2-norm
        obj=objf+c*objl12;     
        
        % Update the threshold for choosing J, based on the current stepsize.
        if step >.001
            ups=max(1e-4,ups/10);
        else
            if step < .000001
                ups=min(.9,50*ups);
                step=.001;
            end
        end
        
        if  (mod(k,100)<50) && (k>=10)
            flagl=0;		% run BFGS acceleration
        else
            flagl=1;		% run CGD iteration
        end
        
        cgdk=cgdk+1;
        
        %pause;
        
    else 		% L-BFGS acceleration iteration
        
        %compute sigx using the Hessian diagonals
        [sigx,t,maxRr]=signx(c,a,g,h);
        
        if (maxRr <= tol)
            stop=1;  break;
        end
        
        asigx=abs(sigx);
        q=asigx.*g+c*sigx;   	%q is the obj gradient w.r.t. x_i "far from 0".
        
        if norm(q)>tol
            
            if kl==0
                d = -q;
            else
                ql = q;
                for i = 1 : kl
                    alpha(i) = S(:,i)'*ql/rho(i);
                    ql = ql - alpha(i)*Y(:,i);
                end
                r = ql*max(1e-6,rho(1)/(Y(:,1)'*Y(:,1)));
                for i = kl: -1: 1
                    betal = Y(:,i)'*r/rho(i);
                    r = r + S(:,i)*(alpha(i)-betal);
                end
                d=-r.*asigx;    %|x_i|>t whenever d_i not=0
                
                % steepest descent safeguard to ensure convergence
                if q'*d>-1e-20*(q'*q) || norm(d)>1e6*norm(q)
                    d = -q;
%                     fprintf('steepest desc. safeguard');
                end
            end
            
            % Armijo stepsize rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            stepl=1;
            newobjf = fnc(x-D*(a+stepl*d));
            nsx = sign(a+stepl*d);
            psx = sign(a);
            dirderiv = q'*d;
            nfl = 1;
            while (newobjf-objf)/stepl+c*(ones(1,n)*((nsx-psx).*a+stepl*nsx.*d))/stepl > sigma*dirderiv
                stepl = stepl*beta;
                newobjf = fnc(x-D*(a+stepl*d));
                nsx = sign(a+stepl*d);
                
                if  (stepl<1e-30)
%                     fprintf('L-BFGS stepsize small! Check roundoff error in newobj-obj.\n');
                    stop=1;
                    break;
                end
                
                % check the rounding error
                %lbfgscheck=abs(((newobj-obj)/stepl)-q'*d)
                %pause;
                
                nfl = nfl+1;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            s = stepl*d;
            a = a + s;
            objf = newobjf;
            a1 = a(1:size(a,1)/2,:);
            a2 = a(size(a,1)/2+1:end,:);
            temp_a = [a1 a2];
            objl12 = norm_l12(temp_a);          % l1,2-norm
            obj=objf+c*objl12;   
            
            % Check if objective is further improved by zeroing out those x_i with sigx_i=0.
            % (This did not seem to help).
            %       if c>0 & hessian < tol*1e3
            %         newx = abs(sigx).*x;
            %         newobj = fnc(c,newx)
            %         nfl = nfl+1;
            %         if (newobj<obj)
            %           s = newx-x;
            %           x = newx;
            %           obj=newobj;
            %         end
            %       end
            
            spitl(k+1)=nfl;
            
            if (mod(k,100)<50) %| stepl> 0.1
                flagl=0;		% run L-BFGS acceleration
            else
                flagl=1;		% run CGD iteration
            end
            
            lbfgsk=lbfgsk+1;
            
        else
            
            flagl=1; 		% run CGD iteration
            
        end
        
    end
    
    %time=etime(clock,t2);
    %fprintf(' Armijo time= %g  \n',time);
    %t2=clock;
    
    gold=g;
    g=grad(a);
    y=g-gold;		%save change in g for L-BFGS.
    
    %time=etime(clock,t2);
    %fprintf(' g time= %g  \n',time);
    
    %Check to save s and y for L-BFGS.
    if y'*y>1e-40
        sy=s'*y;
        gammal = sy/(y'*y);		%gammal = O(1/max_eigenvalue(Hessian))
        if gammal*max(h) > 1e-10
            kl = min(ml,kl+1);
            S = [s, S];
            Y = [y, Y];
            rho = [sy, rho];
            S = S(:,1:ml);
            Y = Y(:,1:ml);
            rho = rho(:,1:ml);
        end
    end
    
    %time=etime(clock,t2);
    %fprintf(' resid time= %g  \n',time);
    %pause
    
    % Acceleration step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (flagl==1)&&(mod(k,10)==1)&&(kl>0)
        
        % Use a rank-1 approx. of Hessian of f in quadratic model
        % and minimize with respect to all coordinates:
        %     min  g'd + |b'd|^2/2 + c||x+d||_1
        % The quadratic Hessian bb' is only psd, so it may have no minimum.
        
        % Choose b to satisfy (bb')s = y
        b=Y(:,1)/sqrt(rho(1));
        
        % If there exists i with b_i=0 and |g_i|>c, then subproblem has no optimal soln
        
        sgc = find(c>=abs(g));
        absb = abs(b);
        maxb = max(absb);
        absb(sgc) = maxb;
        if min( absb ) > maxb*1e-6
            
            %    d = -x + db_i*e_i,
            % where e_i is ith unit coordinate vector, db_i solves
            %    qb_i =  min  (g_i - (b'x)b_i)db_i + b_i^2 db_i^2/2 + c|db_i|,
            % and i is chosen to minimize qb_i .
            
            b(find(abs(b)<=0))=1e-20;  %perturb the zero entries of b to be nonzero
            bt=b';
            gb=(g-(bt*a)*b)';
            db=-median([ zeros(1,n) ; (gb+c)./(bt.^2) ; (gb-c)./(bt.^2) ]);
            qb=gb.*db + (bt.*db).^2/2 + c*abs(db);
            [minqb,index]=min(qb);
            d=-a;
            d(index)=d(index)+db(index);
            
            
            ad = a+d;
            a1 = ad(1:size(a,1)/2,:);
            a2 = ad(size(a,1)/2+1:end,:);
            temp_a = [a1 a2];
            dirderiv = g'*d+gamma2*(bt*d)^2+c*(norm_l12(temp_a)-objl12);
            if dirderiv < (1e-8)*norm(d)
                
                %if the direction is descent, increase the iteration and check the termination
                k=k+1;
                normd(k)=maxRr;
                objval(k+1)=obj;
                %         h=min(max(hessian(x),1e-12),1e12);
                %         R=-median([ x' ; (g'+c)./h ; (g'-c)./h ]);
                %         absR=abs(R);
                %         maxRr=max(h.*absR);
                R=-median([ a' ; (g'+c) ; (g'-c) ]);
                absR=abs(R);
                maxRr=max(absR);
                
                if (maxRr <= tol)
                    stop=1;  break;
                end
                
                % Armijo rule for acceleration step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                sp = 1;
                newobjf = fnc(x-D*(a+sp*d));
                nsx = sign(a+sp*d);
                psx = sign(a);
                nf = 1;
                while (newobjf-objf)/sp+c*(ones(1,n)*((nsx-psx).*a+sp*nsx.*d))/sp > sigma*dirderiv
                    sp = sp*beta;
                    newobjf = fnc(x-D*(a+sp*d));
                    nsx = sign(a+sp*d);
                    nf = nf+1;
                    
                    if  (sp<1e-30)
%                         fprintf('rank-1 accel. stepsize small!  Check roundoff error in newobj-obj.\n');
                        stop=1;
                        break;
                    end
                    
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                spit(k+1)=nf;
                
                s = sp*d;           %save change in x for L-BFGS.
                a = a + s;
                objf = newobjf;
                objl1 = norm(a,1);
                a1 = a(1:size(a,1)/2,:);
                a2 = a(size(a,1)/2+1:end,:);
                temp_a = [a1 a2];
                objl12 = norm_l12(temp_a);          % l1,2-norm
                obj=objf+c*objl12;   
                
                gold=g;
                g=grad(a);
                y=g-gold;		%save change in g for L-BFGS.
                
                %Check to save s and y for L-BFGS.
                if y'*y>1e-40
                    sy=s'*y;
                    gammal = sy/(y'*y);
                    if gammal*max(h) > 1e-10		%h is the Hessian diagonal at previous x.
                        kl = min(ml,kl+1);
                        S = [s, S];
                        Y = [y, Y];
                        rho = [sy, rho];
                        S = S(:,1:ml);
                        Y = Y(:,1:ml);
                        rho = rho(:,1:ml);
                    end
                end
                rank1k=rank1k+1;
            end
        end
        
        %pause;
        
    end
    
    k=k+1;
    normd(k)=maxRr;
    objval(k+1)=obj;
    
    if mod(k,500)==0
        
        r1=norm(median([ a' ; (g'+c) ; (g'-c) ]))/r0;   %diagnostic
        
%         fprintf('obj= %g  maxRr= %g r= %g CGD step= %g  LBFGS step= %g  #iter= %g  #nzr= %g  t= %g \n',obj,full(maxRr),full(r1),step,stepl,k,sum(full(asigx)),full(t));
        %pause;
        
    end
    
    
    
end

nonz=nz(a,n);           %number of nonzero components of the computed stationary point

% fprintf(' iter= %g cgditer= %g lbfgsiter= %g rank1iter= %g ob=%g maxRr= %g cpu= %g nonzero=%g \n',k,cgdk,lbfgsk,rank1k,objval(end),full(maxRr),cputime-t1,nonz);