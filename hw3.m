% Shivi Kalra
% Homework 3 Macro 696G
beta = 0.9932;
sigma =1.5;
b=0.5;
ys = [1,b];
P = [0.97, 0.03; 0.5,0.5]; %transition matrix
%PStar = (P(2,1)/(1-P(1,1)+P(2,1))....
    %+(1-P(2,1))/(1-P(1,1)+P(2,1))); %invariant distribution
%asset vector
alow = -2;
ahigh = 5;
numa = 10;
a = linspace(alow,ahigh, numa);
%Initial guess for q
qmin = 0.99;
qmax = 1;
qguess = (qmin+qmax)/2;
% Iterate over asset prices
aggsavings =1;
while abs(aggsavings)>=0.01;
    consp = bsxfun(@minus,a',qguess*a);
    consp = bsxfun(@plus,consp,permute(ys,[1 3 2]));
    ret =(consp.^(1-sigma))./(1-sigma);
    ret(consp<0)=-Inf;
    
end


    % Initial value fn guess
    v_guess = zeros(2,numa);
    
    % Value fn Iteration
    v_tol = 1;
    while v_tol>0.0001;
        %construct total return fn
        v_mat = ret + beta*....
            repmat(permute(P*v_guess,(321)),[numa,1]);
        
     % choose highest value
     [vfn,pol_indx] = max(v_mat,[],2);
     v_tol = abs(max(v_guess(:)-vfn(:)));
     v_guess = vfn; % update value fn
    end;
    
    % Keep decision rule
    polfn = a(pol_indx);
    
    % Set up initial distribution
     Mu = zeros (z,numa);
     Mu(1,11) =1; % initial guess; everyone employed and zero assets
     
     % Iterate over distributions 
     
     mu_tol = 1;
     while mu_tol > 1e-0.8;
         [emp_ind, a_ind, mass] = find (Mu>0);
         Munew = zeros (size(Mu));
         for ii =1:length(emp_ind)
             apr_ind = pol_indx(emp_ind(ii),a_ind(ii));
             Munew(:,aprind)= Munew(:,aprind)+....
                 (P(emp_ind(ii),:)* Mu(emp_ind(ii),a_ind(ii)))';
         end
         
     end
     mu_tol = max (abs(Munew(:)- Mu(:)));
     Mu = Munew;
     
     
     
     
    
    
    



