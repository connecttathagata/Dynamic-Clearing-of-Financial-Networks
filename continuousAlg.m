function [time,V,A,c,varargout]=continuousAlg(dt0,T,c0,mu,sigma,dL,A0)

tol = 1e-12; %tolerance
%tol=1e-8
n = length(c0); %number of banks in the financial network

msize = round(1.3*T/dt0);
time = zeros(1,msize);
V = zeros(n,msize); %Wealth
A = zeros(n,n,msize); %Nominal Liabilities
c = zeros(n,msize); %Cash Flow
faroff = zeros(1,msize);
Z = randn(n,msize);

V0 = c0; %initial  condition
V(:,1) = V0;
if exist('A0','var')
    A(:,:,1) = A0;
else
    A(:,:,1) = dL(0)./repmat(max(sum(dL(0),2),1e-8),[1 n]);
end

ii = 1;
Lam = diag(V(:,1) < 0);
c(:,1) = c0;
faroff(1) = norm(V(:,1) - (c(:,1) + A(:,:,1).'*Lam*V(:,1)));
while time(ii) < T
    ii = ii+1;
    

    Lam0 = Lam+1;
    
    %make sure we hit dt0 time intervals if we subdivide intervals
    %previously
    dt=dt0-mod(time(ii-1)/dt0,1)*dt0;

    
    if ii-1 == msize

        mult = 1.1;
        time(msize+1:ceil(mult*msize)) = zeros(1,ceil(mult*msize)-msize);
        V(:,msize+1:ceil(mult*msize)) = zeros(n,ceil(mult*msize)-msize);
        A(:,:,msize+1:ceil(mult*msize)) = zeros(n,n,ceil(mult*msize)-msize);
        c(:,msize+1:ceil(mult*msize)) = zeros(n,ceil(mult*msize)-msize);
        faroff(msize+1:ceil(mult*msize)) = zeros(1,ceil(mult*msize)-msize);
        Z(:,msize+1:ceil(mult*msize)) = randn(n,ceil(mult*msize)-msize);
        msize = ceil(mult*msize);
    end

    drift = mu(time(ii-1),c(:,ii-1));
    diffusion = sigma(time(ii-1),c(:,ii-1))*Z(:,ii-1);
    while ~isequal(Lam,Lam0)
        Lam0 = Lam;

        mubar = (eye(n,n) - A(:,:,ii-1).'*Lam)\(drift -...
            sum(dL(time(ii-1)),1).' + A(:,:,ii-1).'*sum(dL(time(ii-1)),2));
        sigmabar = (eye(n,n) - A(:,:,ii-1).'*Lam)\(diffusion);
        %Checking solvency
        for bank = 1:n
            if sigmabar(bank)^2 - 4*mubar(bank)*V(bank,ii-1) >= 0
                if V(bank,ii-1) > tol
                    if mubar(bank) < -tol
                        dt = min(dt , ((-sigmabar(bank) - sqrt(sigmabar(bank)^2 - 4*mubar(bank)*V(bank,ii-1)))/(2*mubar(bank)))^2);
                    elseif (abs(mubar(bank)) < tol) && (sigmabar(bank) < -tol)
                        dt = min(dt , (V(bank,ii-1)/sigmabar(bank))^2);
                    end
                elseif V(bank,ii-1) < -tol
                    if abs(mubar(bank)) >= tol
                        dt = min(dt , ((-sigmabar(bank) + sqrt(sigmabar(bank)^2 - 4*mubar(bank)*V(bank,ii-1)))/(2*mubar(bank)))^2);
                    elseif (abs(mubar(bank)) < tol) && (sigmabar(bank) > tol)
                        dt = min(dt , (V(bank,ii-1)/sigmabar(bank))^2);
                    end
                end
            end
            if mubar(bank)*sigmabar(bank) < -tol
                dt = min(dt , (sigmabar(bank)/mubar(bank))^2);
            end
        end
        
        DELTA = mubar*dt + sigmabar*sqrt(dt);


            
        % Set of Recovering and Defaulting Banks
        Recovering = find((abs(V(:,ii-1)) < tol) & (DELTA > tol) & (diag(Lam) == 1)); %Were negative, trending upward at 0
        Defaulting = find((abs(V(:,ii-1)) < tol) & (DELTA < -tol) & (diag(Lam) == 0)); % Were positive, trending downward, hit 0
        
        % know Lambda changes for these (default -> positive and vice
        % versa)
        Lam((n+1)*(Recovering-1)+1) = 0; 
        Lam((n+1)*(Defaulting-1)+1) = 1;
        
    end

    time(ii) = time(ii-1) + dt;
    
    % Banks <0 need to use Gamma in A calculation; can't change Lam until
    % next step- in the ACTUAL continuous time model, Gamma = Lam
    % PURELY for numerical reasons
    Gamma = diag(V(:,ii-1) < -1e-8); 
    
    %update A, V and c
    A(:,:,ii)=Gamma*(A(:,:,ii-1)-diag(min(V(:,ii-1),-1e-8))\(dL(time(ii-1))-A(:,:,ii-1).*(dL(time(ii-1))*ones(n,n)))*dt) + ...
        (eye(n,n)-Gamma)*dL(time(ii-1))./repmat(max(sum(dL(time(ii-1)),2),1e-8),[1 n]);

    V(:,ii) = V(:,ii-1) + mubar*dt + sigmabar*sqrt(dt); %If true dc term than remove the *dt
    c(:,ii) = c(:,ii-1) + drift*dt + diffusion*sqrt(dt);
    faroff(ii) = norm(V(:,ii) - (c(:,ii) + A(:,:,ii).'*Lam*V(:,ii)));
    
end
time = time(1:ii);
V = V(:,1:ii);
A = A(:,:,1:ii);
c = c(:,1:ii);
varargout{1} = faroff(1:ii);



end