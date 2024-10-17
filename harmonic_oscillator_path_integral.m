clear
% A simulation of a 1-dimensional quantum harmonic oscillator is done
% through a Monte Carlo- Metropolis algoritm with 3 possible approximated
% actions. The coorelation function is measured and, from that, the first
% energy gap il plotted.

[acceptance_ratio,elapsed_time]=main;

fprintf(['Total elapsed time is ',num2str(elapsed_time),' seconds. \n',...
         'The acceptance ratio for the metropolis algoritm is ',num2str(acceptance_ratio),' (should be 0.5 ± 0.1). \n'])

function [acceptance_ratio,elapsed_time]=main(setup)

    entries={      'N (Lattice size)',...
                   'N_corr (Number of sweeps before measurement) ',...
                   'N_cf (Number of measurements)',...                  
                   'epsylon (Epsylon of randomness for the update)',...
                   'a (Lattice spacing)',...
                   'mode (1=Approximation of 1st derivative, 2=Apporximation of 2nd derivative, 3=Improved approximation of 2nd derivative, 4=All simulations)'};
    if(nargin==0)
        default={  '20',...
                   '20',...
                   '1000',...
                   '1.4',...
                   '0.5',...
                   '1'};     
    else
        default=setup;        
    end
  
    addopt.Resize='on';
    addopt.WindowStyle='normal';
    addopt.Interpreter = 'none';

    header = 'Harmonic oscillator simulation';
    lines = 1;
    setup = inputdlg(entries, header, lines, default, addopt);
    
    N=str2num(setup{1});            % Lattice size
    N_corr=str2num(setup{2});       % Number of sweeps before measurement
    therm = N_corr*5 ;              % Thermalization steps
    epsylon= str2num(setup{4});     % Epsylon of randomness for the update             
    N_cf=str2num(setup{3});         % Number of measurements
    a=str2num(setup{5});            % Lattice spacing       
    mode=str2num(setup{6});         % 1=Approximation of 1st derivative, 2=Apporximation of 2nd derivative, 3=Improved approximation of 2nd derivative, 4=All simulations

    
    tic;
    if mode==4
        plotted=1:3;
        for i=1:3
        [d_E,acceptance_ratio]=simulation(therm,N,epsylon,a,i,N_corr,N_cf);
        [teo,plotted(i)]=plotting(d_E,N,a);
        end
        legend([teo plotted],{'Exact result','Approximation of 1st derivative','Apporximation of 2nd derivative','Improved approximation of 2nd derivative'})
    else
        [d_E,acceptance_ratio]=simulation(therm,N,epsylon,a,mode,N_corr,N_cf);
        [teo,plotted]=plotting(d_E,N,a);
        legend([teo plotted],{'Exact result','Simulation'})
    end

    elapsed_time=toc;
end

function [d_E,acceptance_ratio]=simulation(therm,N,epsylon,a,mode,N_corr,N_cf)
    % Function that does the entire simulation
    % It returns the energy gap and the acceptance_ratio
    
    % therm is the number of thermalization steps
    % N is the lattice size
    % epsylon is the randomness for the update
    % a is the lattice spacing
    % N_corr is the number of correlation steps skipped
    % N_cf is the number of measurements

    par=0;                          % Number of metropolis steps accepted
    tot=0;                          % Total steps
  
    x  = zeros(N,1);		        % Cold start
    G =  zeros(N_cf,N+1);           % Collector of all correlation function
    
    for j = 1:therm             
        [x,par,tot]=update_x(N,x,epsylon,a,par,tot,mode);
    end
    % Simulation and measurement
    for alpha= 1:N_cf         
         for j=1 : N_corr                
            [x,par,tot]=update_x(N,x,epsylon,a,par,tot,mode);    
         end
            
         for k=0:N
            G(alpha,k+1)=compute_G(x,k,N);
         end       
    end
   
    avg_G=sum(G(:,:))/N_cf;         % Averege of G
    acceptance_ratio = 1- par/tot;                     
    d_E= delta_E(avg_G,a);

end

function [teo,plotted]=plotting(d_E,N,a)
    % Function that plots the energy gap and the teoretical prediction
    % Also sets the options for the figure
    % It returns the plot's labels

    %d_E are the energy gaps measured
    % N is the lattice size
    % a is the lattice spacing

    xx=1:1:N;
    xx=xx*a;
    plotted=plot(xx,d_E','*');
    hold on;
    teo=plot(linspace(0,N),linspace(0,N)*0+1,'k');
    
    ylim([0, 2]);
    xlim([0, N/5]);

end

function dE=delta_E(G1,a)   
    % Function that calculates the energy gap from the correlation funcion
    % It returns the energy gap

    % G1 is the correlation function
    % a is the lattice spacing

    G2=G1;
    G1(end)=[];
    G2(1)=[];
    dE=log(abs(G1./G2))/a;
end

function [x,partial,total]=update_x(N,x,epsylon,a,partial,total,mode)  
    % Function that updates the lattice through the metropolis algoritm
    % It returns the updated lattice and the acceptance ratio
    
    % x is the lattice to update
    % N is the lattice size
    % epsylon is the randomness for the update
    % a is the lattice spacing
    % partial and total are for the acceptance ratio
    % mode is to choose the action to use

    x_old=x;
   
    % Metropolis algoritm
    for k = 1:N
        x(k)=x_old(k) + epsylon*(2*rand-1);
        dS=delta_S(x,x_old,k,a,N,mode);
        if dS>0 && exp(-dS)<rand  
            x(k)=x_old(k);
            partial=partial+1;
        end
        total=total+1;
    end
end

function dS=delta_S(x,x_old,k,a,N,mode)
    % Function that calculates the difference in action brought by a sweep
    % It returns the delta of energy

    % x is the updated lattice
    % x_old is the lattice before the update
    % k is the position of the point in the lattice
    % N is the lattice size
    % a is the lattice spacing
    % mode is to choose the action to use

    if mode==1
    S       = 0.5*a*x(k)*x(k)         + x(k)*(x(k) - x(mod(k-2,N)+1) - x(mod(k,N)+1))/a;
    S_old   = 0.5*a*x_old(k)*x_old(k) + x_old(k)*(x_old(k) - x_old(mod(k-2,N)+1) - x_old(mod(k,N)+1))/a;
    dS = S-S_old;

    elseif mode==2
    A=eye(N);
    D2=(eye(N,N)*(-2)+circshift(eye(N,N), [0, 1])+circshift(eye(N,N), [0, -1]))/(a^2);
    z1=A(:,k);
    z2=z1+A(:,mod(k-2,N)+1)+A(:,mod(k,N)+1);
                              
    S     = (-0.5*a)*(z2.*x)'*(D2*x)         + (a*0.5)*(z1.*x)'*x;
    S_old = (-0.5*a)*(z2.*x_old)'*(D2*x_old) + (a*0.5)*(z1.*x_old)'*x_old;
    dS = S-S_old;

    elseif mode==3
    A=eye(N);
    D2=(eye(N,N)*(-2)+circshift(eye(N,N), [0, 1])+circshift(eye(N,N), [0, -1]))/(a^2);
    z1=A(:,k);
    z2=z1+A(:,mod(k-2,N)+1)+A(:,mod(k,N)+1);
    z3=z2+A(:,mod(k-3,N)+1)+A(:,mod(k+1,N)+1);
                
    S     = (-0.5*a)*(z3.*x)'*(D2-(a^2/12)*D2^2)*x         + (a*0.5)*(z1.*x)'*x*(1+a^2/12);
    S_old = (-0.5*a)*(z3.*x_old)'*(D2-(a^2/12)*D2^2)*x_old + (a*0.5)*(z1.*x_old)'*x_old*(1+a^2/12);
    dS = S-S_old;

    end    
    
end

function G = compute_G(x,n,N)
    % Function that calculates the nth order correlation function for a lattice x
    % It returns the correlation functtion

    % x is the lattice
    % n is the order of the correlation function
    % N is the lattice size

    A = eye(N);
    B = circshift(A, [0, n]); 
    G = x' * B * x/ N;

end
  
