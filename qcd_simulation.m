% A N^4 quantum lattice is created and link operators built on that.
% The link operators are made evolve through a Metropolis algoritm based on
% 1x1 wilson loop action, with the addition of 2x1 wilson loops for the
% improved action. The actionis purely gluonic
% The average sive of 1x1 and 2x1 wilson loops are measured through the
% simulation and printed at the end of it.

clear

[tot_avg1_link,err_avg1_link,tot_avg2_link,err_avg2_link,acceptance_ratio,elapsed_time]=main;

fprintf(['Average plaquette size is ',num2str(tot_avg1_link),' ± ',num2str(err_avg1_link),'.\n',...
         'Average rectangular size is ',num2str(tot_avg2_link),' ± ',num2str(err_avg2_link),'.\n',...
         'Total elapsed time is ',num2str(elapsed_time),' seconds. \n',...
         'The acceptance ratio for the metropolis algoritm is ',num2str(acceptance_ratio),' (should be 0.5 ± 0.1). \n'])

function [tot_avg1_link,err_avg1_link,tot_avg2_link,err_avg2_link,acceptance_ratio,elapsed_time]=main(setup)
    tic;
    entries={      'N (lattice size)',...
                   'N_COR (Number of sweeps before measurement) ',...
                   'N_CF (Number of measurements)',...
                   'N_MA (Metropolis steps for each link)',...
                   'eps (Epsylon of randomness for the SU_3 matrices)',...
                   'NUM_MAT (Number of SU(3) matrices and their inverses used)',...
                   'Imp_act (false=simplified action, true=improved action)'};
    if(nargin==0)
        default={  '8',...
                   '50',...
                   '2',...
                   '10',...
                   '0.34',...
                   '100',...
                   'true'};     
    else
        default=setup;        
    end
  
    addopt.Resize='on';
    addopt.WindowStyle='normal';
    addopt.Interpreter = 'none';

    header = 'QCD simulation';
    lines = 1;
    setup = inputdlg(entries, header, lines, default, addopt); 
  
    N  = str2num(setup{1});         % Lattice size
    N_COR = str2num(setup{2});      % Number of sweeps before measurement
    term = N_COR*2;                 % Temralization sweeps before the start
    N_CF  = str2num(setup{3});      % Number of measurements
    N_MA = str2num(setup{4});       % Metropolis steps for each link
    eps  = str2num(setup{5});       % Epsylon of randomness for the SU_3 matrices
    NUM_MAT = str2num(setup{6});    % Number of SU(3) matrices and their inverses used
    imp_act  = str2num(setup{7});   % false=simplified action, true=improved action
    dim=3;                          % Dimension of the SU(3) matrices
    
    par=0;                          % Number of metropolis steps acccepted
    tot=0;                          % Total steps

    MU0=0.797;                      % Averege link size
    c_1=5/(3*MU0^4);                % Coefficient for plaquette
    c_2=1/(12*MU0^6);               % Coefficient for rectangle operator
    if imp_act
        BETA = 1.719;               % =6/g^2 coefficient for the action      
    else
        BETA = 5.5;                 % =6/g^2 coefficient for the action (include MU0)
    end

    plaq_size=zeros(N,N,N,N,6);     % Collector where plaquette size are stacked
    rect2_size=zeros(N,N,N,N,6);    % Collector where rectangular size are stacked
    avg_link1=zeros(1,N_CF);        % Collector of the averege linksize for the order 1 wilson loop
    avg_link2=zeros(1,N_CF);        % Collector of the averege linksize for the order 2 wilson loop

    hw = waitbar(0,'Running...');   

    SU3_LIST=create_SU3_list(NUM_MAT,dim,eps);

    L=create_lattice(N,dim);

    % Sweeps of termalization before to start measuring
    for i=1:term
        [L,par,tot]=sweep(L,N_MA,N,SU3_LIST,BETA,par,tot,NUM_MAT,dim,c_1,c_2,imp_act);
    end
    % Simulation and measurment
    for i=1:N_CF
    
        for j=1:N_COR
            [L,par,tot]=sweep(L,N_MA,N,SU3_LIST,BETA,par,tot,NUM_MAT,dim,c_1,c_2,imp_act);
        end

        plaq_size=measure_plaq(L,plaq_size,N);      % Calculate all the plaquette
        avg_link1(1,i) = mean(plaq_size,'all');     % Calculate the averege size of the plaquettes

        rect2_size=measure_rect2(L,rect2_size,N);   % Calculate all the 2x1 vertical wilson loops
        avg_link2(1,i) = mean(rect2_size,'all');    % Calculate the averge size of the 2x1 loops

        waitbar(i/N_CF);
    end

    % Mean and standard deviation of the plaquettes
    tot_avg1_link=mean(avg_link1,'all');    
    err_avg1_link=std(avg_link1);
    % Mean and standard deviation of the 2x1 loops
    tot_avg2_link=mean(avg_link2,'all');
    err_avg2_link=std(avg_link2);
    % Acceptance ratio for the metropolis
    acceptance_ratio = par/tot;

    close(hw);
    elapsed_time=toc;
end

function SU3_LIST=create_SU3_list(NUM_MAT,dim,eps_herm)
 % Function that creates a list of random SU(3) matrices (each inverse
 % matrix is also present
 % It returns the list

 % NUM_MAT is the size of the list (halved since we have also inverse)
 % dim is the dimension of the matrices
 % eps_herm is the epsylon of randomness for the SU_3 matrices

    SU3_LIST(1:NUM_MAT*2)={eye(dim)};

    for i=1:NUM_MAT
        SU3_LIST{i}=generate_suN(dim,eps_herm);       
        SU3_LIST{NUM_MAT*2-i+1}=SU3_LIST{i}';     % Add also the inverse        
    end

end

function [L,par,tot]=sweep(L,N_ma,N,SU3_LIST,BETA,par,tot,NUM_MAT,dim,c_1,c_2,imp_act)
 % Function that updates the lattice
 % This routine does N_ma metropolis steps for each lattice sites and link (N_ma*N^4*4 in total)
 % It gives back the updated lattice and the acceptance ratio par/tot
 
 % L contains all the links of the lattice
 % N_ma steps of thermalization for each link
 % N lattice size
 % BETA,c_1 and c_2 are coefficients for the action
 % par and tot needed for acceptance ratio
 % SU3_LIST,NUM_MAT and dim are needed for the random SU(3) matrices
 staple_rect=0;
 staple_rect_conj=0;
    for x=1:N
        for y=1:N
            for z=1:N
                for t=1:N
                    for n=1:4
                        [staple,staple_conj]=plaq(L,x,y,z,t,n,N);
                        if imp_act
                            [staple_rect,staple_rect_conj]=rect(L,x,y,z,t,n,N);
                        end

                        for i=1:N_ma
                            M=SU3_LIST{randi(NUM_MAT*2)};
                            delta_S = dS(L,x,y,z,t,n,staple,staple_conj,staple_rect,staple_rect_conj,M,BETA,dim,c_1,c_2,imp_act);
                            
                            if delta_S < 0 || rand < exp(-delta_S)
                                L(:,:,x,y,z,t,n) = M*L(:,:,x,y,z,t,n);
                                par=par+1;
                            end
                            tot=tot+1;
                        end
                    end
                end
            end
        end
    end

end

function L=create_lattice(N,dim)
  % Function that create a lattice of links initializated to identities
  % First 2 entries of L are the one for the 4*4 matrix that describe the link, 
  % 3,4,5,6 are the coordinates for the N^4 lattice (3=x, 4=y, 5=z, 6=t),  
  % last entry (7) is the direction of the link (x,y,z,t) = (1,2,3,4)
  % For example L(:,:,1,4,3,3,2) is the link in the point (1,4,3,3) in direction y (y=2)
  % It returns the created lattice

  % N is the size of the lattice
  % dim is the size of the matrices

L(1:dim,1:dim,1:N,1:N,1:N,1:N,1:4)=0;

    for i=1:N
        for j=1:N
            for k=1:N
                for l=1:N
                    for m=1:4
                        L(:,:,i,j,k,l,m)=eye(dim);
                    end
                end
            end
        end
    end

end

function SU_N=generate_suN(N, eps_herm)
    % Function that creates a random SU(3) matrix
    % N is the size of the matrix
    % eps_herm is the range of random values of the hermitian matrix
    % It returns the matrix

    H = generate_hermitian(N, eps_herm);    % Hermitian
    U = expm(1i * H);                       % Unitary
    SU_N = U / (det(U)^(1/N));              % Normalized
end

function H=generate_hermitian(N, eps_herm)
    % Function that generates a random NxN hermitian matrix
    % N is size of the matrix
    % eps_herm is the range of random values of the hermitian matrix
    % It returns the matrix
    
    real_part = (rand(N) - 0.5) * 2 * eps_herm; 
    imag_part = (rand(N) - 0.5) * 2 * eps_herm;
    H = real_part + 1i * imag_part;  

    H = (H + H') / 2;
end

function [staple,staple_conj]=plaq(L,x,y,z,t,n,N)
    % Routine that calculates all the plaquettes involved in the change of
    % action for the link L(x,y,z,t,n).
    % The link itself is not included in the routine, so, in order to get
    % the full plaquette, one need to multiply by the link.
    % Also the contribution for the conjugate link are calculated and
    % stored in staple_conj
    % Anti-clockwise direction
    % It returns the two staples
    
    % L contains all the links of the lattice
    % x,y,z,t are the coordinates of the lattice while n is the direction
    % N size of the lattice

    M=10;
    a=zeros(1,M);
    a(n)=1;
    
    staple=0;
    staple_conj=0;

    for i=1:(4-n)
        % Staple construction for (n=1) xy,xz,xt
        %                         (n=2) yz,yt
        %                         (n=3) zt
        staple = staple + ...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)),m(z,N,a(3)),t,n+i)*...
        L(:,:,x,m(y,N,a(m(1,M,-i+1))),m(z,N,a(m(2,M,-i+1))),m(t,N,a(m(3,M,-i+1))),n)'*...
        L(:,:,x,y,z,t,n+i)';

        staple_conj= staple_conj +...  
        L(:,:,x,m(y,N,-a(m(1,M,-i+1))),m(z,N,-a(m(2,M,-i+1))),m(t,N,-a(m(3,M,-i+1))),n+i)'*...
        L(:,:,x,m(y,N,-a(m(1,M,-i+1))),m(z,N,-a(m(2,M,-i+1))),m(t,N,-a(m(3,M,-i+1))),n)*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)-a(m(1,M,-i+1))),m(z,N,a(3)-a(m(2,M,-i+1))),m(t,N,-a(m(3,M,-i+1))),n+i);

    end

    for i=1:(n-1)
        % Staple construction for (n=2) xy
        %                         (n=3) xz,yz
        %                         (n=4) xt,yt,zt
        staple = staple + ...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,a(2)-a(i+n-2)),m(z,N,a(3)-a(m(i,M,n-3))),m(t,N,a(4)),i)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(i+n-2)),m(z,N,-a(m(i,M,n-3))),t,n)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(i+n-2)),m(z,N,-a(m(i,M,n-3))),t,i);

        staple_conj= staple_conj +...
        L(:,:,x,y,z,t,i)*...
        L(:,:,m(x,N,+a(i+n-1)),m(y,N,+a(i+n-2)),m(z,N,+a(m(i,M,n-3))),t,n)*...
        L(:,:,x,m(y,N,a(2)),m(z,N,a(3)),m(t,N,a(4)),i)';

    end


end

function z=m(x,N,change)
 % Function that mimics the boundary condition of the lattice
 % It's similar to z=mod(x,N) but start at 1 and not 0
 % Only integer numbers can be used
 % It returns the updated number
 
 % x the variable that has to be changed
 % change the added number (has to be integer)
 % N the boundary
 
    if     change==0
        z=x;

    elseif change==1
        if x==N
            z=1;
        else
            z=x+1;
        end

    elseif change==-1
        if x==1
            z=N;
        else
            z=x-1;
        end
    
    elseif change > 1
        z=m(m(x,N,1),N,change-1);

    elseif change < -1
        z=m(m(x,N,-1),N,change+1);

    end

end

function delta_S=dS(L,x,y,z,t,n,staple,staple_conj,staple_rect,staple_rect_conj,M,BETA,dim,c_1,c_2,imp_act)
    % Routine that calculates the the difference in action given by a link
    % L(x,y,z,t,n), his staple and staple_conj and the update SU(3) matrix M
    % It returns the delta in action

    % L contains all the links of the lattice
    % x,y,z,t are the coordinates of the lattice while n is the direction
    % staple and staple_conj are the quantities that multiplied by the link
    % and his inverse respectively gives the plaquette
    % staple_rect and staple_rect_conj are just like the aboves but for 2x1 
    % loop
    % M is the SU(3) randim matrix
    % BETA,c_1 and c_2 are coefficients for the action
    % dim is the dimension of the matrices
    
    U = L(:,:,x,y,z,t,n);
    if imp_act
        delta_S = (BETA/3)*real(trace( (eye(dim)-M)*U* (c_1*staple-c_2*staple_rect) +...
                                    U'*(eye(dim)-M')* (c_1*staple_conj-c_2*staple_rect_conj) ) );
    else    
        delta_S = (BETA/3)*real(trace( (eye(dim)-M)*U*staple + U'*(eye(dim)-M')*staple_conj ) );
    end
end

function plaq_size=measure_plaq(L,plaq_size,N)
 % Function that calculates the wilson loop for the plaquette operator for
 % all the lattice points and all the direction and mean it.
 % Wilson loop is the product of the links in a loop.
 % It returns the measured loop

 % L contains all the links of the lattice
 % N size of the lattice
 % plaq_size is a (N,N,N,N,6) vector that contains all the plaquettes on a
 % lattice

    for x=1:N
        for y=1:N
            for z=1:N
                for t=1:N
                    
                    plaq_size(x,y,z,t,1) = real( trace( L(:,:,x,y,z,t,1)*L(:,:,m(x,N,1),y,z,t,2)*...
                                                      L(:,:,x,m(y,N,1),z,t,1)'*L(:,:,x,y,z,t,2)' ) );
                    plaq_size(x,y,z,t,2) = real( trace( L(:,:,x,y,z,t,1)*L(:,:,m(x,N,1),y,z,t,3)*...
                                                      L(:,:,x,y,m(z,N,1),t,1)'*L(:,:,x,y,z,t,3)' ) );
                    plaq_size(x,y,z,t,3) = real( trace( L(:,:,x,y,z,t,1)*L(:,:,m(x,N,1),y,z,t,4)*...
                                                      L(:,:,x,y,z,m(t,N,1),1)'*L(:,:,x,y,z,t,4)' ) );
                    plaq_size(x,y,z,t,4) = real( trace( L(:,:,x,y,z,t,2)*L(:,:,x,m(y,N,1),z,t,3)*...
                                                      L(:,:,x,y,m(z,N,1),t,2)'*L(:,:,x,y,z,t,3)' ) );
                    plaq_size(x,y,z,t,5) = real( trace( L(:,:,x,y,z,t,2)*L(:,:,x,m(y,N,1),z,t,4)*...
                                                      L(:,:,x,y,z,m(t,N,1),2)'*L(:,:,x,y,z,t,4)' ) );
                    plaq_size(x,y,z,t,6) = real( trace( L(:,:,x,y,z,t,3)*L(:,:,x,y,m(z,N,1),t,4)*...
                                                      L(:,:,x,y,z,m(t,N,1),3)'*L(:,:,x,y,z,t,4)' ) );

                end
            end
        end
    end

    plaq_size=plaq_size/3; % We divide it by 3 for the definition of wilson Loop

end

function rect2_size=measure_rect2(L,rect2_size,N)
 % Function that calculates the wilson loop for a rectangular vertical size
 % 2x1 operator for all the lattice point and all the direction and mean it.
 % Wilson loop is the product of the links in a loop.
 % It returns the measured loop

 % L contains all the links of the lattice
 % N size of the lattice
 % rect2_size is a (N,N,N,N,6) vector that contains all the 2x1 loops on a
 % lattice

    for x=1:N
        for y=1:N
            for z=1:N
                for t=1:N
                    
                    rect2_size(x,y,z,t,1) = real( trace( ...
                        L(:,:,x,y,z,t,1)                *   L(:,:,m(x,N,1),y,z,t,2)*...
                        L(:,:,m(x,N,1),m(y,N,1),z,t,2)  *   L(:,:,x,m(y,N,2),z,t,1)'*...
                        L(:,:,x,m(y,N,1),z,t,2)'        *   L(:,:,x,y,z,t,2)' ) );                                                 
                    rect2_size(x,y,z,t,2) = real( trace( ...
                        L(:,:,x,y,z,t,1)                *   L(:,:,m(x,N,1),y,z,t,3)*...
                        L(:,:,m(x,N,1),y,m(z,N,1),t,3)  *   L(:,:,x,y,m(z,N,2),t,1)'*...
                        L(:,:,x,y,m(z,N,1),t,3)'        *   L(:,:,x,y,z,t,3)' ) );
                    rect2_size(x,y,z,t,3) = real( trace( ...
                        L(:,:,x,y,z,t,1)                *   L(:,:,m(x,N,1),y,z,t,4)*...
                        L(:,:,m(x,N,1),y,z,m(t,N,1),4)  *   L(:,:,x,y,z,m(t,N,2),1)'*...
                        L(:,:,x,y,z,m(t,N,1),4)'        *   L(:,:,x,y,z,t,4)' ) );
                    rect2_size(x,y,z,t,4) = real( trace( ...
                        L(:,:,x,y,z,t,2)                *   L(:,:,x,m(y,N,1),z,t,3)*...
                        L(:,:,x,m(y,N,1),m(z,N,1),t,3)  *   L(:,:,x,y,m(z,N,2),t,2)'*...
                        L(:,:,x,y,m(z,N,1),t,3)'        *   L(:,:,x,y,z,t,3)' ) );
                    rect2_size(x,y,z,t,5) = real( trace( ...
                        L(:,:,x,y,z,t,2)                *   L(:,:,x,m(y,N,1),z,t,4)*...
                        L(:,:,x,m(y,N,1),z,m(t,N,1),4)  *   L(:,:,x,y,z,m(t,N,2),2)'*...
                        L(:,:,x,y,z,m(t,N,1),4)'        *   L(:,:,x,y,z,t,4)' ) );
                    rect2_size(x,y,z,t,6) = real( trace( ...
                        L(:,:,x,y,z,t,3)*L(:,:,x,y,m(z,N,1),t,4)*...
                        L(:,:,x,y,m(z,N,1),m(t,N,1),4)  *   L(:,:,x,y,z,m(t,N,2),3)'*...
                        L(:,:,x,y,z,m(t,N,1),4)'        *   L(:,:,x,y,z,t,4)' ) );

                end
            end
        end
    end

    rect2_size=rect2_size/3; % We divide it by 3 for the definition of wilson Loop

end

function [staple_rect,staple_rect_conj]=rect(L,x,y,z,t,n,N)
    % Routine that calculates all the 2x1 wilson loop involved in the change of
    % action for the link L(x,y,z,t,n)
    % The link itself is not included in the routine, so, in order to get
    % the full plaquette, one need to multiply by the link
    % Also the contribution for the conjugate link are calculated and
    % stored in staple_conj
    % Clockwise direction
    % It returns the staples

    % L contains all the links of the lattice
    % x,y,z,t are the coordinates of the lattice while n is the direction
    % N size of the lattice

    M=10;
    a=zeros(1,M);
    a(n)=1;
    
    staple_rect=0;
    staple_rect_conj=0;

    for i=1:(4-n)
        % Staple construction for (n=1) xy,xz,xt
        %                         (n=2) yz,yt
        %                         (n=3) zt
        staple_rect = staple_rect + ...
        ... % orizzontal rectangle oriented to the bottom-right (clockwise) (mu,nu)        
        L(:,:,m(x,N,a(1)),m(y,N,a(2)),m(z,N,a(3)),t,n)*...
        L(:,:,m(x,N,2*a(1)),m(y,N,2*a(2)-a(m(2,M,-i))),m(z,N,2*a(3)-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n+i)'*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)-a(m(2,M,-i))),m(z,N,a(3)-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n)'*...
        L(:,:,x,m(y,N,-a(m(2,M,-i))),m(z,N,-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n)'*...
        L(:,:,x,m(y,N,-a(m(2,M,-i))),m(z,N,-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n+i)+...
        ... % orizzontal rectangle oriented to the bottom-left (clockwise) (mu,nu)
        L(:,:,m(x,N,a(1)),m(y,N,a(2)-a(m(2,M,-i))),m(z,N,a(3)-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n+i)'*...
        L(:,:,x,m(y,N,-a(m(2,M,-i))),m(z,N,-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n)'*...
        L(:,:,m(x,N,-a(1)),m(y,N,-a(2)-a(m(2,M,-i))),m(z,N,-a(3)-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n)'*...
        L(:,:,m(x,N,-a(1)),m(y,N,-a(2)-a(m(2,M,-i))),m(z,N,-a(3)-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n+i)*...
        L(:,:,m(x,N,-a(1)),m(y,N,-a(2)),m(z,N,-a(3)),t,n)+...
        ... % vertical rectangle oriented to the top (anti-clockwise) (nu,mu)
        L(:,:,m(x,N,a(1)),m(y,N,a(2)),m(z,N,a(3)),t,n+i)*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)+a(m(2,M,-i))),m(z,N,a(3)+a(m(3,M,-i))),m(t,N,a(m(4,M,-i))),n+i)*...
        L(:,:,x,m(y,N,2*a(m(2,M,-i))),m(z,N,2*a(m(3,M,-i))),m(t,N,2*a(m(4,M,-i))),n)'*...
        L(:,:,x,m(y,N,a(m(2,M,-i))),m(z,N,a(m(3,M,-i))),m(t,N,a(m(4,M,-i))),n+i)'*...
        L(:,:,x,y,z,t,n+i)';

        staple_rect_conj= staple_rect_conj +...  
        ... % orizzontal rectangle oriented to the top-right (clockwise) (mu,nu)        
        L(:,:,x,y,z,t,n+i)*...
        L(:,:,x,m(y,N,+a(m(2,M,-i))),m(z,N,+a(m(3,M,-i))),m(t,N,+a(m(4,M,-i))),n)*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)+a(m(2,M,-i))),m(z,N,a(3)+a(m(3,M,-i))),m(t,N,+a(m(4,M,-i))),n)*...
        L(:,:,m(x,N,2*a(1)),m(y,N,2*a(2)),m(z,N,2*a(3)),t,n+i)'*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)),m(z,N,a(3)),t,n)'+...
        ... % orizzontal rectangle oriented to the top-left (clockwise) (mu,nu)
        L(:,:,m(x,N,-a(1)),m(y,N,-a(2)),m(z,N,-a(3)),t,n)'*...
        L(:,:,m(x,N,-a(1)),m(y,N,-a(2)),m(z,N,-a(3)),t,n+i)*...
        L(:,:,m(x,N,-a(1)),m(y,N,-a(2)+a(m(2,M,-i))),m(z,N,-a(3)+a(m(3,M,-i))),m(t,N,a(m(4,M,-i))),n)*...
        L(:,:,x,m(y,N,a(m(2,M,-i))),m(z,N,a(m(3,M,-i))),m(t,N,a(m(4,M,-i))),n)*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)),m(z,N,a(3)),t,n+i)'+...
        ... %  vertical rectangle oriented to bottom (anti-clockwise) (nu,mu)
        L(:,:,x,m(y,N,-a(m(2,M,-i))),m(z,N,-a(m(3,M,-i))),m(t,N,-a(m(4,M,-i))),n+i)'*...
        L(:,:,x,m(y,N,-2*a(m(2,M,-i))),m(z,N,-2*a(m(3,M,-i))),m(t,N,-2*a(m(4,M,-i))),n+i)'*...
        L(:,:,x,m(y,N,-2*a(m(2,M,-i))),m(z,N,-2*a(m(3,M,-i))),m(t,N,-2*a(m(4,M,-i))),n)*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)-2*a(m(2,M,-i))),m(z,N,a(3)-2*a(m(3,M,-i))),m(t,N,a(4)-2*a(m(4,M,-i))),n+i)*...
        L(:,:,m(x,N,a(1)),m(y,N,a(2)-a(m(2,M,-i))),m(z,N,a(3)-a(m(3,M,-i))),m(t,N,a(4)-a(m(4,M,-i))),n+i);
   
    end

    for i=1:(n-1)
        % Staple construction for (n=2) xy
        %                         (n=3) xz,yz
        %                         (n=4) xt,yt,zt
        staple_rect = staple_rect + ...
        ... % orizzontal rectangle oriented to the right (clockwise) (mu,nu)
        L(:,:,x,m(y,N,a(2)),m(z,N,a(3)),m(t,N,a(4)),i)*...
        L(:,:,m(x,N,a(i+n-1)),m(y,N,a(2)+a(i+n-2)),m(z,N,a(3)+a(m(i,M,n-3))),m(t,N,a(4)),i)*...
        L(:,:,m(x,N,2*a(i+n-1)),m(y,N,2*a(i+n-2)),m(z,N,2*a(m(i,M,n-3))),t,n)'*...
        L(:,:,m(x,N,a(i+n-1)),m(y,N,a(i+n-2)),m(z,N,a(m(i,M,n-3))),t,i)'*...
        L(:,:,x,y,z,t,i)'+...
        ... % vertical rectangle oriented to the top-left (anti-clockwise) (nu,mu)         
        L(:,:,x,m(y,N,a(2)),m(z,N,a(3)),m(t,N,a(4)),n)*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,2*a(2)-a(i+n-2)),m(z,N,2*a(3)-a(m(i,M,n-3))),m(t,N,2*a(4)),i)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,a(2)-a(i+n-2)),m(z,N,a(3)-a(m(i,M,n-3))),m(t,N,a(4)),n)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(i+n-2)),m(z,N,-a(m(i,M,n-3))),t,n)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(i+n-2)),m(z,N,-a(m(i,M,n-3))),t,i)+...
        ... % vertical rectangle oriented to the bottom-left (anti-clockwise) (nu,mu)        
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,a(2)-a(i+n-2)),m(z,N,a(3)-a(m(i,M,n-3))),m(t,N,a(4)),i)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(i+n-2)),m(z,N,-a(m(i,M,n-3))),t,n)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(2)-a(i+n-2)),m(z,N,-a(3)-a(m(i,M,n-3))),m(t,N,-a(4)),n)'*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(2)-a(i+n-2)),m(z,N,-a(3)-a(m(i,M,n-3))),m(t,N,-a(4)),i)*...
        L(:,:,x,m(y,N,-a(2)),m(z,N,-a(3)),m(t,N,-a(4)),n);

        staple_rect_conj= staple_rect_conj +...
        ... % orizzontal rectangle oriented to the left (clockwise) (mu,nu)
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,-a(i+n-2)),m(z,N,-a(m(i,M,n-3))),t,i)'*...
        L(:,:,m(x,N,-2*a(i+n-1)),m(y,N,-2*a(i+n-2)),m(z,N,-2*a(m(i,M,n-3))),t,i)'*...
        L(:,:,m(x,N,-2*a(i+n-1)),m(y,N,-2*a(i+n-2)),m(z,N,-2*a(m(i,M,n-3))),t,n)*...
        L(:,:,m(x,N,-2*a(i+n-1)),m(y,N,a(2)-2*a(i+n-2)),m(z,N,a(3)-2*a(m(i,M,n-3))),m(t,N,a(4)),i)*...
        L(:,:,m(x,N,-a(i+n-1)),m(y,N,a(2)-a(i+n-2)),m(z,N,a(3)-a(m(i,M,n-3))),m(t,N,a(4)),i)+...
        ... % vertical rectangle oriented to the top-right (anti-clockwise) (nu,mu)         
        L(:,:,x,y,z,t,i)*...
        L(:,:,m(x,N,a(i+n-1)),m(y,N,a(i+n-2)),m(z,N,a(m(i,M,n-3))),t,n)*...
        L(:,:,m(x,N,a(i+n-1)),m(y,N,a(2)+a(i+n-2)),m(z,N,a(3)+a(m(i,M,n-3))),m(t,N,a(4)),n)*...
        L(:,:,x,m(y,N,2*a(2)),m(z,N,2*a(3)),m(t,N,2*a(4)),i)'*...
        L(:,:,x,m(y,N,a(2)),m(z,N,a(3)),m(t,N,a(4)),n)'+...
        ... % vertical rectangle oriented to the bottom-right (anti-clockwise) (nu,mu)         
        L(:,:,x,m(y,N,-a(2)),m(z,N,-a(3)),m(t,N,-a(4)),n)'*...
        L(:,:,x,m(y,N,-a(2)),m(z,N,-a(3)),m(t,N,-a(4)),i)*...
        L(:,:,m(x,N,a(i+n-1)),m(y,N,-a(2)+a(i+n-2)),m(z,N,-a(3)+a(m(i,M,n-3))),m(t,N,-a(4)),n)*...
        L(:,:,m(x,N,a(i+n-1)),m(y,N,a(i+n-2)),m(z,N,a(m(i,M,n-3))),t,n)*...
        L(:,:,x,m(y,N,a(2)),m(z,N,a(3)),m(t,N,a(4)),i)';

    end

end








