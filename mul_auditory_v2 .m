
function mul_auditory_v2

% define param
N_E = 100;
N_I = 100;
tau_E = 10^-3;
tau_I = 10^-3;
tau_ref_E = 3*10^-3;
tau_ref_I = 3*10^-3;
tau_rec = 0.8;
U = 0.5;
J_EE = [6,4.5 * 10^-2,1.5 * 10^-2];
J_EI = -4;
J_IE = [0.5,3.5 * 10^-3,1.5 * 10^-3];
J_II = -0.5;
P = 15;

delta_left = 5;
delta_right = 5;
lamda_c = 0.25;
alpha = 2;

switch_fig = [4];


%  [E,I,x,y]
v0 =[zeros(N_E,1,P);zeros(N_I,1,P);zeros(N_E,1,P);zeros(N_I,1,P)]; 

% define background inputs,  uniform distribution,
e_E = sort(unifrnd(-10,10,N_E,P));
e_E = repmat(unifrnd(-10,10,N_E,1),1,P);
e_I = sort(unifrnd(-10,10,N_I,P));
e_E = repmat(unifrnd(-10,10,N_I,1),1,P);


if (find(switch_fig == 1))

A = zeros(1,P);
A(8)=4;
s = sensory_input(A);
[mOE,mOI,mOx,mOy,tt,~] = ode_loop(s);

clf;
figure(1);
subplot(2,2,1); 
plot(tt,mOE);
xlabel('t'); ylabel('E')
subplot(2,2,2); 
plot(tt,mOI);
xlabel('t'); ylabel('I')
subplot(2,2,3); 
plot(tt,mOx);
xlabel('t'); ylabel('x')
subplot(2,2,4); 
plot(tt,mOy);
xlabel('t'); ylabel('y')
end

%Response of the model to a sustained pure tone input
if (find(switch_fig ==2))
    A = zeros(1,P);
    A(8)=4;
    s = sensory_input(A);
    [mOE,~,~,~,tt,~] = ode_loop(s);
    clf;
    figure(2);
    for col = 6:10 
        subplot(3,5,col-5);  
        plot(tt,mOE(:,col));
        ylim([0,100]);

        subplot(3,5,col); 
        plot(tt,s(:,col));
    end
subplot(3,5,[11,12,13,14,15]);
plot(6:10,s(6:10),'k');
end

%Localization of PS induced activity.
if (find(switch_fig == 3))
    clf;
    figure(3);
    for s_value = 7:-1:1 
        A = zeros(1,P);
        A(8)=s_value;
        s = sensory_input(A);
        [mOE,~,~,~,tt,~] = ode_loop(s);
        subplot(8,1,-s_value+8);
        imagesc(1:P,tt,mOE);
        %imagesc(mOE);
        set(gca,'ydir','normal')
        ylim([0 0.05])
        colormap;
        
        if (s_value == 2 || s_value == 4 || s_value == 7)
            subplot(8,1,8);
            plot(1:P,s);
            hold on;
        end
    end
    

end

if (find(switch_fig == 4))
    clf;
    figure(4)
    
    %Temporal co-tuning of excitation and inhibition
    %1.weak pure tone inputs
    for s_value = [2 8]
    A = zeros(1,P);
    A(8)=s_value;
    s = sensory_input(A);
    [mOE,mOI,~,~,tt,~] = ode_loop(s);
    tt = tt*1000;
    subplot(3,1,1*(s_value ==2)+2*(s_value ==8));
    plot(tt,mOE(:,8),':k');
    hold on;
    plot(tt,mOI(:,8),'k');
    xlim([-10,100]);
    end

    for s_value = 2:0.2:10
        A(8)=s_value;
        s = sensory_input(A);
        %[mOE,mOI,~,~,tt,~] = ode_loop(s);
    
    end


end




function [mOE,mOI,mOx,mOy,tt,xx] = ode_loop(s)

s0 = zeros(1,P);
[tt,xx] = ode45(@ACode,[-5:0.001:0],v0,[],s0);

tspan = [0:0.001:1];


v0 = xx(end,:);
[tt2,xx2] = ode45(@ACode,tspan,v0,[],s);

xx = [xx; xx2];
tt = [tt;tt2];

% mean each column
OE = xx(:,1:N_E*P);
rOE = reshape(OE,[],N_E,P);
mOE = squeeze(mean(rOE,2));
OI = xx(:,N_E*P+1:N_E*P+N_I*P);
rOI = reshape(OI,[],N_I,P);
mOI = squeeze(mean(rOI,2));
Ox = xx(:,N_E*P+N_I*P+1: 2*N_E*P+N_I*P);
rOx = reshape(Ox,[],N_E,P);
mOx = squeeze(mean(rOx,2));
Oy = xx(:,2*N_E*P+N_I*P+1:end);
rOy = reshape(Oy,[],N_I,P);
mOy = squeeze(mean(rOy,2));

end




% define ode function
function ACode = ACode(t,v,s)
% WHAT is v, what is its shape?

% what is E? How does this extract E from v?
    E = reshape(v(1:N_E*P),N_E,P);
    I = reshape(v(N_E*P+1:N_E*P+N_I*P),N_I,P);
    x = reshape(v(N_E*P+N_I*P+1: 2*N_E*P+N_I*P),N_E,P);
    y = reshape(v(2*N_E*P+N_I*P+1:end),N_I,P);
    
    % single column,s = 0
    % s represents sensory input
%    s = 0; % no input at all
    % lets model a sudden change from no input to something...
    % another way
%     s = 10*(t>=1);
    % problems with this:
    % % (1) the times and magnitude are hard coded, whjat if
    % you want to change them? 
    % % (2)Numerical methods sometimgs have problems with sudden changes
    % liek this
    % （3） alternative hthis:
    % run for some time to get to unsrtimulated staeady state
    % then change the inputs and run again
    
    % gain function [z] = max(z,0)
    
    % As to E
    E_gain = zeros(N_E,P);
    I_gain = zeros(N_I,P);
    for Q = 1:P 
        if (Q == 1) 
            R_new = 0:2;
        elseif (Q == 2)
            R_new = -1:2;
        elseif (Q == P)
            R_new = -2:0;  
        elseif (Q == P-1)
            R_new = -2:1;  
        else 
            R_new = -2:2;  
        end
        
        E_gain_left = 0;
        I_gain_left = 0;
        for R = R_new
            E_sum1 = (J_EE(abs(R)+1)/N_E)*sum(U*x(:,Q+R).*E(:,Q+R));
            E_gain_left = E_gain_left + E_sum1;

            I_sum1 = (J_IE(abs(R)+1)/N_E)*sum(E(:,Q+R));
            I_gain_left = I_gain_left + I_sum1;
        end

        % define sensory input function
        % each tone is presented with an amplitude
        % tones can be at any of 15 columns (frequencies corresponding to
        % those columns)
        % each tone has a temporal envelope - typically it is just off and
        % then on
        % each tone spreads across columns from the "best frequency" with
        % exponential decay
        % So every column gets a sum of inputs from any tones, each reduced
        % by the exponential decay
        % Imagine two tones, one at column 4 and one at column 8
        % This is A = [0 0 0 4 0 0 0 8 0 0 0 0 0 0 0]
        % these spread, the higher amplitude one spreads further and we get
        % an overall set of inputs something like
        % s = [0.1 1 2 4.1 2.5 5 6.5 8 6 3 1 0.1 0 0 0]

        E_gain(:,Q) = E_gain_left + J_EI/N_I*sum(U*y(:,Q).*I(:,Q))+e_E(:,Q) + s(:,Q);
        I_gain(:,Q) = I_gain_left + J_II/N_I*sum(I(:,Q))+e_I(:,Q);

    end


    ACode = [reshape((-E + (1-tau_ref_E.*E).*max(E_gain,0))/tau_E,N_E*P,1); ...
             reshape((-I + (1-tau_ref_I.*I).*max(I_gain,0))/tau_I,N_E*P,1); ...
             reshape((1-x)/tau_rec-(U.*x.*E),N_E*P,1); ...
             reshape((1-y)/tau_rec-(U.*y.*I),N_E*P,1);
             ];

end

% define sensory input function,for loop 
function s1 = sensory_input(A)
      s1 = zeros(1,P);
      for Q = 1:P
          s = 0;
          for M = 1:P 
            if (Q<M)
                delta = delta_left;
            elseif (Q>=M)
                delta = delta_right;
            end
            
            s = s + A(M)*exp(-abs(Q-M)/(lamda_c + heaviside(A(M)-alpha)*(A(M)-alpha)/delta));
            
          end
      s1(Q) = s;
      end
    
end

end

















