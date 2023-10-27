
function mul_auditory_v3

rng(10);
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

switch_fig = [7];


%  [E,I,x,y]
v0 =[zeros(N_E,1,P);zeros(N_I,1,P);zeros(N_E,1,P);zeros(N_I,1,P)]; 

% define background inputs,  uniform distribution,
e_E = sort(unifrnd(-10,10,N_E,P));
e_E = repmat(sort(unifrnd(-10,10,N_E,1)),1,P);
e_I = sort(unifrnd(-10,10,N_I,P));
e_I = repmat(sort(unifrnd(-10,10,N_I,1)),1,P);

 e_on = repmat([zeros(1,40),ones(1,60)],1,P);
% i_on = repmat([ones(1,41),zeros(1,59)]',1,P);
% e_E = e_on.*e_E;
% e_I = i_on.*e_I;



% It takes about three minutes to run the ftc function once, 
% and about 30 minutes to run the entire Figure 5
if (find(switch_fig == 5))
    %clf(figure(5));
    figure(5);
    subplot(2,2,1);
    ftc = frequency_tuning_curve;
    plot(1:P,ftc,'k');
    hold on;
    %scatter(1:P,ftc,'filled');
    fprintf('2.%s\n',datetime('now','TimeZone','local','Format','y-MM-dd HH:mm:ss:SSS'));
    J_EE = [6,0,0];
    J_IE = [0.5,0,0];
    ftc = frequency_tuning_curve;
    plot(1:P,ftc,'k:+');
    legend('with inter-column connections', 'without inter-column connections','Location','SouthEast');
    ylim([0,9]);
    xlim([1,15]);
    ylabel('Input amplitude [Hz]');
% %     
%     % reset 
    J_EE = [6,4.5 * 10^-2,1.5 * 10^-2];
    J_IE = [0.5,3.5 * 10^-3,1.5 * 10^-3];
%     
    subplot(2,2,3);
    line_style = ['k','r','g','b'];
    delta_v = [5 8 16 32];
    for d = 1:length(delta_v)
        delta_left = delta_v(d);
        ftc = frequency_tuning_curve;
        plot(ftc,line_style(d));
        hold on;
    end
    legend('d_l=5', 'd_l=8', 'd_l=16', 'd_l=32','Location','SouthEast');
    ylim([0,9]);
    xlim([1,15]);
    xlabel('Column index');
    ylabel('Input amplitude [Hz]');
    
   % 
   delta_left = 5;
   subplot(2,2,2);
   J_II = J_II*0.1;
   J_EI = J_EI*0.1;
   ftc = frequency_tuning_curve;
   plot(ftc,'k');
   J_II = J_II*2;
   J_EI = J_EI*2;
   hold on;
   ftc = frequency_tuning_curve;
   plot(ftc,'k:');
   legend('decreased recurrent inh', 'increased recurrent inh','Location','SouthEast');
   ylim([0,9]);
   xlim([1,15]);

   
   J_EI = -4;
   J_II = -0.5;
   subplot(2,2,4);
   e_E = repmat(sort(unifrnd(-10,10,N_E,1)),1,P)+2;
   e_I = repmat(sort(unifrnd(-10,10,N_I,1)),1,P)+2;
   ftc = frequency_tuning_curve;
   plot(ftc,'k');
   hold on;
   e_E = repmat(sort(unifrnd(-10,10,N_E,1)),1,P)-2;
   e_I = repmat(sort(unifrnd(-10,10,N_I,1)),1,P)-2;
   ftc = frequency_tuning_curve;
   plot(ftc,'k:');
   legend('increased background input', 'decreased background input','Location','SouthEast');
   ylim([0,9]);
   xlim([1,15]);
   xlabel('Column index');


end

% Forward masking: dynamics of recovery
% about 4 minutes to run 
if (find(switch_fig == 6))
    
     [tt0,xx0,v00] = ode_loop_first;


    clf(figure(6));
    figure(6);
    subplot(3,1,1);
    ISI_v = [0.125 0.25 0.5 1 2 4];

    for isi = 1:length(ISI_v)
    
   % ts = [0 0.1 ISI_v(isi) 5];
    ts = [0 0.05 ISI_v(isi)+0.05 ISI_v(isi)+0.1 5];
    A = zeros(length(ts)-1,P);
    A(1,8)=10; 
    A(2,8)=0;
    A(3,8)=10;

    [mOE,~,~,~,tt,~] = ode_loop(A,ts,xx0,tt0,v00);
    plot(tt,mOE(:,8),'k');
    xlim([-0.1,5]);
    
    hold on;
    end
    ylabel('Network Activity[Hz]');

    % B
    subplot(3,1,2);
    P2P1 = zeros(3,5)*nan;
    input_A = [4 5.5 10];
    line_style = ["k:+" "k--*" "k-o"];
    for s_value = 1:length(input_A)
        ISI_v = [0.125 0.25 0.5 1 2 4];
        [tt0,xx0,v00] = ode_loop_first;
        for isi = 1:length(ISI_v)
        
        ts = [0 0.05 ISI_v(isi)+0.05 ISI_v(isi)+0.1 5];
        A = zeros(length(ts)-1,P);
        A(1,8)=input_A(s_value); 
        A(3,8)=input_A(s_value);
    
        [mOE,~,~,~,tt,~] = ode_loop(A,ts,xx0,tt0,v00);

        [~,locs]=findpeaks(mOE(:,8),'MinPeakHeight', 5);
        if size(locs,1)>=2
           P2P1(s_value,isi) = mOE(locs(2),8)/mOE(locs(1),8);
           continue;
        end

        end
        plot(ISI_v,P2P1(s_value,:),line_style(s_value),'Linewidth', 1,'MarkerSize', 3);
        hold on;
    end
    legend('input Amplitude 4', 'input Amplitude 5.5','input Amplitude 10','Location','SouthEast');
    ylim([0.1 1]);
    ylabel('P2/P1 Ratio');

      % C
    subplot(3,1,3);
    P2P1 = zeros(5,5)*nan;
    column = [3 4 5 6 8];
    line_style = ["c:","g--","b-.","r:","k-"];
    for col = 1:length(column)
        ISI_v = [0.125 0.25 0.5 1 2 4];
        [tt0,xx0,v00] = ode_loop_first;
        for isi = 1:length(ISI_v)
        
        ts = [0 0.05 ISI_v(isi)+0.05 ISI_v(isi)+0.1 5];
        A = zeros(length(ts)-1,P);
        A(1,column(col))=10; 
        A(3,column(col))=10;
    
        [mOE,~,~,~,tt,~] = ode_loop(A,ts,xx0,tt0,v00);

        [~,locs]=findpeaks(mOE(:,8),'MinPeakHeight', 5);
        if size(locs,1)>=2
           P2P1(col,isi) = mOE(locs(2),8)/mOE(locs(1),8);
           continue;
        end

        end
        plot(ISI_v,P2P1(col,:),line_style(col),'Linewidth', 1);
        hold on;
    end
    legend('Column 3', 'Column 4','Column 5', 'Column 6','Column 8','Location','SouthEast');
    ylim([0.1 1]);
    ylabel('P2/P1 Ratio');
    xlabel('Inter-stimulus Interval');

end

if (find(switch_fig == 7))
    clf(figure(7));
    figure(7);
    ISI_v = [1/8 1/4 1/2 1 2 4]*tau_rec;
    ISI_v_str = ["1/8*{\tau_rec}" "1/4*{\tau_rec}" "1/2*\tau_rec" "1*\tau_rec" "2*\tau_rec" "4*\tau_rec"];
    y_start = [1.5 1.6 1.8 2.0 2.4 3.2];
    y_end = [1.75 2.1 2.3 2.6 3.1 3.9];

    for isi = 1:length(ISI_v)
    
    
    ts = [0 0.05 ISI_v(isi)+0.05 ISI_v(isi)+0.1 5];
    A = zeros(length(ts)-1,P);
    A(1,4)=10; 
    A(3,4)=10;

    [tt0,xx0,v00] = ode_loop_first;
    [mOE,~,~,~,tt,~] = ode_loop(A,ts,xx0,tt0,v00);
    subplot(7,1,7-isi);
    imagesc(1:P,tt,mOE);
    set(gca,'ydir','normal')
    title([ISI_v_str(isi)],'FontSize',6);
    %plot(tt,mOE);
%     subplot(7,2,isi+6);
%     plot(tt,mOE(:,4:8));
    %ylim([y_start(isi) y_end(isi)]);
    ylim([ISI_v(isi)+0.05 ISI_v(isi)+0.1]);
    colorbar;

    end
    subplot(7,1,7);
    imagesc(1:P,tt,mOE);
    set(gca,'ydir','normal')
    ylim([0 0.05]);

end


function ftc = frequency_tuning_curve
    warning('off','all');

    ts = [0 0.05];
    ftc = zeros(P,1);
    [tt0,xx0,v00] = ode_loop_first;
   
    for s_value1 = 1:0.2:10
        A = zeros(length(ts)-1,P);
        A(1,8)=s_value1;
       
        [mOE,~,~,~,tt,~] = ode_loop(A,ts,xx0,tt0,v00);
      
        [~,locs]=findpeaks(mOE(:,8),tt,MinPeakProminence=10);
        if size(locs)>=1
           ftc(8) = s_value1;
           break
        else
           ftc(8) = nan;
        end
    end
   
    
    a = [-1 1];
    b = [1 P];
    
    for i = 1:2 
    start = ftc(8);
    for col1 = 8+a(i):a(i):b(i)
        
        for s_value1 = start:0.2:10

            A = zeros(length(ts)-1,P);
            A(1,col1)=s_value1;
           
            [mOE,~,~,~,tt,~] = ode_loop(A,ts,xx0,tt0,v00);
            
           
            [~,locs]=findpeaks(mOE(:,8),tt,MinPeakProminence=10);
            %fprintf('2.%s\n',datetime('now','TimeZone','local','Format','y-MM-dd HH:mm:ss:SSS'));
            if size(locs)>=1
               ftc(col1) = s_value1;
               start = s_value1;
               break
            else
               ftc(col1) = nan;
            end
        end
    end
    end 

end

function [tt0,xx0,v00] = ode_loop_first

% xx = v0(:)';
% tt = ts(1);

s = zeros(1,P);
[ttemp,xtemp] = ode45(@ACode,[-5:0.001:0],v0,[],s);
v00 = xtemp(end,:);
tt0 = ttemp;
xx0 = xtemp;
%     cur_date = date;
%     cur_time = fix(clock);
fprintf('first\n');
fprintf('1.%s\n',datetime('now','TimeZone','local','Format','y-MM-dd HH:mm:ss:SSS'));


end


function [mOE,mOI,mOx,mOy,tt,xx] = ode_loop(A,ts,xx0,tt0,v00)

% xx = v0(:)';
% tt = ts(1);
xx = [];
tt = [];
% xx = xx0;
% tt = tt0;
for i  = 1:length(ts)-1
     
    s = sensory_input(A(i,:));
    %[ttemp,xtemp] = ode45(@ACode,[ts(i):0.001:ts(i+1)],v00,[],s);
    %fprintf('1.%s\n',datetime('now','TimeZone','local','Format','y-MM-dd HH:mm:ss:SSS'));
    [ttemp,xtemp] = ode45(@ACode,[ts(i):0.001:ts(i+1)],v00,[],s);
    
    v00 = xtemp(end,:);
    xx = [xx;xtemp(2:end,:)];
    tt = [tt;ttemp(2:end)];
%     xx = [xx;xtemp];
%     tt = [tt;ttemp];

    

end

 
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

    e_on = reshape(e_on,N_E,P);
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

        E_gain(:,Q) = E_gain_left + J_EI/N_I*sum(U*y(:,Q).*I(:,Q))+e_E(:,Q) + s(:,Q).*e_on(:,Q);
        I_gain(:,Q) = I_gain_left + J_II/N_I*sum(I(:,Q))+e_I(:,Q);

    end


    ACode = [reshape((-E + (1-tau_ref_E.*E).*max(E_gain,0))/tau_E,N_E*P,1); ...
             reshape((-I + (1-tau_ref_I.*I).*max(I_gain,0))/tau_I,N_I*P,1); ...
             reshape((1-x)/tau_rec-(U.*x.*E),N_E*P,1); ...
             reshape((1-y)/tau_rec-(U.*y.*I),N_I*P,1);
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

















