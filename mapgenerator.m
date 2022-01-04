%% Instruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is To generate the id-iq map according to e-mot characters.
% update:
% JH 20211212 first edition release.
% JH 20211214 integrate the display function and improve the generating speed;
% JH 20211228 bugfix for low torque area
% powered by:
% jason.huang3@magna.com (MAGNA_PT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% configure and init
close all;clc;clear;
t=0;
tic;
draw_flg = 0;
map_command = 1;
% draw_flg = 1;
% map_command = 0;
motor_idx = 2;



%% constants of E-motor

global r np psi_f li lad laq Ld Lq imax umax cp slt_rat beta;

%
switch motor_idx
    case 1
        r = 0.004;      %phase reistance
        np = 4;             %pole number
        psi_f = 0.055 ;    %permanent magnet flux
        li = 20e-6;         %leakage inductance
        lad = 98e-6;      %Magnetizing inductance
        laq = 300e-6;
        Ld = li + lad;
        Lq = li + laq;
        imax = 550;
        umax = 136;
        cp = -psi_f/Ld;    %critical point
        slt_rat = Lq / Ld;        %Salient rate
        beta = slt_rat - 1;
        
    case 2
        
        r = 0.04;      %phase reistance
        np = 3;             %pole number
        psi_f = 0.178 ;    %permanent magnet flux
        
        Ld = 1e-3;
        Lq = 1.7e-3;
        
        imax = 250;
        umax = 350;
        cp = -psi_f/Ld;    %critical point
        slt_rat = Lq / Ld;        %Salient rate
        beta = slt_rat - 1;
end


%% E-motr operation points

id = imax * (-1:0.01:0)';
iq = imax* (0:0.01:1)';
nid = length(id);
te = zeros(nid,nid);
we = zeros(nid,nid);

for it3 = 1:nid
    % X.*Y denotes element-by-element multiplication.
    % here we calculate every id & iq and corresponding tq
    te(:,it3) = 1.5 *np *iq.*(psi_f+(Ld-Lq)*id(it3));
    %consider the Umax equation, calculate the max speed that can reach(depends on the load)
    we(:,it3) =umax ./ ( sqrt(  (Lq*iq).*(Lq*iq) + (psi_f+Ld*id(it3)).*(psi_f+Ld*id(it3))  ) );
end

wr = we2rpm(we);
telist = [ 84 210 40];
wrlist = [10400 6490 2285 3427 4000];
%pick up all the (id, iq) which results in the number of telist. ('ShowText','off')

% ----------------the trajectory on current circle calc--------------
theta = pi/2:0.01:pi;
id_ilim = imax.*cos(theta);
iq_ilim = imax.*sin(theta);
i_lim = [id_ilim',iq_ilim'];
T= 1.5 *np *iq_ilim.*(psi_f+(Ld-Lq).*id_ilim);
[M,I]=max(T);

% -------------------------------MTPA clac-----------------------

TE=0:1:M;

if slt_rat ==1
    id_mtpa = zeros(1,length(TE));
    iq_mtpa = linspace(-1.25*imax,1.25*imax,length(TE));
    
else
    %unitization
    ib = psi_f/beta/Ld;
    teb = 1.5*np*psi_f*ib;
    teN = TE/teb;
    
    h = teN.*teN;
    d = (12*h.*(768*h+81).^0.5-108*h).^(1/3);
    f = sqrt(9+6*d-288*h./d);
    idN = 0.75-f/12-sqrt(3-d+48*h./d+9./f)*sqrt(6)/12;
    iqN = teN./(1-idN);
    id_mtpa=idN*ib;
    iq_mtpa=iqN*ib;
end

[id_mtpa,iq_mtpa]=ilim(id_mtpa,iq_mtpa,imax);
mtpa_pts = unique([id_mtpa,iq_mtpa],'row');
id_mtpa = mtpa_pts(:,1);
iq_mtpa = mtpa_pts(:,2);
% scatter(id_mtpa,iq_mtpa);
% hold on

%------------------------- MTPV (MFPT)  calc------------------------
RPM = [2000:5:2500,2500:100:12000,12000:1000:15000];

if slt_rat ==1
    id_mtpv = zeros(1,length(RPM))-psi_f/Ld;
    iq_mtpv = linspace(0.3*imax,1.25*imax,length(RPM));
    
else
    
    WE=rpm2we(RPM);
    subst_1 = (Ld - Lq)/Lq./WE;
    Uq_minf = (psi_f - sqrt(psi_f^2 + umax^2*8*(subst_1.^2)) )./(-4 * subst_1);
    Ud_minf = sqrt(umax^2-Uq_minf.^2);
    id_mtpv = (Uq_minf./WE-psi_f)/Ld;
    iq_mtpv = Ud_minf/Lq./WE;
    
end

[id_mtpv,iq_mtpv]=ilim(id_mtpv,iq_mtpv,imax);
mtpv_pts = unique([id_mtpv,iq_mtpv],'row');
id_mtpv = mtpv_pts(:,1);
iq_mtpv = mtpv_pts(:,2);
mtpv_endx = min(id_mtpv);
mtpv_endy = max(iq_mtpv);
try
    mtpv_endrpm = pt2wr(mtpv_endx,mtpv_endy);
catch
    mtpv_endrpm = 99999;
end

%% map generator

if map_command == 1
    tq_axis = 5:5:M;
    %     rpm_axis = [5:50:mtpv_endrpm,mtpv_endrpm:100:6000,6000:500:8000,8000:1000:12000];
    rpm_axis = 5:100:14000;
    [~,rol_num]=size(tq_axis);
    [~,col_num]=size(rpm_axis);
    id_map = zeros(rol_num,col_num)*nan;
    iq_map = zeros(rol_num,col_num)*nan;
    rol_cont=1;
    wr_kp2_alt=99999;
    
    f = waitbar(0,'initial...','Name','Progress');
    set(f,'color','w');
    
    for it_tq = tq_axis
        ilim_mode = 0;
        waitbar(rol_cont/(rol_num+1),f,sprintf('%s',['Einen Moment bitte ... ',num2str(round(100*rol_cont/(rol_num+1))),'%']));
        col_cont=1;
        tq_temp=[it_tq,it_tq];
        tq_pts = contour(id, iq, te,tq_temp','ShowText','off','Color','y');
        hold on;
        tq_pts=(tq_pts(:,2:end))';
        id_consttq = tq_pts(:,1);
        iq_consttq = tq_pts(:,2);
        [kp1_x,kp1_y] = kreuz(id_consttq',iq_consttq',id_mtpa',iq_mtpa',0.1,'2-low');
        wr_kp1 = pt2wr(kp1_x,kp1_y);
        [kp2_x,kp2_y] = kreuz(id_consttq',iq_consttq',id_mtpv',iq_mtpv',0.1,'2-low');
        wr_kp2 = pt2wr(kp2_x,kp2_y);
        [kp3_x,kp3_y] = kreuz(id_consttq',iq_consttq',id_ilim',iq_ilim',0.1,'2-low');
        wr_kp3 = pt2wr(kp3_x,kp3_y);
        
        
        if kp2_x >= kp3_x      % for situation not reaching i limit, find all the pts > kp2 in mtpv
            phase1 = (mtpv_pts(:,1)>kp2_x).*mtpv_pts;
            
            if kp2_y>kp1_y
                phase1=[0,0];  % low tq area
            end
            phase1(all(phase1==0,2),:) = [];  %const_tq > k2 && const_tq<k1
            
            phase2 = ((tq_pts(:,1)>kp2_x) .* (tq_pts(:,1)<kp1_x)).*tq_pts;
            phase2(all(phase2==0,2),:) = [];
            tq_track = [phase1;phase2];        % conbine
            
        else
            ilim_mode = 1;
            phase1 = mtpv_pts;             % use all the mptv track
            %         scatter(phase1(:,1),phase1(:,2));
            %         hold on;
            
            phase2 = ((tq_pts(:,1)>kp3_x) .* (tq_pts(:,1)<kp1_x)).*tq_pts; % const_tq > k3 && const_tq < k1
            phase2(all(phase2==0,2),:) = [];
            %         scatter(phase2(:,1),phase2(:,2));
            %         hold on;
            
            phase3 = ((i_lim(:,1)>kp2_x) .* (i_lim(:,1)<kp3_x)).*i_lim; % ilim < k3 && ilim > k2
            phase3(all(phase3==0,2),:) = [];
            %         scatter(phase3(:,1),phase3(:,2));
            %         hold on;
            tq_track = [phase1;phase2;phase3];
        end
        %      scatter(tq_track(:,1),tq_track(:,2));
        %      hold on;
        for it_rev = rpm_axis
            if it_rev <= wr_kp1
                id_opt = kp1_x;
                iq_opt = kp1_y;
                if col_cont==1
                    scatter(id_opt,iq_opt,'filled');
                    hold on
                end
                
            elseif ((it_rev > wr_kp2_alt)&& (rol_cont>1)&&(kp2_y<kp1_y))||((ilim_mode==1)&&(it_rev > wr_kp3_alt))
                id_map(rol_cont,col_cont:end)=id_map(rol_cont-1,col_cont:end);
                iq_map(rol_cont,col_cont:end)=id_map(rol_cont-1,col_cont:end);
                break
            else
                rev_temp = [it_rev,it_rev];
                rev_pts = contour(id, iq, wr,rev_temp','ShowText','off','Color','g','Linewidth',0.1,'LineStyle',':');
                hold on;
                rev_pts = (rev_pts(:,2:end))';
                
                chra_pts = -psi_f/Ld; % To finish the compare, the ellips shall be devided into 2 parts in oder to get the unique condition
                maxrev_p1 = (rev_pts(:,1)>chra_pts).*rev_pts;
                maxrev_p1(all(maxrev_p1==0,2),:) = [];
                maxrev_p2 = (rev_pts(:,1)<=chra_pts).*rev_pts;
                maxrev_p2(all(maxrev_p2==0,2),:) = [];
                id_maxrev_p1 = maxrev_p1(:,1);
                iq_maxrev_p1 = maxrev_p1(:,2);
                id_maxrev_p2 = maxrev_p2(:,1);
                iq_maxrev_p2 = maxrev_p2(:,2);
                
                [iq_opt,id_opt]=kreuz((tq_track(:,2))',(tq_track(:,1))',iq_maxrev_p1',id_maxrev_p1',0.1,'Nan');  % ! when wr cross with tq , a y-x interp1 shall be used !
                if isnan(iq_opt)
                    [iq_opt,id_opt]=kreuz((tq_track(:,2))',(tq_track(:,1))',iq_maxrev_p2',id_maxrev_p2',0.1,'1-hi');
                end
                scatter(id_opt,iq_opt,'filled');
                hold on
            end
            id_map(rol_cont,col_cont)=id_opt;
            iq_map(rol_cont,col_cont)=iq_opt;
            col_cont=col_cont+1;
            %             scatter(id_opt,iq_opt,'filled');
            %             hold on
        end
        rol_cont=rol_cont+1;
        if kp2_y<kp1_y
            wr_kp2_alt = wr_kp2;
        else  % for low torque area
            wr_kp2_alt = it_rev;
        end
        wr_kp3_alt = wr_kp3;
    end
    close(f);
    save motordata.mat tq_axis rpm_axis id_map iq_map
    
end
%}


%% motor charactor figure drawing
%
if draw_flg == 1
    figure(1);
    set(figure(1),'position',[  5  361  560  420]);
    set(gca,'Box','on');
    title('E-mot opration points  (current)');
    axis equal;
    xlabel('id');
    ylabel('iq');
    hold on;
    %x-axis
    plot([-1.25*imax 1.25*imax],[0,0],'LineStyle','--','Color','k');
    %y-axis
    plot([0,0],[1.25*imax -1.25*imax],'LineStyle','--','Color','k');
    %critical point line
    plot([cp,cp],[1.25*imax -1.25*imax],'LineStyle',':','Color','k');
    
    %------------------id -iq -torque map------------------------------
    
    contour(id, iq, wr,wrlist','ShowText','off');
    hold on;
    contour(id,iq,te,telist','--','ShowText','on');
    hold on
    %-----------------draw some reference lines------------------------
    
    [id_exp_1,iq_exp_1]=example(id,iq,te, 419.8,'k');
    [id_exp_2,iq_exp_2]=example(id,iq,te,150,'k');
    
    % draw the current and voltage if needed
    % plot_s(-322.6,445.4,6000);
    
    % draw  current limit
    f = @(x,y) x.^2 + y.^2 - imax^2;
    fimplicit(f,[-1.25*imax 1.25*imax -1.25*imax 1.25*imax],'Color','r','Linewidth',2);
    hold on
    
    
    
    % draw power factor = 1
    z = @(x,y) (psi_f+Ld*x).*x+Lq*y.^2;
    fimplicit(z,[-1.25*imax 1.25*imax -1.25*imax 1.25*imax],'Color','g','Linewidth',1);
    hold on
    
    %---------------------draw MTPA , MTPV-----------------------------
    plot(id_mtpa,iq_mtpa,'Color','m','Linewidth',2);
    hold on;
    plot(id_mtpv,iq_mtpv,id_mtpv,-iq_mtpv,'Color','b','Linewidth',2);
    
    
    %----------------rpm - torque -power pic.----------------------------
    figure(2)
    title('speed VS torque');
    xlabel('rpm');
    ylabel('Nm');
    set(figure(2),'position',[  471  361  560  420])
    set(gca,'Box','on');
    hold on
    
    figure(3)
    set(figure(3),'position',[  981  361  560  420])
    set(gca,'Box','on');
    title('speed VS power');
    xlabel('rpm');
    ylabel('Power/KW');
    hold on
    
    trac(id_mtpv,iq_mtpv,'b');
    trac(id_mtpa,abs(iq_mtpa),'m');
    trac(id_ilim,iq_ilim,'r');
    trac(id_exp_1,iq_exp_1,'k');
    trac(id_exp_2,iq_exp_2,'k');
end
%}

%-------------------------fprintf ------------------------------
t = t + toc;

fprintf(' the character current is %f A \n',psi_f/Ld);

fprintf(' the MAX torque is %f N.m \n correspondind \n id = %f A \n iq = %f A \n\n',M, id_ilim(I),iq_ilim(I));

fprintf(' total time %f s \n',t);


%% function area

function rpm=we2rpm(w)
% ATTENTION electric w to rpm

global np;
rpm=w/2/pi*60/np;
end

function w=rpm2we(rpm)
% ATTENTION rpm to electric w
global np;
w=2*pi*rpm/60*np;
end

function trac(x,y,color)
% X Y is id,iq, T is torque can get from this id iq
global np psi_f Ld Lq umax ;
T= 1.5 *np *y.*(psi_f+(Ld-Lq).*x);
N =umax ./ ( sqrt(  (Lq*y).*(Lq*y) + (psi_f+Ld*x).*(psi_f+Ld*x)  ) );


figure(2);
N=we2rpm(N);
plot(N,T,'Color',color,'Linewidth',1);
hold on
figure(3);
plot(N,N.*T/9549,'Color',color);
hold on;
end

function [res1,res2]=example(x,y,z,v,color)
figure(1);
cnt=contour(x,y,z, 'LevelList',v,'Color',color,'Linewidth',1,'LineStyle','--');
hold on;
szc = size(cnt);
idz = 1;
while idz<szc(2)
    izi = cnt(2,idz);
    % subst the nonvalid infomation as some exist points
    cnt(2,idz) = cnt(2,end);
    cnt(1,idz) = cnt(1,end);
    idz = idz+izi+1;
end
res1=cnt(1,2:end);
res2=cnt(2,2:end);
end

function plot_s(x,y,w) %#ok
%plot speciffic id-iq and ud-uq, @ speed w
global Lq psi_f Ld np
w = 2*pi.*w/60*np;
X = -w*Lq*y;
Y = w*psi_f+w*Ld*x;

h= compass(x,y);
set(h,'Color','r','LineWidth',2,'LineStyle',':');
hold on
H= compass(X,Y);
set(H,'Color','b','LineWidth',2,'LineStyle',':');
hold on
end

function [res_x,res_y] = kreuz(L1x,L1y,L2x,L2y,aufloe,def)

if isempty(L2x)==1 || (isempty(L1x)==1) % if one is empty
    cros_flg = false;
    if isempty(L2x)~=1
        low=min(L2x);
        hi=max(L2x);
    elseif isempty(L1x)~=1
        low=min(L1x);
        hi=max(L1x);
    else
        res_x = nan;
        res_y = nan;
        return
    end
else  % both valid
    low = max(min(L1x),min(L2x));
    hi = min(max(L1x),max(L2x));
    % to find the crossing point of 2 lines
    init = interp1(L1x,L1y,low)-interp1(L2x,L2y,low);
    sign_alt = sign(init);
    cros_flg = false;
    for it_x = low:aufloe:hi
        comp1 = interp1(L1x,L1y,it_x);
        comp2 = interp1(L2x,L2y,it_x);
        if abs(comp1-comp2)<=eps
            cros_flg = true;
            res_x = it_x;
            res_y = comp1;
            break
        else
            sign_temp = sign(comp1-comp2);
            if sign_temp ~= sign_alt
                cros_flg = true;
                res_x = it_x-aufloe/2;
                res_y = (comp1+comp2)/2;
                break
            end
        end  % end of cross detection
    end  %end of for-loop
end% end of first if

if cros_flg == false
    switch def  % set the default value if no cross detected
        case '1-low'
            res_x =low;
            res_y =interp1(L1x,L1y,low);
        case '1-hi'
            res_x =hi;
            res_y =interp1(L1x,L1y,hi);
        case '2-low'
            res_x =low;
            res_y =interp1(L2x,L2y,low);
        case '2-hi'
            res_x =hi;
            res_y =interp1(L2x,L2y,hi);
        case 'Nan'
            res_x = nan;
            res_y = nan;
    end
end
end  % end of function

function [x_out,y_out] =ilim(x_in,y_in,imax)
checkpts = [x_in,y_in];
[rol,~] =  size(checkpts);
if rol == 1
    checkpts = [x_in',y_in'];
    [rol,~] =  size(checkpts);
end
for it_rol = 1:rol
    ck_x=checkpts(it_rol,1);
    ck_y=checkpts(it_rol,2);
    if ck_x^2+ck_y^2 > imax^2
        checkpts(it_rol,:)=NaN;
    elseif (ck_x>0||ck_x<(-imax)||(ck_y>imax)||(ck_y<0))
        checkpts(it_rol,:)=NaN;
    end
end
checkpts(any(isnan(checkpts)'),:) = [];
x_out = checkpts(:,1);
y_out = checkpts(:,2);
end

function rpm = pt2wr(pt_x,pt_y)
% to calculate the max rpm in this operation point.
global umax Lq psi_f Ld
we_temp = umax / ( sqrt(  (Lq*pt_y)^2 + (psi_f+Ld*pt_x)^2 ));
rpm = we2rpm(we_temp);
end
