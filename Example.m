% Creating the signal Ex2_Nonstationary_NonSin_waves_plus_trend
%
%  Ref. Example 2 in A. Cicone, H. Zhou. 'Multidimensional Iterative Filtering method 
%      for the decomposition of high-dimensional non-stationary signals'.
%      Preprint ArXiv http://arxiv.org/abs/1507.07173
% 

z2=[];

iMax=8;
for i=0:iMax
    T=2-2*i/10;
    t=0:0.02:T;
    z21=t/T;
    z22=fliplr(z21);
    if not(i==iMax) && not(i==0) 
        z2t=[z21 z22(2:end-1)];
    elseif i==iMax
        z2t=[z21 z22(2:end)];
    else
        z2t=[z21(ceil(end/2):end) z22(2:end-1)];
    end
    z2=[z2 z2t];    
end
z2t=fliplr(z2);
z2=[z2 z2t(2:end)];
z2((length(z2)+1)/2)=0;
z2=z2-0.5;

%%
f1=plot_2D_v1(z2((length(z2)+1)/2:end),200);

%%

t2=linspace(0,10,length(z2));


z3=sin(pi*0.7*t2);

%%

f2=plot_2D_v1(z3((length(z3)+1)/2:end),200);

%%

x=linspace(0,10+1/2,200*2+1);
f3=0.5*x.'*ones(1,length(x));

%%

f=f1+f2+f3;
%% 
figure
h=surf(f);
set(h, 'edgecolor','none')
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar
set(gca,'fontsize', 25);
axis([1 size(f,1) 1 size(f,2) floor(min(min(f))) ceil(max(max(f)))])

%% Section of the signal in the middle
figure
plot((f(:,(end+1)/2)),'k','Linewidth',2)
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(f,1) floor(min(min(f))) ceil(max(max(f)))])

%% First IMF
close all
clc

opts=decompSettings_MIF_2D_v01('plots',0,'MIF.delta',0.01);
[IMF,SDlog,havelog] = Decomp_MIF_2D_v10(f,opts,1);

figure
h=surf(IMF(:,:,1));
set(h, 'edgecolor','none')
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar
set(gca,'fontsize', 25);
axis([1 size(IMF(:,:,1),1) 1 size(IMF(:,:,1),2) floor(min(min(IMF(:,:,1)))) ceil(max(max(IMF(:,:,1))))])


%% Section of the IMF in the middle
figure
plot((f1(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF(:,(end+1)/2,1)),'k','Linewidth',2)
legend('Ground truth','IMF_1')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF(:,:,1),1) floor(min(min(IMF(:,:,1)))) ceil(max(max(IMF(:,:,1))))])


%% Section of the Remainder along the diagonal
figure
plot((f2(:,(end+1)/2)+f3(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF(:,(end+1)/2,2)),'k','Linewidth',2)
legend('Ground truth','Remainder')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF(:,:,2),1) floor(min(min(IMF(:,:,2)))) ceil(max(max(IMF(:,:,2))))])

%% Second IMF
close all

opts=decompSettings_MIF_2D_v01('plots',0,'MIF.delta',0.01);
[IMF2,SDlog1,havelog1] = Decomp_MIF_2D_v10(IMF(:,:,end),opts,1);

%% IMF using mask length (1000x2+1)

figure
h=surf(IMF2(:,:,1));
set(h, 'edgecolor','none')
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar
set(gca,'fontsize', 25);
axis([1 size(IMF2(:,:,1),1) 1 size(IMF2(:,:,1),2) floor(min(min(IMF2(:,:,1)))) ceil(max(max(IMF2(:,:,1))))])

%% Section of the IMF in the middle
figure
plot((f2(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF2(:,(end+1)/2,1)),'k','Linewidth',2)
legend('Ground truth','IMF_2')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF2(:,:,1),1) floor(min(min(IMF2(:,:,1)))) ceil(max(max(IMF2(:,:,1))))])

%% Section of the Remainder along the diagonal
figure
plot((f3(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF2(:,(end+1)/2,2)),'k','Linewidth',2)
legend('Ground truth','Remainder')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF2(:,:,2),1) floor(min(min(IMF2(:,:,2)))) ceil(max(max(IMF2(:,:,2))))])



%% Difference between the IMF and the original nonstationary oscillations

h=surf(IMF(:,:,1)-f1);
set(h, 'edgecolor','none')


%% Difference between the obtained trend and the trend plus bump we add in the signal

h=surf(IMF(:,:,1)-f1+IMF2(:,:,1)-f2);
set(h, 'edgecolor','none')







