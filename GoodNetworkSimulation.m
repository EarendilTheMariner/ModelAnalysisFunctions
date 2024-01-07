%%
clear
MyNetwork = NetworkParameters();

MyNetwork.AddCelltype(Celltype.V2a_2,'bi','','','','ipsi','',6000,0.1); %Exc
MyNetwork.AddCelltype(Celltype.V2a_1,'loc','','','flex','ipsi','',300,0.2);%Exc

MyNetwork.AddCelltype(Celltype.V3,'bi','','','ext','contra','',5000,0.2); % Exc
MyNetwork.AddCelltype(Celltype.DI6,'bi','','','','contra','',1500,1); %Inh

MyNetwork.AddCelltype(Celltype.V1,'bi','','','flex','ipsi','',1500,1); %Inh
MyNetwork.AddCelltype(Celltype.V2b,'bi','','','ext','ipsi','',2000,1); %Inh

MyNetwork.AddCelltype(Celltype.V0v,'bi','','','','contra','',4000,0.2); %Exc
MyNetwork.AddCelltype(Celltype.V0d,'loc','','','','contra','',300,1); %Inh

MyNetwork.AddCelltype(Celltype.MN,'loc','','','','ipsi','',200,0.2); %Exc
MyNetwork.AddCelltype(Celltype.V1r,'bi','','','','ipsi','',1500,1); %Inh

MyNetwork.AddCelltype('NF4','loc','','','ext','ipsi',{'5SpL','6SpL','7Sp','8Sp'},200,0.5) %Exc
MyNetwork.AddCelltype('NF5','loc','','','ext','ipsi',{'5SpM','6SpM','7Sp','Q9','Ad9'},200,0.5) %Exc 

MyNetwork.SetCellBiasPop('All','NF4',0);
MyNetwork.SetCellBiasPop('All','NF5',0);
MyNetwork.SetCellBiasPop('All','V1r',0);
MyNetwork.SetCellBiasPop('All','DI6',0);

MyNetwork.SetCellBiasPop('V0v','MN',0);
MyNetwork.SetCellBiasPop('V0d','MN',0);
MyNetwork.SetCellBiasPop('V0v','DI6',1);

MyNetwork.SetCellBiasPop('MN','All',0);
MyNetwork.SetCellBiasPop('MN','NF5',1);
MyNetwork.SetCellBiasPop('MN','NF4',1);
MyNetwork.SetCellBiasPop('MN','V3',1);
MyNetwork.SetCellBiasPop('MN','V1r',1);
MyNetwork.SetCellBiasPop('MN','DI6',0);

MyNetwork.SetCellBiasPop('V1r','All',0);
MyNetwork.SetCellBiasPop('V1r','MN',1);
MyNetwork.SetCellBiasPop('V1r','DI6',1);

MyNetwork.SetCellBiasPop('V3','V0v',0);
MyNetwork.SetCellBiasPop('V3','V0d',0);
MyNetwork.SetCellBiasPop('V3','DI6',1);

MyNetwork.SetCellBiasPop('DI6','V0v',0);
MyNetwork.SetCellBiasPop('DI6','V0d',0);
MyNetwork.SetCellBiasPop('DI6','DI6',1);

MyNetwork.SetCellBiasPop('V2a-2','V3',0);
MyNetwork.SetCellBiasPop('V2a-1','V3',0);
MyNetwork.SetCellBiasPop('V2a-2','DI6',1);
MyNetwork.SetCellBiasPop('V2a-1','DI6',1);

N = Network(['T12','T13', append('L',string([1:6])),append('S',string([1:4]))],0.75,MyNetwork,'Balance',1);
%%
t_steps = 6000;
threshold = ones(size(N.ConnMat,1),t_steps)*20;
threshold(N.Types == Celltype.NF5,:) = 5;
threshold(N.Types == Celltype.NF4,:) = 5;

I_e = ones(size(N.ConnMat,1),t_steps)*20;

gain = ones(size(N.ConnMat,1),t_steps)*0.95;
N.Simulate(t_steps,'I_e',I_e,'threshold',threshold,'gain',gain);

%%
%load("03_Utilities\ChannelMaps\neuropixPhase3B2_kilosortChanMap.mat");
%[NOI,ElectrodePos] = PositionElectrodeAndSample3D([xcoords,ycoords],[400 10000 0],-2200,MouseNet.Position,'rotate',[35,0,0],'DistThresh',1000);
%%

Diff = GeneratePlot(N);
%% Convolves Kernel with sine modulated impulse trains and cosine modulated impulse 
figure
for k = 0:0.1:10
    combsin = zeros(3000,1);
    combcos = zeros(3000,1);
    pad = combsin;
    combsin(1:1:end) = sin(k*2*pi*(1:length(combsin))/(length(combsin)));%.*(1:-1/(length(combsin)-1):0);
    combsin = [pad;combsin;pad];
    combcos(1:1:end) = cos(k*2*pi*(1:length(combcos))/(length(combcos)));%.*(1:-1/(length(combcos)-1):0);
    combcos = [pad;combcos;pad];
    sigsin = movmean(conv(combsin,Diff),50);
    sigcos = movmean(conv(combcos,Diff),50);
    sigtps = sigsin(((length(Diff)/2)+length(pad)):(length(sigsin)-(length(Diff)/2+length(pad))));
    sigtpc = sigcos(((length(Diff)/2)+length(pad)):(length(sigcos)-(length(Diff)/2+length(pad))));

    %polarplot((0:2*pi/length(sigtp):2*pi)',[nan; sigtp]);
    subplot(1,2,1)
    plot(sigtps,1:length(sigtps))
    hold on
    subplot(1,2,2)
    plot(sigtpc,1:length(sigtpc));
    hold on

end
%vline(1)

%% Pre Synaptic Neurons
Pop = 'NF5';

% Yellow = Pop
% Blue = Pre/post

figure
%Kernel function
L = round(range(N.Position(:,2)),-2);
edges = [-L:10:L];
plot(edges(1:end-1),Diff);
grid();
title('Kernel function');
xlabel('Position from soma (µm)');
ylabel('Exc/Inh diff');


figure 
ixoi = contains(string(N.Types),Pop) ;
Pre = sum(N.ConnMat(ixoi' & (N.Latera > 0)',:),1);
Cioi = logical(Pre); %& (N.Transmit < 0)' ;
scatter3(N.Geometry.Position(N.Geometry.Type=="WM",1),N.Geometry.Position(N.Geometry.Type=="WM",2),N.Geometry.Position(N.Geometry.Type=="WM",3),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none','MarkerFaceColor',Colors().BergGray02);
hold on 
scatter3(N.Position(Cioi,1),N.Position(Cioi,2),N.Position(Cioi,3),abs(Pre(Cioi))*5,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none','MarkerFaceColor',Colors().BergBlue);
scatter3(N.Position(ixoi'&(N.Latera>0)',1),N.Position(ixoi'&(N.Latera>0)',2),N.Position(ixoi'&(N.Latera>0)',3),20,"filled",'diamond')
axis equal
title('Presynaptic Neurons Distribution')

% Post Synaptic Neurons
figure 
ixoi = contains(string(N.Types),Pop) ;
Post = sum(N.ConnMat(:,ixoi' & (N.Latera > 0)'),2);
Cioi = logical(Post) ;
scatter3(N.Geometry.Position(N.Geometry.Type=="WM",1),N.Geometry.Position(N.Geometry.Type=="WM",2),N.Geometry.Position(N.Geometry.Type=="WM",3),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeColor','none','MarkerFaceColor',Colors().BergGray02);
hold on 
scatter3(N.Position(Cioi,1),N.Position(Cioi,2),N.Position(Cioi,3),abs(Post(Cioi))*5,'filled','MarkerFaceAlpha',0.7,'MarkerEdgeColor','none','MarkerFaceColor',Colors().BergBlue);
scatter3(N.Position(ixoi'&(N.Latera>0)',1),N.Position(ixoi'&(N.Latera>0)',2),N.Position(ixoi'&(N.Latera>0)',3),20,"filled",'diamond')
axis equal
title('Postsynaptic Neurons Distribution')


%%
function difference = GeneratePlot(N)
screensize = get(groot,'ScreenSize');
ConnMat = N.ConnMat;
Ex =  N.Transmit > 0;
In =  N.Transmit < 0;
Dist = bsxfun(@minus,N.Position(:,2),N.Position(:,2)');
ProjEx = ConnMat(:,Ex); 
ProjIn = ConnMat(:,In);
DistEx = Dist(:,Ex);
DistIn = Dist(:,In);
[vec,bevr] = eig(ConnMat,'vector');
vec = real(vec);
[brevr,bixev] = sort(real(bevr),'descend');
bievr = imag(bevr(bixev));
[~,ind] = unique(brevr,'stable');
L = round(range(N.Position(:,2)),-2);
edges = [-L:10:L];
yin = discretize(DistIn(ProjIn~= 0), edges);
yex = discretize(DistEx(ProjEx~= 0), edges);
[GnIn,Min,Cin] = grpstats(ProjIn(ProjIn~=0), yin,["gname","mean","numel"]);
[GnEx,Mex,Cex] = grpstats(ProjEx(ProjEx~=0), yex,["gname","mean","numel"]);
MInZ = zeros(length(edges)-1,1);
MExZ = zeros(length(edges)-1,1);
CInZ = zeros(length(edges)-1,1);
CExZ = zeros(length(edges)-1,1);
MInZ(str2double(GnIn)) = Min;
MExZ(str2double(GnEx)) = Mex;
CInZ(str2double(GnIn)) = Cin;
CExZ(str2double(GnEx)) = Cex;
fig = figure(Position=screensize);
subplot (2,3,1)
bar(edges(1:end-1), abs(MInZ.*CInZ),'FaceColor',[Colors().BergOrange],'FaceAlpha',0.5);
hold on
bar(edges(1:end-1), abs(MExZ.*CExZ),'FaceColor',[Colors().BergElectricBlue],'FaceAlpha',0.5);
legend({'Inhib','Excit'});
set(gca, 'xtick', round(edges(1:end-1),0));
set(gca, 'xticklabels',round(edges(1:end-1),0));
xlabel("Distance from Soma [um]");
difference = abs(MExZ.*CExZ)-abs(MInZ.*CInZ);
difplus = difference;
difneg = difference;
difplus(difplus < 0) = 0;
difneg(difneg > 0) = 0;
difneg = abs(difneg);
col = [normalize(difneg,'range').*Colors().BergOrange+normalize(difplus,'range').*Colors().BergElectricBlue];
subplot (2,3,2)
b = bar(edges(1:end-1),abs(MExZ.*CExZ)-abs(MInZ.*CInZ),'FaceColor','flat');
for k = 1:size(MExZ,1)
    b.CData(k,:) = col(k,:);
end
set(gca, 'xtick', round(edges(1:100:end-1),0));
set(gca, 'xticklabels',round(edges(1:100:end-1),0));
legend({'Diff Excit Inihib'});
xlabel("Distance from Soma [um]");
subplot (2,3,3)
scatter(brevr,bievr,'filled','MarkerFaceColor',Colors().BergBlack,'MarkerFaceAlpha',0.5)
vline(1)
axis equal
axis([-2 2 -2 2])
subplot(2,3,4:6)
plot(N.Rates(1:end,:));
end