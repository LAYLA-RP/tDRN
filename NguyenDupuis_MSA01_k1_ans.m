%% MSA��������k�ĵ�����Ϊ����Ȩ�أ����_k�޸���ͼ����ʾ������function G0_Plot(G_temp),�޸�function G_FlowPlot(G_temp)���޸�%����Ķ�
%����od�������ӻ���(G=G_append_OD(G0)���ڸĽ�G0)û���ϣ�ֱ�Ӹ�node��line������
% ���ӷֳ���
%% ·��������������
clear;clc
%%  step1�� ·����������¼��
% GIS���excel���ݵ��벢������ͼG�ĸ�ʽ�洢;
[EdgeTable,NodeTable]=topotable();%���˽ṹ�������ݵ��룻% load topo_basic  %%֮ǰ�ֶ��������
% ���ɻ�������·��
G0=digraph(EdgeTable,NodeTable);%node������name�еĻ���������ʾnode��name% G=digraph(P);
% G0=G_append_OD(G0);%����OD��ʾ��û�б�Ҫ��ֻҪԭʼ��ͼ����OD�ڵ�Ϳ���
%%  step2����ͨ�������ʱ��t�� OD����¼��.�ֶ�¼���Ǳ��뱣֤ԭ�б��node������˳��
[OD_input_p,OD_input_f,OD_struct] = OD_matr(2005,10,NodeTable,380,930);%ȱ������,�����Զ�����Ϊ�˼�������ʵ����
% OD_input=OD_struct(1,1).volume;
% [mon_t]=daycal(year);
% for i=1:12
%     [OD_struct,OD_input_p,OD_pair] = OD_matr(year,mon_t(i,1));
% end
%% step4��������־�������¼�룬ͬ�ϣ���ɱ�񲢶�ȡ
%�������룺ʩ����ͨ��֯������ż���Ӧ����
%�ƻ�ʩ��ʱ�估��Ӧ���������λ�ã�׮�ţ�/������շ�������
%���������ʱ����·������ʱ���
% �ź�����matrx_mon
% matrx_mon=timeseries([0:11],[0:11]);
%net=xlsread('maintain.xlsx','A2:C100');  
%��Ӧʱ�̵���һ�������ع��ؼ�λ�ú���
% T=3;%PCI���²�ѯ����/��
%��Σ���ʱ�ν��棩�����ط���
%t:ά�ֵ�ʱ��;programme:��������(0/1/2�ֱ��ʾδ��������/����A����/����B����);paremater:��������Ҫ�ṩ�Ĳ���
[time_event,parameter]=xlsread('maintain.xlsx','A2:D100'); 
%% ��ʼ״̬
erfa1=0;%���гɱ�Ԥ��ϵ��,�ú������Է���erfa�Ե���Ч���Ƿ�������
theta_input=1.6;
beita1=1;epsilon = 0.0001;%���ó�ʼ����������
[H0{1}] =MSA_G(theta_input,G0,OD_input_p,erfa1,beita1,epsilon);
%��ʼ״̬ͼ
% G0_Plot(G0)
% G_FlowPlot1(H0)
save G_basic0_2%�洢�����������
%% ����·�������ط���
% G_update=update_topu_meantime(time_event,parameter,G0);
% G_temp=MSA_G(theta_input,G_update,OD_input_p);%���˽ṹ��������Ľ�ͨ�ط���
% G_flow=flow_lane(time_event,G_temp);%time_event:t/programme(0��ʾû����������/1��2�ֱ��ʾ���ַ���)

%% ��������
clear;clc
load G_basic0_2.mat
x=1;%������ţ������Լ�����ʱʶ��ͼ
theta_input=1.6;
erfa=0.8;%���гɱ�Ԥ��ϵ��,�ú������Է���erfa�Ե���Ч���Ƿ�������
beita=0.8;epsilon = 0.0001;%���ó�ʼ����������
[G_LaneProperty,paramater_plus,H_MSA_Property]=G_dup(G0,time_event,parameter,theta_input,OD_input_p,erfa,beita,epsilon );
save G_basic_7

%% ·����ʷ��������ͳһ����·�β������������
clear;clc;
load G_basic_7.mat
[G_Flow,G_Property_Separator]=UpdateFlow_LinkSeparator3(G_LaneProperty,paramater_plus,time_event);%��ͬʱ��·��·�λ���ͳһ��������ȡ�������������ں������(Xi��PCI)
G_FlowSum=G_FlowSum_MSA(12,G_Flow,time_event);%����PCI�������ڣ�������Ӧ��·�������ܺͣ�eg.PCI��������Ϊ3m����ôFlowSum��������ҲΪ3m
save G_basic_8

%% �����ͼ
clear;clc;
load G_basic_8.mat
%% ����ͼ���ӻ�
G0_Plot(G0)
%ԭʼ�ṹ��̬������
% G_FlowPlot1(H0)
% G_FlowPlot2(H0)
% G0_Plot(G_LaneProperty{2})
% G_FlowPlot1(G_Flow)
G_FlowPlot_T1(G_LaneProperty)
G_LaneProperty1=G_LaneProperty;
G_LaneProperty1{2}.Edges.X1i(39)=762;
G_LaneProperty1{2}.Edges.X2i(39)=762;
G_LaneProperty1{2}.Edges.X1i(41)=762;
G_LaneProperty1{2}.Edges.X2i(41)=762;
G_FlowPlot_T1(G_LaneProperty1)
G_FlowPlot1(G_FlowSum)
% G_FlowPlot01(H_MSA_Property)
% G_FlowPlot2(G_LaneProperty)%����೵��ͬʱ��ʾ
% G_FlowPlot1(G_FlowSum)

% nn=length(G_FlowSum);
% m=size(G_FlowSum{1}.Edges,1);
% FlowSum1=zeros(m,nn);
% FlowSum2=zeros(m,nn);
% % figure;
% % for i=1:nn
% %     subplot(2,3,i);
% %     plot(G_FlowSum{i});
% % end
% for j=1:nn
%     FlowSum1(:,j)=G_FlowSum{j}.Edges.X1i;
%     FlowSum2(:,j)=G_FlowSum{j}.Edges.X2i;
% end
% % subplot(2,1,1);
% % bar3(FlowSum1)
% % 
% % subplot(2,1,2);
% % bar3(FlowSum2)

%% Function
%% ·�����ݸ���
function [G_LaneProperty,paramater_plus,H_MSA_Property]=G_dup(G0,time_event,parameter,theta_input,OD_input_p,erfa1,beita1,epsilon)%·�����˸���-��������
n=size(time_event,1);
H_TopuUpdate=cell(n,1);
H_MSA_Property=cell(n,1);
G_LaneProperty=cell(n,1);
paramater_plus=cell(n,1);
for i=1:n
    if i==1
        G_temp=G0;
    else
        G_temp=G_LaneProperty{i-1};
    end
    parameter_temp=str2num(parameter{i,1});
    time_event_temp=time_event(i,:);
    if time_event_temp(2)==0
        H_TopuUpdate{i}=G_temp;%���˱���Ĵ洢
        H_MSA_Property{i}=MSA_G(theta_input,H_TopuUpdate{i},OD_input_p,erfa1,beita1,epsilon);%���˽ṹ��������Ľ�ͨ�ط��䣬ֻ���ڼ��㣬���Բ��������
        G_LaneProperty{i}=H_MSA_Property{i};%·���������Ը���
    elseif time_event_temp(2)==1
        G_save_temp=GraphChange_A(G_temp,parameter_temp);%[18,19,300,500]����A��·�ξֲ���գ�����ָ����򿪷���������;
        H_TopuUpdate{i}=G_save_temp{1};
        paramater_plus{i}=G_save_temp{2};
        H_MSA_Property{i}=MSA_G(theta_input,G_save_temp{1},OD_input_p,erfa1,beita1,epsilon);%���˽ṹ��������Ľ�ͨ�ط���
        G_LaneProperty{i}=G_LaneProperty_Change_A(H_MSA_Property{i},G_temp,parameter_temp);%[18,19,300,500]����A
    elseif time_event_temp(2)==2
         parameter_temp_node=str2num(parameter{i,2});
         G_save_temp=GraphChange_B(G_temp,parameter_temp,parameter_temp_node(1,1));%[27,18;34,29],28����B���ѵ����%����A��·�ξֲ���գ�����ָ����򿪷���������;
         H_TopuUpdate{i}=G_save_temp{2};
         H_MSA_Property{i}=MSA_G(theta_input,G_save_temp{2},OD_input_p,erfa1,beita1,epsilon);%���˽ṹ��������Ľ�ͨ�ط���
         G_LaneProperty{i}=G_LaneProperty_Change_B(H_MSA_Property{i},G_temp,G_save_temp{4},G_save_temp{3});%����B
    end
end
end

%% ���ӻ�function
%G0���ƣ����OD��
function G0_Plot(G_temp)
G=G_temp;
XData0=G.Nodes.XDate;
YData0=G.Nodes.YDate;
nn=length(G.Nodes.YoN);
j=0;
for i=1:nn
    if G.Nodes.YoN(i)==1
        j=j+1;
        hightlight_node(j)=i;%hightlight node��ȡ
%         NodeCData(i)=1;
        MarkerSize(i)=6;
        Markercell{i}='o';
    else
%         NodeCData(i)=0;
        MarkerSize(i)=4;
        Markercell{i}='s';
    end
end
figure
P=plot(G);
% P.NodeCData=NodeCData;
Marker=string(Markercell);
P.MarkerSize=MarkerSize;
P.Marker=Marker;
% P.NodeColor='k';
highlight(P,hightlight_node,'NodeColor','k')%ͻ����ʾ
str0=['��ʼ·���ṹ'];
title(str0)
P.XData=XData0';
P.YData=YData0';
% colorbar
end      

%����pcu��·���ķֲ�G_FlowPlot
function G_FlowPlot02(G_temp)%G_FlowPlot2������ȥ����
n=length(G_temp);
max_Xi=[0,0];
max_Xi_temp=zeros(n,2);
for i=n
    G=G_temp{i};
    for j=1:2
        if j==1
            max_Xi_temp(i,j)=max(G.Edges.X1i);
        elseif j==2
            max_Xi_temp(i,j)=max(G.Edges.X2i);
        end
        if max_Xi(j)<max_Xi_temp(i,j)
            max_Xi(j)=max_Xi_temp(i,j);
        end
    end
end
for i=1:n
    G=G_temp{i};
    XData0=G.Nodes.XDate;
    YData0=G.Nodes.YDate;
    X1i=G.Edges.X1i;
    X1i(X1i(:,1)==0,1)=0.0001;
    X2i=G.Edges.X2i;
    X2i(X2i(:,1)==0,1)=0.0001;
    LWidths{1} = 10*X1i/max_Xi(1);
    LWidths{2} = 10*X2i/max_Xi(2);
%     EdgeCData{1} = X1i/max_Xi(1);
%     EdgeCData{2} = X2i/max_Xi(2);
    EdgeCData{1} = X1i;
    EdgeCData{2} = X2i;
    nn=length(G.Nodes.YoN);
    jj=0;
    for j=1:nn
        if G.Nodes.YoN(j)==1
            jj=jj+1;
            hightlight_node(jj)=j;%hightlight node��ȡ
    %         NodeCData(j)=1;
            MarkerSize(j)=6;
            Markercell{j}='o';
        else
    %         NodeCData(j)=0;
            MarkerSize(j)=4;
            Markercell{j}='s';
        end
    end
    P=cell(1,2);
    f=figure;
    f.Position= [680 -50 660 2200];
    for j=1:2
        subplot(2,1,j)
        P{j}=plot(G);
        P{j}.EdgeLabel=G.Edges.X1i;
        P{j}.EdgeCData=EdgeCData{j};
        P{j}.LineWidth=LWidths{j};
        Marker=string(Markercell);
        P{j}.MarkerSize=MarkerSize;
        P{j}.Marker=Marker;
%         P{j}.NodeColor='k';
        highlight(P{j},hightlight_node,'NodeColor','k')
%         P{j}.NodeCData=NodeCData;
        if j==1
            P{j}.EdgeLabel=G.Edges.X1i;
            str1=['��',num2str(i),'�θ���·����೵��·��Сʱpcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['��',num2str(i),'�θ���·���ڲ೵��·��Сʱpcu'];
        end
            title(str1)
%             P{j}.XData=XData0';
%             P{j}.YData=YData0';
            colorbar
    end
end
end
function G_FlowPlot1(G_temp)%ҳ������ʾ1������
n=length(G_temp);
max_Xi=[0,0];
max_Xi_temp=zeros(n,2);
for i=n
    G=G_temp{i};
    for j=1
        if j==1
            max_Xi_temp(i,j)=max(G.Edges.X1i);
        elseif j==2
            max_Xi_temp(i,j)=max(G.Edges.X2i);
        end
        if max_Xi(j)<max_Xi_temp(i,j)
            max_Xi(j)=max_Xi_temp(i,j);
        end
    end
end
for i=1:n
    G=G_temp{i};
    XData0=G.Nodes.XDate;
    YData0=G.Nodes.YDate;
    X1i=G.Edges.X1i;
    X1i(X1i(:,1)==0,1)=0.0001;
    X2i=G.Edges.X2i;
    X2i(X2i(:,1)==0,1)=0.0001;
    LWidths{1} = 10*X1i/max_Xi(1);
    LWidths{2} = 10*X2i/max_Xi(2);
%     EdgeCData{1} = X1i/max_Xi(1);
%     EdgeCData{2} = X2i/max_Xi(2);
    EdgeCData{1} = X1i;
    EdgeCData{2} = X2i;
    nn=length(G.Nodes.YoN);
    jj=0;
    for j=1:nn
        if G.Nodes.YoN(j)==1
            jj=jj+1;
            hightlight_node(jj)=j;%hightlight node��ȡ
    %         NodeCData(j)=1;
            MarkerSize(j)=6;
            Markercell{j}='o';
        else
    %         NodeCData(j)=0;
            MarkerSize(j)=4;
            Markercell{j}='s';
        end
    end
    P=cell(1,2);
    f=figure;
    f.Position= [480 0 1120 840];
    for j=1
%         subplot(2,1,j)
        P{j}=plot(G);
        P{j}.EdgeCData=EdgeCData{j};
        P{j}.LineWidth=LWidths{j};
        Marker=string(Markercell);
        P{j}.MarkerSize=MarkerSize;
        P{j}.Marker=Marker;
        P{j}.NodeLabelColor='r';
        P{j}.ArrowSize=10;
        P{j}.EdgeFontSize=10;
%         P{j}.NodeColor='k';
        highlight(P{j},hightlight_node,'NodeColor','r')
        P{j}.NodeLabelColor='r';
%         P{j}.NodeCData=NodeCData;
        if j==1
            P{j}.EdgeLabel=G.Edges.X1i;
            str1=['��',num2str(i),'�θ���·����೵��·��Сʱpcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['��',num2str(i),'�θ���·���ڲ೵��·��Сʱpcu'];
        end
            title(str1)
            P{j}.XData=XData0';
            P{j}.YData=YData0';
            colorbar
    end
end
end
function G_FlowPlot01(G_temp)%G_FlowPlot1������ȥ����
n=length(G_temp);
max_Xi=[0,0];
max_Xi_temp=zeros(n,2);
for i=n
    G=G_temp{i};
    for j=1
        if j==1
            max_Xi_temp(i,j)=max(G.Edges.X1i);
        elseif j==2
            max_Xi_temp(i,j)=max(G.Edges.X2i);
        end
        if max_Xi(j)<max_Xi_temp(i,j)
            max_Xi(j)=max_Xi_temp(i,j);
        end
    end
end
for i=1:n
    G=G_temp{i};
%     XData0=G.Nodes.XDate;
%     YData0=G.Nodes.YDate;
    X1i=G.Edges.X1i;
    X1i(X1i(:,1)==0,1)=0.0001;
    X2i=G.Edges.X2i;
    X2i(X2i(:,1)==0,1)=0.0001;
    LWidths{1} = 10*X1i/max_Xi(1);
    LWidths{2} = 10*X2i/max_Xi(2);
%     EdgeCData{1} = X1i/max_Xi(1);
%     EdgeCData{2} = X2i/max_Xi(2);
    EdgeCData{1} = X1i;
    EdgeCData{2} = X2i;
    nn=length(G.Nodes.YoN);
    jj=0;
    for j=1:nn
        if G.Nodes.YoN(j)==1
            jj=jj+1;
            hightlight_node(jj)=j;%hightlight node��ȡ
    %         NodeCData(j)=1;
            MarkerSize(j)=6;
            Markercell{j}='o';
        else
    %         NodeCData(j)=0;
            MarkerSize(j)=4;
            Markercell{j}='s';
        end
    end
    P=cell(1,2);
    f=figure;
    f.Position= [480 0 1120 840];
    for j=1
%         subplot(2,1,j)
        P{j}=plot(G);
        P{j}.EdgeCData=EdgeCData{j};
        P{j}.LineWidth=LWidths{j};
        Marker=string(Markercell);
        P{j}.MarkerSize=MarkerSize;
        P{j}.Marker=Marker;
        P{j}.NodeLabelColor='r';
        P{j}.ArrowSize=10;
        P{j}.EdgeFontSize=10;
%         P{j}.NodeColor='k';
        highlight(P{j},hightlight_node,'NodeColor','r')
        P{j}.NodeLabelColor='r';
%         P{j}.NodeCData=NodeCData;
        if j==1
            P{j}.EdgeLabel=G.Edges.X1i;
            str1=['��',num2str(i),'�θ���·����೵��·��Сʱpcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['��',num2str(i),'�θ���·���ڲ೵��·��Сʱpcu'];
        end
            title(str1)
%             P{j}.XData=XData0';
%             P{j}.YData=YData0';
            colorbar
    end
end
end
function G_FlowPlot2(G_temp)%���֣�ҳ������ʾ2������
n=length(G_temp);
max_Xi=[0,0];
max_Xi_temp=zeros(n,2);
for i=n
    G=G_temp{i};
    for j=1:2
        if j==1
            max_Xi_temp(i,j)=max(G.Edges.X1i);
        elseif j==2
            max_Xi_temp(i,j)=max(G.Edges.X2i);
        end
        if max_Xi(j)<max_Xi_temp(i,j)
            max_Xi(j)=max_Xi_temp(i,j);
        end
    end
end
for i=1:n
    G=G_temp{i};
    XData0=G.Nodes.XDate;
    YData0=G.Nodes.YDate;
    X1i=G.Edges.X1i;
    X1i(X1i(:,1)==0,1)=0.0001;
    X2i=G.Edges.X2i;
    X2i(X2i(:,1)==0,1)=0.0001;
    LWidths{1} = 10*X1i/max_Xi(1);
    LWidths{2} = 10*X2i/max_Xi(2);
%     EdgeCData{1} = X1i/max_Xi(1);
%     EdgeCData{2} = X2i/max_Xi(2);
    EdgeCData{1} = X1i;
    EdgeCData{2} = X2i;
    nn=length(G.Nodes.YoN);
    jj=0;
    for j=1:nn
        if G.Nodes.YoN(j)==1
            jj=jj+1;
            hightlight_node(jj)=j;%hightlight node��ȡ
    %         NodeCData(j)=1;
            MarkerSize(j)=6;
            Markercell{j}='o';
        else
    %         NodeCData(j)=0;
            MarkerSize(j)=4;
            Markercell{j}='s';
        end
    end
    P=cell(1,2);
    f=figure;
    f.Position= [680 -50 560 1000];
    for j=1:2
        subplot(2,1,j)
        P{j}=plot(G);
        P{j}.EdgeCData=EdgeCData{j};
        P{j}.LineWidth=LWidths{j};
        Marker=string(Markercell);
        P{j}.MarkerSize=MarkerSize;
        P{j}.Marker=Marker;
        P{j}.NodeLabelColor='r';
        P{j}.ArrowSize=10;
%         P{j}.NodeColor='k';
        highlight(P{j},hightlight_node,'NodeColor','r')
%         P{j}.NodeCData=NodeCData;
        if j==1
            P{j}.EdgeLabel=G.Edges.X1i;
            str1=['��',num2str(i),'�θ���·����೵��·��Сʱpcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['��',num2str(i),'�θ���·���ڲ೵��·��Сʱpcu'];
        end
            title(str1)
            P{j}.XData=XData0';
            P{j}.YData=YData0';
            colorbar
    end
end
end
%����=������
function G_FlowPlot_T1(G_temp)%ҳ������ʾn��������
n=length(G_temp);
max_Xi=0;
max_Xi_temp=zeros(n,1);
for i=n
    G=G_temp{i};
    Xi=G.Edges.X1i+G.Edges.X2i;
    max_Xi_temp(i,1)=max(sum(Xi,2));
    if max_Xi<max_Xi_temp(i,1)
        max_Xi=max_Xi_temp(i,1);
    end
end
for i=1:n
    G=G_temp{i};
    XData0=G.Nodes.XDate;
    YData0=G.Nodes.YDate;
    X1i=sum(G.Edges.X1i,2);
    X2i=sum(G.Edges.X2i,2);
    Xi=X1i+X2i;
    Xi(Xi(:,1)==0,1)=0.0001;    
    LWidths= 10*Xi/max_Xi;
%     EdgeCData{1} = X1i/max_Xi(1);
%     EdgeCData{2} = X2i/max_Xi(2);
    EdgeCData= X1i+X2i;
    nn=length(G.Nodes.YoN);
    jj=0;
    for j=1:nn
        if G.Nodes.YoN(j)==1
            jj=jj+1;
            hightlight_node(jj)=j;%hightlight node��ȡ
    %         NodeCData(j)=1;
            MarkerSize(j)=3;
            Markercell{j}='o';
        else
    %         NodeCData(j)=0;
            MarkerSize(j)=1;
            Markercell{j}='s';
        end
    end
    P=cell(1,2);
    f=figure;
    f.Position= [480 0 1120 840];
    
        P=plot(G);
        P.EdgeCData=EdgeCData;
        P.LineWidth=LWidths;
        Marker=string(Markercell);
        P.MarkerSize=MarkerSize;
        P.Marker=Marker;
        P.NodeLabelColor='r';
        P.ArrowSize=10;
        P.EdgeFontSize=10;
        P.NodeFontSize=10;
%         P.NodeColor='k';
        highlight(P,hightlight_node,'NodeColor','r')
        P.NodeLabelColor='r';
%         P.NodeCData=NodeCData;
        P.EdgeLabel=int32(Xi);  
        str1=['��',num2str(i),'�׶�·��Сʱpcu'];
        title(str1)
        P.XData=XData0';
        P.YData=YData0';
        colorbar
         colormap spring;  %colormap����������Ϊjet��hsv��hot��spring��summer��autumn��winter��gray��bone��copper��pink��lines��
        
end
end
function G1=G_append_OD(G0)
G=G0;
G_Nodes=G.Nodes;
G_Nodes_NODEID=G_Nodes.NODEID;
G_Nodes_YoN=G_Nodes.YoN;
G_Nodes_XDate=G_Nodes.XDate;
G_Nodes_YDate=G_Nodes.YDate;
j=0;
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_YoN(i)==1
        j=j+1;
        nodepossion(j)=i;
        NODEID_od(j)=(i);%���ں��湹��Edges
    end
end
N=size(nodepossion,2);
MaxNodeID=max(G_Nodes.NODEID);
for i1=1:N
    NODEID_append(i1)=MaxNodeID+i1;
    YoN_append(i1)=1;
    XDate_append(i1)=G_Nodes_XDate(nodepossion(i1))-1000;
    YDate_append(i1)=G_Nodes_YDate(nodepossion(i1));  
end
NODEID_new=[G_Nodes_NODEID;NODEID_append'];
YoN_new=[G_Nodes_YoN;YoN_append'];
XDate_new=[G_Nodes_XDate;XDate_append'];
YDate_new=[G_Nodes_YDate;YDate_append'];
G_Nodes_New=table;
G_Nodes_New.NODEID=NODEID_new;
G_Nodes_New.YoN=YoN_new;
G_Nodes_New.XDate=XDate_new;
G_Nodes_New.YDate=YDate_new;

%% Edges
G_Edges=G.Edges;
G_Edges_EndNodes=G_Edges.EndNodes;
G_Edges_Cap_d1=G_Edges.Cap_d1;
G_Edges_Cap_d2=G_Edges.Cap_d2;
G_Edges_Speed_d=G_Edges.Speed_d;
G_Edges_Distance=G_Edges.Distance;
G_Edges_PCI1=G_Edges.PCI1;
G_Edges_PCI2=G_Edges.PCI2;
EndNodes_append=[NODEID_od',NODEID_append';NODEID_append',NODEID_od'];
Edges_adnum=size(EndNodes_append,1);
Cap_d1_append=repelem(G_Edges_Cap_d1(1),Edges_adnum);
Cap_d2_append=repelem(G_Edges_Cap_d2(1),Edges_adnum);
Speed_d_append=repelem(G_Edges_Speed_d(1)*10,Edges_adnum);
Distance_append=repelem(1000,Edges_adnum);
PCI1_append=repelem(G_Edges_PCI1(1),Edges_adnum);
PCI2_append=repelem(G_Edges_PCI2(1),Edges_adnum);
G_Edges_New=table;
G_Edges_New.EndNodes=[G_Edges_EndNodes;EndNodes_append];
G_Edges_New.Cap_d1=[G_Edges_Cap_d1;Cap_d1_append'];
G_Edges_New.Cap_d2=[G_Edges_Cap_d2;Cap_d2_append'];
G_Edges_New.Speed_d=[G_Edges_Speed_d;Speed_d_append'];
G_Edges_New.Distance=[G_Edges_Distance;Distance_append'];
G_Edges_New.PCI1=[G_Edges_PCI1;PCI1_append'];
G_Edges_New.PCI2=[G_Edges_PCI2;PCI2_append'];
%����ͼ
G2=digraph(G_Edges,G_Nodes);
G1=digraph(G_Edges_New,G_Nodes_New);
end%����OD��ʾ

%����������α仯ͼ
function EdgesPlot(G_FlowD,plot_Edges)
plot_vexnum=size(plot_Edges,1);%����
s=fix(plot_vexnum/10);%����̫�࣬10�η�һ��ͼ
n=length(G_FlowD);
possionnumber=zeros(plot_vexnum,n);
for jj=1:s+1
    if jj==s+1
        plot_row=mod(plot_vexnum,10);
    else
        plot_row=10;
    end
    figure;
    for i=1:plot_row
        for j=1:n
            G=G_FlowD{j};
            G_vexnum=length(G.Edges.EndNodes);%����
            %��ѯλ��
            for ii=1:G_vexnum
                if G.Edges.EndNodes(ii,1)==plot_Edges(i,1)&&G.Edges.EndNodes(ii,2)==plot_Edges(i,2)
                    possionnumber(i,j)=ii;
                end
            end
        subplot(plot_row,n,(i-1)*n+j)
        plot(G.Edges.X1i(possionnumber(i,j),:))
        end
    end
end
end

%% ��ʼfunction
function [EdgeTable,NodeTable] = topotable()
%
%   
%Ԥ����line��node����γ�topo��Ҫ���ڽӹ�ϵ���EdgeTable��NodeTable 
line=readtable('line.xls');
% t0 = line.Shape_Length./line.Speed_d*(3600/1000);
% Cap=line.Cap_d.*cap_reduction(Si);
node=readtable('node.xls');
EdgeTable_a =  table([line.FNODE_,line.TNODE_],line.Cap_d1,line.Cap_d2,line.Speed_d,line.Shape_Length,line.PCI1,line.PCI2,'VariableNames',{'EndNodes','Cap_d1','Cap_d2','Speed_d','Distance','PCI1','PCI2'});
EdgeTable_b =  table([line.TNODE_,line.FNODE_],line.Cap_d1,line.Cap_d2,line.Speed_d,line.Shape_Length,line.PCI1,line.PCI2,'VariableNames',{'EndNodes','Cap_d1','Cap_d2','Speed_d','Distance','PCI1','PCI2'});
EdgeTable=[EdgeTable_a ;EdgeTable_b];%����ͼ���ݱ�Ϊ˫��ͼ����
NodeTable =  table(node.FID_node,node.ID,node.XDate,node.YDate,'VariableNames',{'NODEID','YoN','XDate','YDate'});

end
function[OD_input_p,OD_input_f, OD_struct] = OD_matr(year,mon,NodeTable,p,f)
%% ͨ����ȡod�㣬�γ�od_pair��񣬲���ʾ���od��������
k = find(NodeTable.YoN==1);
xlswrite('OD.xlsx',k','OD_input_f','B1');
xlswrite('OD.xlsx',k,'OD_input_f','	A2');
xlswrite('OD.xlsx',k','OD_input_p','B1');
xlswrite('OD.xlsx',k,'OD_input_p','	A2');
% disp('�밴Ҫ�����OD.xlsx��·�α��OD');
% disp('�����������...');
% pause;
% OD_input_p=xlsread('OD.xlsx','OD_input_p','B2:ZZ10000');
% OD_input_f=xlsread('OD.xlsx','OD_input_f','B2:ZZ10000');

%% Ϊ�����ʵ�飬���Զ��������od����ʵ����ʷ�������룬������od_pair���
t1=8;t2=16;%t1,t2�ֱ�Ϊ�߷���/ƽ���ڵĳ���ʱ����/h��
OD_cell=cell(2,12);
A=cell(1,12);
B=cell(1,12);
OD_num=length(k);
OD_Node=zeros(OD_num,1);
for i = 1:OD_num
    OD_Node(i,1)=NodeTable.NODEID(k(i));
end
mon_i=[1,1,1.5,1.5,1,1,1,2,2,1,5,1,1];%��ͬ�·ݣ�����ı仯
for i=1:12
    A{1,i}=round(rand(OD_num)+mon_i(i)*fix(mod(p,100)/10)*10+fix(p/100)*100);%�������OD����io2OD().%��ȡʮλ���Ͱ�λ��fix(mod(p,100)/10)+fix(p/100)*100)
    OD_cell{1,i}=A{1,i}-diag(diag(A{1,i}));%�������OD����io2OD()ƽ����(16)
    B{1,i}=round(rand(OD_num)+mon_i(i)*fix(mod(f,100)/10)*10+fix(f/100)*100);%�������OD����io2OD()
    OD_cell{2,i}=B{1,i}-diag(diag(B{1,i}));%�������OD����io2OD()�߷���(8)
end
OD_name=['plat';'fast'];
T=[t1;t2];
OD_struct=struct('name',[],'volume',OD_cell,'t',[]);
for i=1:12
    OD_struct(1,i)=struct('name',OD_name(1,1),'volume',OD_cell{1,i},'t',T(1,1));
    OD_struct(2,i)=struct('name',OD_name(2,1),'volume',OD_cell{2,i},'t',T(2,1));
end
OD_input_p=OD_struct(1,mon).volume;
OD_input_f=OD_struct(2,mon).volume;
xlswrite('OD.xlsx',OD_input_f,'OD_input_f','B2');
xlswrite('OD.xlsx',OD_input_p,'OD_input_p','B2');
end

%% ·������function
%% G_Topu function
%ͬһʱ�β�ͬ·��·������
function G_update=UpdateTopu_meantime(time_event,parameter,G0)
n=size(time_event,1);
G_temp=cell(n+1,1);
G_temp{1}=G0;
for i=1:n
     parameter_temp=str2num(parameter{i});
     if time_event(i,3)==1
         if time_event(i,2)==1
             G_temp{i+1}=GraphChange_A(G_temp{i},parameter_temp);%[18,19,300,500]����A��·�ξֲ���գ�����ָ����򿪷���������;
         elseif time_event(i,2)==2
             G_temp{i+1}=GraphChange_B(G_temp{i},parameter_temp);%[27,18;34,29],28����B���ѵ����%����A��·�ξֲ���գ�����ָ����򿪷���������;
         end
         elseif time_event(i,3)==0%����AB��Ҫ�޸�Ϊ��Ӧ�ķ�����
             if time_event(i,2)==1
                 G_temp{i+1}=GraphChange_A(G_temp{i},parameter_temp);%[18,19,300,500]����A��·�ξֲ���գ�����ָ����򿪷���������;
             elseif time_event(i,2)==2
                 G_temp{i+1}=GraphChange_B(G_temp{i},parameter_temp);%[27,18;34,29],28����B���ѵ����%����A��·�ξֲ���գ�����ָ����򿪷���������;
             end
     end
end
G_update=G_temp{n+1};
end
%һ��ʱ��ֻ��һ�����������,·�����˸���
%����A
function G_save=GraphChange_A(G,parameter_temp)%���main.k����Ķ�+��������ƽ��100������ʾ
% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
G_origin=G;%ԭʼ��G
node1=parameter_temp(1);node2=parameter_temp(2);fd=parameter_temp(3);length=parameter_temp(4);
%% ȡ��ԭʼG�еĸ�������
G_Edges=G.Edges;
G_Nodes=G.Nodes;
G_Edges_EndNodes=G_Edges.EndNodes;
% G_Edges_Weight=G_Edges.Weight;
G_Edges_Cap_d1=G_Edges.Cap_d1;
G_Edges_Cap_d2=G_Edges.Cap_d2;
G_Edges_Speed_d=G_Edges.Speed_d;
G_Edges_Distance=G_Edges.Distance;
G_Edges_PCI1=G_Edges.PCI1;
G_Edges_PCI2=G_Edges.PCI2;
% G_Edges_Num_lane=G_Edges.Num_lane;
G_Nodes_NODEID=G_Nodes.NODEID;
G_Nodes_YoN=G_Nodes.YoN;
G_Nodes_XDate=G_Nodes.XDate;
G_Nodes_YDate=G_Nodes.YDate;
%% ���ڵ������Ա������½ڵ�
MaxNodeID=max(G_Nodes.NODEID);
%% �ҳ���Ҫ�޸Ľڵ��λ��
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node1&&G_Edges_EndNodes(i,2)==node2
        PosisionNumber=i;
    end
end
%���򳵵���Ӧedgeλ��
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node2&&G_Edges_EndNodes(i,2)==node1
        PosisionNumber_r=i;
    end
end
%% ��Ҫ�޸ĵĲ��ֵ�ǰ��ͺ��涼ȡ����
if PosisionNumber==1
    G_Edges_EndNodes_Before=[];
%     G_Edges_Weight_Before=[];
    G_Edges_Cap_d1_Before=[];
    G_Edges_Cap_d2_Before=[];
    G_Edges_Speed_d_Before=[];
    G_Edges_Distance_Before=[];
    G_Edges_PCI1_Before=[];
    G_Edges_PCI2_Before=[];
%     G_Edges_Num_lane=[];
    G_Edges_EndNodes_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Weight_After=G_Edges_Weight(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d1_After=G_Edges_Cap_d1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d2_After=G_Edges_Cap_d2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Speed_d_After=G_Edges_Speed_d(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Distance_After=G_Edges_Distance(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI1_After=G_Edges_PCI1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI2_After=G_Edges_PCI2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Num_lane_After=G_Edges_Num_lane(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
elseif PosisionNumber==size(G_Edges_EndNodes,1)
    G_Edges_EndNodes_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
%     G_Edges_Weight_Before=G_Edges_Weight(1:PosisionNumber-1,:);
    G_Edges_Cap_d1_Before=G_Edges_Cap_d1(1:PosisionNumber-1,:);
    G_Edges_Cap_d2_Before=G_Edges_Cap_d2(1:PosisionNumber-1,:);
    G_Edges_Speed_d_Before=G_Edges_Speed_d(1:PosisionNumber-1,:);
    G_Edges_Distance_Before=G_Edges_Distance(1:PosisionNumber-1,:);
    G_Edges_PCI1_Before=G_Edges_PCI1(1:PosisionNumber-1,:);
    G_Edges_PCI2_Before=G_Edges_PCI2(1:PosisionNumber-1,:);
%     G_Edges_Num_lane_Before=G_Edges_Num_lane(1:PosisionNumber-1,:);
    G_Edges_EndNodes_After=[];
%     G_Edges_Weight_After=[];
    G_Edges_Cap_d1_After=[];
    G_Edges_Cap_d2_After=[];
    G_Edges_Speed_d_After=[];
    G_Edges_Distance_After=[];
    G_Edges_PCI1_After=[];
    G_Edges_PCI2_After=[];
%     G_Edges_Num_lane_After=[];
else
    G_Edges_EndNodes_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
%     G_Edges_Weight_Before=G_Edges_Weight(1:PosisionNumber-1,:);
    G_Edges_Cap_d1_Before=G_Edges_Cap_d1(1:PosisionNumber-1,:);
    G_Edges_Cap_d2_Before=G_Edges_Cap_d2(1:PosisionNumber-1,:);
    G_Edges_Speed_d_Before=G_Edges_Speed_d(1:PosisionNumber-1,:);
    G_Edges_Distance_Before=G_Edges_Distance(1:PosisionNumber-1,:);
    G_Edges_PCI1_Before=G_Edges_PCI1(1:PosisionNumber-1,:);
    G_Edges_PCI2_Before=G_Edges_PCI2(1:PosisionNumber-1,:);
%     G_Edges_Num_lane_Before=G_Edges_Num_lane(1:PosisionNumber-1,:);
    G_Edges_EndNodes_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Weight_After=G_Edges_Weight(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d1_After=G_Edges_Cap_d1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d2_After=G_Edges_Cap_d2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Speed_d_After=G_Edges_Speed_d(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Distance_After=G_Edges_Distance(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI1_After=G_Edges_PCI1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI2_After=G_Edges_PCI2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Num_lane_After=G_Edges_Num_lane(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
end
%% ��������ͼ��Nodes������������½ڵ㣬MaxNodeID+1,MaxNodeID+2
%Ŀ��·�νڵ�����
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node1
        nodepossion1=i;
    end
end
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node2
        nodepossion2=i;
    end
end
XDate_node1=G_Nodes_XDate(nodepossion1);
XDate_node2=G_Nodes_XDate(nodepossion2);
YDate_node1=G_Nodes_YDate(nodepossion1);
YDate_node2=G_Nodes_YDate(nodepossion2);
%�����ڵ�����
a=[XDate_node1,YDate_node1];
b=[XDate_node2,YDate_node2];
L_node12=norm(b-a);
XDate_node3= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
YDate_node3= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
XDate_node4= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
YDate_node4= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
%��������ƽ��
dertaX=100;
dertaY=dertaX*XDate_node1/YDate_node2;
%������ͼ
G_Nodes_New=table;
G_Nodes_NODEID_New=[G_Nodes_NODEID;MaxNodeID+1;MaxNodeID+2];
G_Nodes_YoN_New=[G_Nodes_YoN;0;0];%������YoN�趨Ϊ0
G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node3+dertaX;XDate_node4+dertaX];%����XDate
G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node3+dertaY;YDate_node4+dertaY];%����YDate
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% ��������ͼ��Edges���Բ����½ڵ㣬ʹ��һ�α������
G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),MaxNodeID+1;
    MaxNodeID+1,MaxNodeID+2;
    MaxNodeID+2,G_Edges_EndNodes(PosisionNumber,2)];
% G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber)];
G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
    G_Edges_Cap_d2(PosisionNumber_r);
    G_Edges_Cap_d1(PosisionNumber)];%�����Ƿ���Ҫ�����龳�ۼ�
G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
    0;
    G_Edges_Cap_d2(PosisionNumber)];%�����Ƿ���Ҫ�����龳�ۼ�
G_Edges_Speed_d_Append=[G_Edges_Speed_d(PosisionNumber);
    G_Edges_Speed_d(PosisionNumber);
    G_Edges_Speed_d(PosisionNumber)];
G_Edges_Distance_Append=[fd;
    G_Edges_Distance(PosisionNumber)-fd-length;
    length];
G_Edges_PCI1_Append=[G_Edges_PCI1(PosisionNumber);
    G_Edges_PCI2(PosisionNumber_r);
    G_Edges_PCI1(PosisionNumber)];
G_Edges_PCI2_Append=[G_Edges_PCI2(PosisionNumber);
    0;
    G_Edges_PCI2(PosisionNumber)];
% G_Edges_Num_lane_Append=[G_Edges_Num_lane(PosisionNumber);
%     G_Edges_Num_lane(PosisionNumber);
%     G_Edges_Num_lane(PosisionNumber)]*0.5;%��������
%% �������ӵĽڵ���֮ǰû�仯��ƴ������
G_Edges_EndNodes_New=[G_Edges_EndNodes_Before;
    G_Edges_EndNodes_Append;
    G_Edges_EndNodes_After];
% G_Edges_Weight_New=[G_Edges_Weight_Before;
%     G_Edges_Weight_Append;
%     G_Edges_Weight_After];
G_Edges_Cap_d1_New=[G_Edges_Cap_d1_Before;
    G_Edges_Cap_d1_Append;
    G_Edges_Cap_d1_After];
G_Edges_Cap_d2_New=[G_Edges_Cap_d2_Before;
    G_Edges_Cap_d2_Append;
    G_Edges_Cap_d2_After];
G_Edges_Speed_d__New=[G_Edges_Speed_d_Before;
    G_Edges_Speed_d_Append;
    G_Edges_Speed_d_After];
G_Edges_Distance_New=[G_Edges_Distance_Before;
    G_Edges_Distance_Append;
    G_Edges_Distance_After];
G_Edges_PCI1_New=[G_Edges_PCI1_Before;
    G_Edges_PCI1_Append;
    G_Edges_PCI1_After];
G_Edges_PCI2_New=[G_Edges_PCI2_Before;
    G_Edges_PCI2_Append;
    G_Edges_PCI2_After];
% G_Edges_Num_lane_New=[G_Edges_Num_lane_Before;
%     G_Edges_Num_lane_Append;
%     G_Edges_Num_lane_After];
%% �������Ը��赽������ͼ��Edges������
G_Edges_New=table;
G_Edges_New.EndNodes=G_Edges_EndNodes_New;
% G_Edges_New.Weight=G_Edges_Weight_New;
G_Edges_New.Cap_d1=G_Edges_Cap_d1_New;
G_Edges_New.Cap_d2=G_Edges_Cap_d2_New;
G_Edges_New.Speed_d=G_Edges_Speed_d__New;
G_Edges_New.Distance=G_Edges_Distance_New;
G_Edges_New.PCI1=G_Edges_PCI1_New;
G_Edges_New.PCI2=G_Edges_PCI2_New;
% G_Edges_New.Num_lane=G_Edges_Num_lane_New;
%% �����µ�����ͼ
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% ��ȡ�����Ľڵ�
add_node1=max(G_new.Nodes.NODEID)-1;
r_parameter=[parameter_temp,add_node1];
%% ���¾�����ͼ�ŵ�Ԫ��������
% G1{1}=G_origin;
G1{1}=G_new;
G1{2}= r_parameter;
G_save=G1;

end
function G_save=GraphChange_AA(G,parameter_temp)%���A(����A)�������˶���ֶΣ����Ҫʹ�ã�Ӧ�ò���G_LaneProperty_renew_A�޸�G_LaneProperty_Change_A����Ϊ19���Ѿ���A_r�б�ɾ��
G_save_temp=GraphChange_A(G,parameter_temp);
G_new=GraphChange_A_r(G_save_temp{1},parameter_temp);
G_save={G_new,G_save_temp{2}};
end
function G_new=GraphChange_A_r(G,parameter_temp)%���main.k����Ķ�%��GraphChange_A���������Ӷ���19��18�ֶ�
% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
G_origin=G;%ԭʼ��G
node1=parameter_temp(2);node2=parameter_temp(1);fd=parameter_temp(4);length=parameter_temp(3);
%% ȡ��ԭʼG�еĸ�������
G_Edges=G.Edges;
G_Nodes=G.Nodes;
G_Edges_EndNodes=G_Edges.EndNodes;
% G_Edges_Weight=G_Edges.Weight;
G_Edges_Cap_d1=G_Edges.Cap_d1;
G_Edges_Cap_d2=G_Edges.Cap_d2;
G_Edges_Speed_d=G_Edges.Speed_d;
G_Edges_Distance=G_Edges.Distance;
G_Edges_PCI1=G_Edges.PCI1;
G_Edges_PCI2=G_Edges.PCI2;
% G_Edges_X1i=G_Edges.X1i;
% G_Edges_X2i=G_Edges.X2i;
% G_Edges_Num_lane=G_Edges.Num_lane;
G_Nodes_NODEID=G_Nodes.NODEID;
G_Nodes_YoN=G_Nodes.YoN;
G_Nodes_XDate=G_Nodes.XDate;
G_Nodes_YDate=G_Nodes.YDate;
%% ���ڵ������Ա������½ڵ�
MaxNodeID=max(G_Nodes.NODEID);
%% �ҳ���Ҫ�޸Ľڵ��λ��
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node1&&G_Edges_EndNodes(i,2)==node2
        PosisionNumber=i;
    end
end
%% ��Ҫ�޸ĵĲ��ֵ�ǰ��ͺ��涼ȡ����
if PosisionNumber==1
    G_Edges_EndNodes_Before=[];
%     G_Edges_Weight_Before=[];
    G_Edges_Cap_d1_Before=[];
    G_Edges_Cap_d2_Before=[];
    G_Edges_Speed_d_Before=[];
    G_Edges_Distance_Before=[];
    G_Edges_PCI1_Before=[];
    G_Edges_PCI2_Before=[];
%    G_Edges_X1i_Before=[];
%     G_Edges_X2i_Before=[];
%     G_Edges_Num_lane=[];
    G_Edges_EndNodes_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Weight_After=G_Edges_Weight(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d1_After=G_Edges_Cap_d1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d2_After=G_Edges_Cap_d2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Speed_d_After=G_Edges_Speed_d(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Distance_After=G_Edges_Distance(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI1_After=G_Edges_PCI1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI2_After=G_Edges_PCI2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_X1i_After=G_Edges_X1i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_X2i_After=G_Edges_X1i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Num_lane_After=G_Edges_Num_lane(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
elseif PosisionNumber==size(G_Edges_EndNodes,1)
    G_Edges_EndNodes_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
%     G_Edges_Weight_Before=G_Edges_Weight(1:PosisionNumber-1,:);
    G_Edges_Cap_d1_Before=G_Edges_Cap_d1(1:PosisionNumber-1,:);
    G_Edges_Cap_d2_Before=G_Edges_Cap_d2(1:PosisionNumber-1,:);
    G_Edges_Speed_d_Before=G_Edges_Speed_d(1:PosisionNumber-1,:);
    G_Edges_Distance_Before=G_Edges_Distance(1:PosisionNumber-1,:);
    G_Edges_PCI1_Before=G_Edges_PCI1(1:PosisionNumber-1,:);
    G_Edges_PCI2_Before=G_Edges_PCI2(1:PosisionNumber-1,:);
%     G_Edges_X1i_Before=G_Edges_X1i(1:PosisionNumber-1,:);
%     G_Edges_X2i_Before=G_Edges_X1i(1:PosisionNumber-1,:);
%     G_Edges_Num_lane_Before=G_Edges_Num_lane(1:PosisionNumber-1,:);
    G_Edges_EndNodes_After=[];
%     G_Edges_Weight_After=[];
    G_Edges_Cap_d1_After=[];
    G_Edges_Cap_d2_After=[];
    G_Edges_Speed_d_After=[];
    G_Edges_Distance_After=[];
    G_Edges_PCI1_After=[];
    G_Edges_PCI2_After=[];
%     G_Edges_X1i_After=[];
%     G_Edges_X2i_After=[];
%     G_Edges_Num_lane_After=[];
else
    G_Edges_EndNodes_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
%     G_Edges_Weight_Before=G_Edges_Weight(1:PosisionNumber-1,:);
    G_Edges_Cap_d1_Before=G_Edges_Cap_d1(1:PosisionNumber-1,:);
    G_Edges_Cap_d2_Before=G_Edges_Cap_d2(1:PosisionNumber-1,:);
    G_Edges_Speed_d_Before=G_Edges_Speed_d(1:PosisionNumber-1,:);
    G_Edges_Distance_Before=G_Edges_Distance(1:PosisionNumber-1,:);
    G_Edges_PCI1_Before=G_Edges_PCI1(1:PosisionNumber-1,:);
    G_Edges_PCI2_Before=G_Edges_PCI2(1:PosisionNumber-1,:);
%     G_Edges_X1i_Before=G_Edges_X1i(1:PosisionNumber-1,:);
%     G_Edges_X2i_Before=G_Edges_X2i(1:PosisionNumber-1,:);
%     G_Edges_Num_lane_Before=G_Edges_Num_lane(1:PosisionNumber-1,:);
    G_Edges_EndNodes_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Weight_After=G_Edges_Weight(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d1_After=G_Edges_Cap_d1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d2_After=G_Edges_Cap_d2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Speed_d_After=G_Edges_Speed_d(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Distance_After=G_Edges_Distance(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI1_After=G_Edges_PCI1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI2_After=G_Edges_PCI2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_X1i_After=G_Edges_X1i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_X2i_After=G_Edges_X2i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Num_lane_After=G_Edges_Num_lane(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
end
%% ��������ͼ��Nodes������������½ڵ㣬MaxNodeID+1,MaxNodeID+2
%Ŀ��·�νڵ�����
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node1
        nodepossion1=i;
    end
end
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node2
        nodepossion2=i;
    end
end
XDate_node1=G_Nodes_XDate(nodepossion1);
XDate_node2=G_Nodes_XDate(nodepossion2);
YDate_node1=G_Nodes_YDate(nodepossion1);
YDate_node2=G_Nodes_YDate(nodepossion2);
%�����ڵ�����
a=[XDate_node1,YDate_node1];
b=[XDate_node2,YDate_node2];
L_node12=norm(b-a);
XDate_node3= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
YDate_node3= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
XDate_node4= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
YDate_node4= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
%������ͼ
G_Nodes_New=table;
G_Nodes_NODEID_New=[G_Nodes_NODEID;MaxNodeID+1;MaxNodeID+2];
G_Nodes_YoN_New=[G_Nodes_YoN;0;0];%������YoN�趨Ϊ0
G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node3;XDate_node4];%����XDate
G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node3;YDate_node4];%����YDate
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% ��������ͼ��Edges���Բ����½ڵ㣬ʹ��һ�α������
G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),MaxNodeID+1;
    MaxNodeID+1,MaxNodeID+2;
    MaxNodeID+2,G_Edges_EndNodes(PosisionNumber,2)];
% G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber)];
G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
    G_Edges_Cap_d2(PosisionNumber);
    G_Edges_Cap_d1(PosisionNumber)];%�����Ƿ���Ҫ�����龳�ۼ�
G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
    0;
    G_Edges_Cap_d2(PosisionNumber)];%�����Ƿ���Ҫ�����龳�ۼ�
G_Edges_Speed_d_Append=[G_Edges_Speed_d(PosisionNumber);
    G_Edges_Speed_d(PosisionNumber);
    G_Edges_Speed_d(PosisionNumber)];
G_Edges_Distance_Append=[fd;
    G_Edges_Distance(PosisionNumber)-fd-length;
    length];
G_Edges_PCI1_Append=[G_Edges_PCI1(PosisionNumber);
    G_Edges_PCI2(PosisionNumber);
    G_Edges_PCI1(PosisionNumber)];
G_Edges_PCI2_Append=[G_Edges_PCI2(PosisionNumber);
    0;
    G_Edges_PCI2(PosisionNumber)];
% G_Edges_X1i_Append=[G_Edges_X1i(PosisionNumber);
%     G_Edges_X1i(PosisionNumber)+G_Edges_X2i(PosisionNumber);
%     G_Edges_X1i(PosisionNumber)];
% G_Edges_X2i_Append=[G_Edges_X2i(PosisionNumber);
%     0;
%     G_Edges_X2i(PosisionNumber)];
% G_Edges_Num_lane_Append=[G_Edges_Num_lane(PosisionNumber);
%     G_Edges_Num_lane(PosisionNumber);
%     G_Edges_Num_lane(PosisionNumber)]*0.5;%��������
%% �������ӵĽڵ���֮ǰû�仯��ƴ������
G_Edges_EndNodes_New=[G_Edges_EndNodes_Before;
    G_Edges_EndNodes_Append;
    G_Edges_EndNodes_After];
% G_Edges_Weight_New=[G_Edges_Weight_Before;
%     G_Edges_Weight_Append;
%     G_Edges_Weight_After];
G_Edges_Cap_d1_New=[G_Edges_Cap_d1_Before;
    G_Edges_Cap_d1_Append;
    G_Edges_Cap_d1_After];
G_Edges_Cap_d2_New=[G_Edges_Cap_d2_Before;
    G_Edges_Cap_d2_Append;
    G_Edges_Cap_d2_After];
G_Edges_Speed_d__New=[G_Edges_Speed_d_Before;
    G_Edges_Speed_d_Append;
    G_Edges_Speed_d_After];
G_Edges_Distance_New=[G_Edges_Distance_Before;
    G_Edges_Distance_Append;
    G_Edges_Distance_After];
G_Edges_PCI1_New=[G_Edges_PCI1_Before;
    G_Edges_PCI1_Append;
    G_Edges_PCI1_After];
G_Edges_PCI2_New=[G_Edges_PCI2_Before;
    G_Edges_PCI2_Append;
    G_Edges_PCI2_After];
% G_Edges_X1i_New=[G_Edges_X1i_Before;
%     G_Edges_X1i_Append;
%     G_Edges_X1i_After];
% G_Edges_X2i_New=[G_Edges_X2i_Before;
%     G_Edges_X2i_Append;
%     G_Edges_X2i_After];
% G_Edges_Num_lane_New=[G_Edges_Num_lane_Before;
%     G_Edges_Num_lane_Append;
%     G_Edges_Num_lane_After];
%% �������Ը��赽������ͼ��Edges������
G_Edges_New=table;
G_Edges_New.EndNodes=G_Edges_EndNodes_New;
% G_Edges_New.Weight=G_Edges_Weight_New;
G_Edges_New.Cap_d1=G_Edges_Cap_d1_New;
G_Edges_New.Cap_d2=G_Edges_Cap_d2_New;
G_Edges_New.Speed_d=G_Edges_Speed_d__New;
G_Edges_New.Distance=G_Edges_Distance_New;
G_Edges_New.PCI1=G_Edges_PCI1_New;
G_Edges_New.PCI2=G_Edges_PCI2_New;
% G_Edges_New.X1i=G_Edges_X1i_New;
% G_Edges_New.X2i=G_Edges_X2i_New;
% G_Edges_New.Num_lane=G_Edges_Num_lane_New;
%% �����µ�����ͼ
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% ��ȡ�����Ľڵ�
add_node1=max(G_new.Nodes.NODEID)-1;
r_parameter=[parameter_temp,add_node1];
%% ���¾�����ͼ�ŵ�Ԫ��������
% G1{1}=G_origin;
G1{1}=G_new;
G1{2}= r_parameter;
G_save=G1;

end

%����B
function G_save=GraphChange_B(G,node1_3,node2)
%node1��node2��node3�Ƕ�Ӧ·�����з����Եģ�node1��node2��node3����Ϊǰ������
%node1_3����һ����node1��Ӧ�������ڶ�����node3��Ӧ����������node1_3=[2,3;15,4]����Ϊ���2��3,15��4����
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
G_origin=G;%ԭʼ��G
%% ȡ��ԭʼG�еĸ�������
G_Edges=G.Edges;
G_Nodes=G.Nodes;
G_Edges_EndNodes=G_Edges.EndNodes;
% G_Edges_Weight=G_Edges.Weight;
G_Edges_Cap_d1=G_Edges.Cap_d1;
G_Edges_Cap_d2=G_Edges.Cap_d2;
G_Edges_Speed_d=G_Edges.Speed_d;
G_Edges_Distance=G_Edges.Distance;
G_Edges_PCI1=G_Edges.PCI1;
G_Edges_PCI2=G_Edges.PCI2;
G_Edges_X1i=G_Edges.X1i;
G_Edges_X2i=G_Edges.X2i;
% G_Edges_Num_lane=G_Edges.Num_lane;
G_Nodes_NODEID=G_Nodes.NODEID;
G_Nodes_YoN=G_Nodes.YoN;
G_Nodes_XDate=G_Nodes.XDate;
G_Nodes_YDate=G_Nodes.YDate;
%% ������node2�����Ľڵ��ж��ٸ�,����Ҫ���ӵĽڵ�����node_num
cnt_node=[];
cnt_edges=[];
for i=1:size(G_Edges_EndNodes,1)
    node_nargin=G_Edges_EndNodes(i,:)-node2;
    if node_nargin(1)*node_nargin(2)==0%������һ��Ϊ0��˵����һ������node2
        cnt_node=[cnt_node;node_nargin(1)+node_nargin(2)];
        cnt_edges=[cnt_edges;G_Edges_EndNodes(i,:)];
    end
end
% node_num=size(unique(cnt_node),1);%��ͬԪ�صĸ���
node_num=size(cnt_node,1);%��ͬԪ�صĸ���
%% ���ڵ������Ա������½ڵ�
MaxNodeID=max(G_Nodes.NODEID);
%% �½�node_num���ڵ�
G_Nodes_New=table;
G_Nodes_NODEID_New=G_Nodes_NODEID;
for i=1:length(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node2
        k=i;
    end 
end
XDate=G_Nodes_XDate(k);
YDate=G_Nodes_YDate(k);
for i=1:node_num
    G_Nodes_NODEID_New=[G_Nodes_NODEID_New;MaxNodeID+i];
end
G_Nodes_YoN_New=G_Nodes_YoN;
G_Nodes_XDate_New=G_Nodes_XDate;
G_Nodes_YDate_New=G_Nodes_YDate;
for i=1:node_num
    G_Nodes_YoN_New=[G_Nodes_YoN_New;0];%������YoN�趨Ϊ0
    G_Nodes_XDate_New=[G_Nodes_XDate_New;XDate];
    G_Nodes_YDate_New=[G_Nodes_YDate_New;YDate];
end
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% ��������ͼ��Edges,�Ѱ�����node2�Ͳ�����node2��Edge�ֿ�
G_Edges_EndNodes_Del=[];
% G_Edges_Weight_Del=[];
G_Edges_Cap_d1_Del=[];
G_Edges_Cap_d2_Del=[];
G_Edges_Speed_d_Del=[];
G_Edges_Distance_Del=[];
G_Edges_PCI1_Del=[];
G_Edges_PCI2_Del=[];
G_Edges_X1i_Del=[];
G_Edges_X2i_Del=[];
% G_Edges_Num_lane_Del=[];
G_Edges_EndNodes_Cha=[];
% G_Edges_Weight_Cha=[];
G_Edges_Cap_d1_Cha=[];
G_Edges_Cap_d2_Cha=[];
G_Edges_Speed_d_Cha=[];
G_Edges_Distance_Cha=[];
G_Edges_PCI1_Cha=[];
G_Edges_PCI2_Cha=[];
G_Edges_X1i_Cha=[];
G_Edges_X2i_Cha=[];
% G_Edges_Num_lane_Cha=[];
for i=1:size(G_Edges,1)
    Del_nargin=0;
    for j=1:size(cnt_edges,1)
        if G_Edges_EndNodes(i,1)==cnt_edges(j,1)&&G_Edges_EndNodes(i,2)==cnt_edges(j,2)
            Del_nargin=1;
            G_Edges_EndNodes_Cha=[G_Edges_EndNodes_Cha;G_Edges_EndNodes(i,:)];
%             G_Edges_Weight_Cha=[G_Edges_Weight_Cha;G_Edges_Weight(i,:)];
            G_Edges_Cap_d1_Cha=[G_Edges_Cap_d1_Cha;G_Edges_Cap_d1(i,:)];
            G_Edges_Cap_d2_Cha=[G_Edges_Cap_d2_Cha;G_Edges_Cap_d2(i,:)];
            G_Edges_Speed_d_Cha=[G_Edges_Speed_d_Cha;G_Edges_Speed_d(i,:)];
            G_Edges_Distance_Cha=[G_Edges_Distance_Cha;G_Edges_Distance(i,:)];
            G_Edges_PCI1_Cha=[G_Edges_PCI1_Cha;G_Edges_PCI1(i,:)];
            G_Edges_PCI2_Cha=[G_Edges_PCI2_Cha;G_Edges_PCI2(i,:)];
            G_Edges_X1i_Cha=[G_Edges_X1i_Cha;G_Edges_X1i(i,:)];
            G_Edges_X2i_Cha=[G_Edges_X2i_Cha;G_Edges_X2i(i,:)];
%             G_Edges_Num_lane_Cha=[G_Edges_Num_lane_Cha;G_Edges_Num_lane(i,:)];
        end
    end
    if Del_nargin==0
        G_Edges_EndNodes_Del=[G_Edges_EndNodes_Del;G_Edges_EndNodes(i,:)];
%         G_Edges_Weight_Del=[G_Edges_Weight_Del;G_Edges_Weight(i,:)];
        G_Edges_Cap_d1_Del=[G_Edges_Cap_d1_Del;G_Edges_Cap_d1(i,:)];
        G_Edges_Cap_d2_Del=[G_Edges_Cap_d2_Del;G_Edges_Cap_d2(i,:)];
        G_Edges_Speed_d_Del=[G_Edges_Speed_d_Del;G_Edges_Speed_d(i,:)];
        G_Edges_Distance_Del=[G_Edges_Distance_Del;G_Edges_Distance(i,:)];
        G_Edges_PCI1_Del=[G_Edges_PCI1_Del;G_Edges_PCI1(i,:)];
        G_Edges_PCI2_Del=[G_Edges_PCI2_Del;G_Edges_PCI2(i,:)];
        G_Edges_X1i_Del=[G_Edges_X1i_Del;G_Edges_X1i(i,:)];
        G_Edges_X2i_Del=[G_Edges_X2i_Del;G_Edges_X2i(i,:)];        
%         G_Edges_Num_lane_Del=[G_Edges_Num_lane_Del;G_Edges_Num_lane(i,:)];
    end
end
%% �Ѻ�node2�����Ľڵ㻻���½ڵ㣬����֮ǰ�ĺϲ�
% node_change=sort(unique(cnt_node+node2));

G_Edges_EndNodes_Append=cnt_edges;

% for i=1:size(cnt_edges,1)
%     for j=1:node_num
%         if node_change(j)==G_Edges_EndNodes_Cha(i,1)
%             G_Edges_EndNodes_Append=[G_Edges_EndNodes_Append;G_Edges_EndNodes_Cha(i,1),MaxNodeID+j];
%         elseif node_change(j)==G_Edges_EndNodes_Cha(i,2)
%             G_Edges_EndNodes_Append=[G_Edges_EndNodes_Append;MaxNodeID+j,G_Edges_EndNodes_Cha(i,2)];
%         end
%     end
% end

for i=1:node_num
    if cnt_edges(i,1)==node2
        G_Edges_EndNodes_Append(i,1)=MaxNodeID+i;
    elseif cnt_edges(i,2)==node2
        G_Edges_EndNodes_Append(i,2)=MaxNodeID+i;
    end
end

% G_Edges_Weight_Append=G_Edges_Weight_Cha;
G_Edges_Cap_d1_Append=G_Edges_Cap_d1_Cha;
G_Edges_Cap_d2_Append=G_Edges_Cap_d2_Cha;
G_Edges_Speed_d_Append=G_Edges_Speed_d_Cha;
G_Edges_Distance_Append=G_Edges_Distance_Cha;
G_Edges_PCI1_Append=G_Edges_PCI1_Cha;
G_Edges_PCI2_Append=G_Edges_PCI2_Cha;
G_Edges_X1i_Append=G_Edges_X1i_Cha;
G_Edges_X2i_Append=G_Edges_X2i_Cha;
% G_Edges_Num_lane_Append=G_Edges_Num_lane_Cha;

% for i=1:size(node_change,1)
%     if node_change(i)==node1
%         break1=i;
%     elseif node_change(i)==node3
%         break3=i;
%     end
% end

%% ���Ӳ�ֳ���֮��ڵ���໥���ӣ���Щ��֮���Ȩ�غ;�����Ϊ0��node1��node2��node2��node3��Ӧ�Ĳ�����

G_Edges_EndNodes_Cross=[];
% G_Edges_Weight_Cross=[];
G_Edges_Cap_d1_Cross=[];
G_Edges_Cap_d2_Cross=[];
G_Edges_Speed_d_Cross=[];
G_Edges_Distance_Cross=[];
G_Edges_PCI1_Cross=[];
G_Edges_PCI2_Cross=[];
G_Edges_X1i_Cross=[];
G_Edges_X2i_Cross=[];
% G_Edges_Num_lane_Cross=[];
for i=1:node_num
    for j=1:node_num
%         if (i-j)^2>0&&(i-break1)^2+(j-break3)^2>0
        if cnt_edges(i,2)==node2&&cnt_edges(j,1)==node2&&cnt_edges(i,1)~=cnt_edges(j,2)
            cross_nargin=0;
            for k=1:size(node1_3,1)
                if cnt_edges(i,1)==node1_3(k,1)&&cnt_edges(j,2)==node1_3(k,2)
                    cross_nargin=1;
                end
            end
            if cross_nargin==0
            G_Edges_EndNodes_Cross=[G_Edges_EndNodes_Cross;MaxNodeID+i MaxNodeID+j];
%             G_Edges_Weight_Cross=[G_Edges_Weight_Cross;0];
            G_Edges_Cap_d1_Cross=[G_Edges_Cap_d1_Cross;800];
            G_Edges_Cap_d2_Cross=[G_Edges_Cap_d2_Cross;800];
            G_Edges_Speed_d_Cross=[G_Edges_Speed_d_Cross;60];
            G_Edges_Distance_Cross=[G_Edges_Distance_Cross;0.0001];
            G_Edges_PCI1_Cross=[G_Edges_PCI1_Cross;100];
            G_Edges_PCI2_Cross=[G_Edges_PCI2_Cross;100];
            G_Edges_X1i_Cross=[G_Edges_X1i_Cross;0.0001];
            G_Edges_X2i_Cross=[G_Edges_X2i_Cross;0.0001];
%             G_Edges_Num_lane_Cross=[G_Edges_Num_lane_Cross;2];
            end
        end
    end
end
%% ���������ַ�һ��
G_Edges_EndNodes_All=[G_Edges_EndNodes_Del;
    G_Edges_EndNodes_Append;
    G_Edges_EndNodes_Cross];
% G_Edges_Weight_All=[G_Edges_Weight_Del;
%     G_Edges_Weight_Append;
%     G_Edges_Weight_Cross];
G_Edges_Cap_d1_All=[G_Edges_Cap_d1_Del;
    G_Edges_Cap_d1_Append;
    G_Edges_Cap_d1_Cross];
G_Edges_Cap_d2_All=[G_Edges_Cap_d2_Del;
    G_Edges_Cap_d2_Append;
    G_Edges_Cap_d2_Cross];
G_Edges_Speed_d_All=[G_Edges_Speed_d_Del;
    G_Edges_Speed_d_Append;
    G_Edges_Speed_d_Cross];
G_Edges_Distance_All=[G_Edges_Distance_Del;
    G_Edges_Distance_Append;
    G_Edges_Distance_Cross];
G_Edges_PCI1_All=[G_Edges_PCI1_Del;
    G_Edges_PCI1_Append;
    G_Edges_PCI1_Cross];
G_Edges_PCI2_All=[G_Edges_PCI2_Del;
    G_Edges_PCI2_Append;
    G_Edges_PCI2_Cross];
G_Edges_X1i_All=[G_Edges_X1i_Del;
    G_Edges_X1i_Append;
    G_Edges_X1i_Cross];
G_Edges_X2i_All=[G_Edges_X2i_Del;
    G_Edges_X2i_Append;
    G_Edges_X2i_Cross];
% G_Edges_Num_lane_All=[G_Edges_Num_lane_Del;
%     G_Edges_Num_lane_Append;
%     G_Edges_Num_lane_Cross];
%% �������Ը��赽������ͼ��Edges������
G_Edges_New=table;
G_Edges_New.EndNodes=G_Edges_EndNodes_All;
% G_Edges_New.Weight=G_Edges_Weight_All;
G_Edges_New.Cap_d1=G_Edges_Cap_d1_All;
G_Edges_New.Cap_d2=G_Edges_Cap_d2_All;
G_Edges_New.Speed_d=G_Edges_Speed_d_All;
G_Edges_New.Distance=G_Edges_Distance_All;
G_Edges_New.PCI1=G_Edges_PCI1_All;
G_Edges_New.PCI2=G_Edges_PCI2_All;
G_Edges_New.X1i=G_Edges_X1i_All;
G_Edges_New.X2i=G_Edges_X2i_All;
% G_Edges_New.Num_lane=G_Edges_Num_lane_All;
%% �����µ�����ͼ
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% ���¾�����ͼ�ŵ�Ԫ��������
node_2plus=[];
for i=1:node_num
    node_2plus=[node_2plus;MaxNodeID+i];
end
G1{1}=G_origin;
G1{2}=G_new;
G1{3}=node2;
G1{4}=node_2plus;
G_save=G1;
end

%% G_Property function
%����A�³������������������ط���
function G_save=G_LaneProperty_Change_A(G_MSA,G,parameter_temp)%����Ķ�+��������ƽ��100
%function [G_new,r_parameter]=GraphChange_A(G,r_parameter)
% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
G_origin=G_MSA;%ԭʼ��G
node1=parameter_temp(1);node2=parameter_temp(2);fd=parameter_temp(3);length=parameter_temp(4);
%% ȡ��ԭʼG�еĸ�������
G_Edges=G_MSA.Edges;
G_Nodes=G_MSA.Nodes;
G_Edges_EndNodes=G_Edges.EndNodes;
% G_Edges_Weight=G_Edges.Weight;
G_Edges_Cap_d1=G_Edges.Cap_d1;
G_Edges_Cap_d2=G_Edges.Cap_d2;
G_Edges_Speed_d=G_Edges.Speed_d;
G_Edges_Distance=G_Edges.Distance;
G_Edges_PCI1=G_Edges.PCI1;
G_Edges_PCI2=G_Edges.PCI2;
% G_Edges_Num_lane=G_Edges.Num_lane;
G_Edges_X1i=G_Edges.X1i;
G_Edges_X2i=G_Edges.X2i;
G_Edges_vi1=G_Edges.vi1;
G_Edges_vi2=G_Edges.vi2;
% G_Edges_Ti=G_Edges.Ti;
G_Nodes_NODEID=G_Nodes.NODEID;
G_Nodes_YoN=G_Nodes.YoN;
G_Nodes_XDate=G_Nodes.XDate;
G_Nodes_YDate=G_Nodes.YDate;
%% ���ڵ������Ա������½ڵ�
MaxNodeID=max(G_Nodes.NODEID);
%% �ҳ���Ҫ�޸Ľڵ��λ��
%��ȡ������򳵵��Ľڵ�[19��18]
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node2&&G_Edges_EndNodes(i,2)==node1
        PosisionNumber=i;
    end
end
%��ȡ����Ľڵ㣬��G_MSA����������ڵ�[47��48]
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==(MaxNodeID-1)&&G_Edges_EndNodes(i,2)==MaxNodeID
        PosisionNumber_1=i;
    end
end
%% ��Ҫ�޸ĵĲ��ֵ�ǰ��ͺ��涼ȡ����
if PosisionNumber==1
    G_Edges_EndNodes_Before=[];
%     G_Edges_Weight_Before=[];
    G_Edges_Cap_d1_Before=[];
    G_Edges_Cap_d2_Before=[];
    G_Edges_Speed_d_Before=[];
    G_Edges_Distance_Before=[];
    G_Edges_PCI1_Before=[];
    G_Edges_PCI2_Before=[];
%     G_Edges_Num_lane_Before=[];
    G_Edges_X1i_Before=[];
    G_Edges_X2i_Before=[];
    G_Edges_vi1_Before=[];
    G_Edges_vi2_Before=[];
%     G_Edges_Ti_Before=[];
    G_Edges_EndNodes_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Weight_After=G_Edges_Weight(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d1_After=G_Edges_Cap_d1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d2_After=G_Edges_Cap_d2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Speed_d_After=G_Edges_Speed_d(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Distance_After=G_Edges_Distance(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI1_After=G_Edges_PCI1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI2_After=G_Edges_PCI2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Num_lane_After=G_Edges_Num_lane(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_X1i_After=G_Edges_X1i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_X2i_After=G_Edges_X2i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_vi1_After=G_Edges_vi1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_vi2_After=G_Edges_vi2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Ti_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    elseif PosisionNumber==size(G_Edges_EndNodes,1)
    G_Edges_EndNodes_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
%     G_Edges_Weight_Before=G_Edges_Weight(1:PosisionNumber-1,:);
    G_Edges_Cap_d1_Before=G_Edges_Cap_d1(1:PosisionNumber-1,:);
    G_Edges_Cap_d2_Before=G_Edges_Cap_d2(1:PosisionNumber-1,:);
    G_Edges_Speed_d_Before=G_Edges_Speed_d(1:PosisionNumber-1,:);
    G_Edges_Distance_Before=G_Edges_Distance(1:PosisionNumber-1,:);
    G_Edges_PCI1_Before=G_Edges_PCI1(1:PosisionNumber-1,:);
    G_Edges_PCI2_Before=G_Edges_PCI2(1:PosisionNumber-1,:);
%     G_Edges_Num_lane_Before=G_Edges_Num_lane(1:PosisionNumber-1,:);
    G_Edges_X1i_Before=G_Edges_X1i(1:PosisionNumber-1,:);
    G_Edges_X2i_Before=G_Edges_X2i(1:PosisionNumber-1,:);
    G_Edges_vi1_Before=G_Edges_vi1(1:PosisionNumber-1,:);
    G_Edges_vi2_Before=G_Edges_vi2(1:PosisionNumber-1,:);
%     G_Edges_Ti_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
    G_Edges_EndNodes_After=[];
%     G_Edges_Weight_After=[];
    G_Edges_Cap_d1_After=[];
    G_Edges_Cap_d2_After=[];
    G_Edges_Speed_d_After=[];
    G_Edges_Distance_After=[];
    G_Edges_PCI1_After=[];
    G_Edges_PCI2_After=[];
%     G_Edges_Num_lane_After=[];
    G_Edges_X1i_After=[];
    G_Edges_X2i_After=[];
    G_Edges_vi1_After=[];
    G_Edges_vi2_After=[];
%     G_Edges_Ti_After=[];    
else
    G_Edges_EndNodes_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);
%     G_Edges_Weight_Before=G_Edges_Weight(1:PosisionNumber-1,:);
    G_Edges_Cap_d1_Before=G_Edges_Cap_d1(1:PosisionNumber-1,:);
    G_Edges_Cap_d2_Before=G_Edges_Cap_d2(1:PosisionNumber-1,:);
    G_Edges_Speed_d_Before=G_Edges_Speed_d(1:PosisionNumber-1,:);
    G_Edges_Distance_Before=G_Edges_Distance(1:PosisionNumber-1,:);
    G_Edges_PCI1_Before=G_Edges_PCI1(1:PosisionNumber-1,:);
    G_Edges_PCI2_Before=G_Edges_PCI2(1:PosisionNumber-1,:);
%     G_Edges_Num_lane_Before=G_Edges_Num_lane(1:PosisionNumber-1,:);
    G_Edges_X1i_Before=G_Edges_X1i(1:PosisionNumber-1,:);
    G_Edges_X2i_Before=G_Edges_X2i(1:PosisionNumber-1,:);
    G_Edges_vi1_Before=G_Edges_vi1(1:PosisionNumber-1,:);
    G_Edges_vi2_Before=G_Edges_vi2(1:PosisionNumber-1,:);
%     G_Edges_Ti_Before=G_Edges_EndNodes(1:PosisionNumber-1,:);    
    G_Edges_EndNodes_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Weight_After=G_Edges_Weight(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d1_After=G_Edges_Cap_d1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Cap_d2_After=G_Edges_Cap_d2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Speed_d_After=G_Edges_Speed_d(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_Distance_After=G_Edges_Distance(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI1_After=G_Edges_PCI1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_PCI2_After=G_Edges_PCI2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Num_lane_After=G_Edges_Num_lane(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_X1i_After=G_Edges_X1i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_X2i_After=G_Edges_X2i(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_vi1_After=G_Edges_vi1(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
    G_Edges_vi2_After=G_Edges_vi2(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
%     G_Edges_Ti_After=G_Edges_EndNodes(PosisionNumber+1:size(G_Edges_EndNodes,1),:);
end
%% ��������ͼ��Nodes������������½ڵ㣬MaxNodeID+1,MaxNodeID+2
%Ŀ��·�νڵ�����
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node1
        nodepossion1=i;
    end
end
for i=1:size(G_Nodes_NODEID)
    if G_Nodes_NODEID(i)==node2
        nodepossion2=i;
    end
end
XDate_node1=G_Nodes_XDate(nodepossion1);
XDate_node2=G_Nodes_XDate(nodepossion2);
YDate_node1=G_Nodes_YDate(nodepossion1);
YDate_node2=G_Nodes_YDate(nodepossion2);
%�����ڵ�����
a=[XDate_node1,YDate_node1];
b=[XDate_node2,YDate_node2];
L_node12=norm(b-a);
XDate_node3= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
YDate_node3= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
XDate_node4= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
YDate_node4= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
%��������ƽ��
dertaX=100;
dertaY=dertaX*XDate_node1/YDate_node2;
%������ͼ
G_Nodes_New=table;
G_Nodes_NODEID_New=[G_Nodes_NODEID;MaxNodeID+1;MaxNodeID+2];
G_Nodes_YoN_New=[G_Nodes_YoN;0;0];%������YoN�趨Ϊ0
G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node4-dertaX;XDate_node3-dertaX];%����XDate
G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node4-dertaY;YDate_node3-dertaY];%����YDate
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% ��������ͼ��Edges���Բ����½ڵ㣬ʹ��һ�α������
G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),MaxNodeID+1;
    MaxNodeID+1,MaxNodeID+2;
    MaxNodeID+2,G_Edges_EndNodes(PosisionNumber,2)];
% G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber)];
G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
    G_Edges_Cap_d1(PosisionNumber);
    G_Edges_Cap_d1(PosisionNumber)];%
G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
    G_Edges_Cap_d2(PosisionNumber);
    G_Edges_Cap_d2(PosisionNumber)];%
G_Edges_Speed_d_Append=[G_Edges_Speed_d(PosisionNumber);
    G_Edges_Speed_d(PosisionNumber);
    G_Edges_Speed_d(PosisionNumber)];
G_Edges_Distance_Append=[G_Edges_Distance(PosisionNumber)-fd-length;
     length;
    fd];
G_Edges_PCI1_Append=[G_Edges_PCI1(PosisionNumber);
    G_Edges_PCI1(PosisionNumber);
    G_Edges_PCI1(PosisionNumber)];
G_Edges_PCI2_Append=[G_Edges_PCI2(PosisionNumber);
    G_Edges_PCI2(PosisionNumber);
    G_Edges_PCI2(PosisionNumber)];
% G_Edges_Num_lane_Append=[G_Edges_Num_lane(PosisionNumber);
%     G_Edges_Num_lane(PosisionNumber);
%     G_Edges_Num_lane(PosisionNumber)];%����������
G_Edges_X1i_Append=[G_Edges_X1i(PosisionNumber);
    G_Edges_X1i(PosisionNumber)*2;
    G_Edges_X1i(PosisionNumber)];
G_Edges_X2i_Append=[G_Edges_X2i(PosisionNumber);
    G_Edges_X2i(PosisionNumber_1)*2;
    G_Edges_X2i(PosisionNumber)];
G_Edges_vi1_Append=[G_Edges_vi1(PosisionNumber);
    G_Edges_vi1(PosisionNumber);
    G_Edges_vi1(PosisionNumber)];
G_Edges_vi2_Append=[G_Edges_vi2(PosisionNumber);
    G_Edges_vi2(PosisionNumber_1);
    G_Edges_vi2(PosisionNumber)];
% G_Edges_Ti_Append=[G_Edges_Ti(PosisionNumber);
%     G_Edges_Ti(PosisionNumber);
%     G_Edges_Ti(PosisionNumber)];

%% �������ӵĽڵ���֮ǰû�仯��ƴ������
G_Edges_EndNodes_New=[G_Edges_EndNodes_Before;
    G_Edges_EndNodes_Append;
    G_Edges_EndNodes_After];
% G_Edges_Weight_New=[G_Edges_Weight_Before;
%     G_Edges_Weight_Append;
%     G_Edges_Weight_After];
G_Edges_Cap_d1_New=[G_Edges_Cap_d1_Before;
    G_Edges_Cap_d1_Append;
    G_Edges_Cap_d1_After];
G_Edges_Cap_d2_New=[G_Edges_Cap_d2_Before;
    G_Edges_Cap_d2_Append;
    G_Edges_Cap_d2_After];
G_Edges_Speed_d_New=[G_Edges_Speed_d_Before;
    G_Edges_Speed_d_Append;
    G_Edges_Speed_d_After];
G_Edges_Distance_New=[G_Edges_Distance_Before;
    G_Edges_Distance_Append;
    G_Edges_Distance_After];
G_Edges_PCI1_New=[G_Edges_PCI1_Before;
    G_Edges_PCI1_Append;
    G_Edges_PCI1_After];
G_Edges_PCI2_New=[G_Edges_PCI2_Before;
    G_Edges_PCI2_Append;
    G_Edges_PCI2_After];
% G_Edges_Num_lane_New=[G_Edges_Num_lane_Before;
%     G_Edges_Num_lane_Append;
%     G_Edges_Num_lane_After];
G_Edges_X1i_New=[G_Edges_X1i_Before;
    G_Edges_X1i_Append;
    G_Edges_X1i_After];
G_Edges_X2i_New=[G_Edges_X2i_Before;
    G_Edges_X2i_Append;
    G_Edges_X2i_After];
G_Edges_vi1_New=[G_Edges_vi1_Before;
    G_Edges_vi1_Append;
    G_Edges_vi1_After];
G_Edges_vi2_New=[G_Edges_vi2_Before;
    G_Edges_vi2_Append;
    G_Edges_vi2_After];
% G_Edges_Ti_New=[G_Edges_X2i_Before;
%     G_Edges_Ti_Append;
%     G_Edges_Ti_After];

%% �ֲ����� 
%��ȡ���ǰ���λ�õ�PCI��Cap[18��19]
for i=1:size(G.Edges.EndNodes,1)
    if G.Edges.EndNodes(i,1)==node1&&G.Edges.EndNodes(i,2)==node2
        PosisionNumber_3=i;
    end
end
%���������/�ٶ�����Ϊ0[47��48]
for i=1:size(G_Edges_EndNodes_New,1)
    if G_Edges_EndNodes_New(i,1)==(MaxNodeID-1)&&G_Edges_EndNodes_New(i,2)==MaxNodeID
        PosisionNumber_2=i;
    end
end
G_Edges_X1i_New(PosisionNumber_2)=0.001;
G_Edges_X2i_New(PosisionNumber_2)=0.001;
G_Edges_vi1_New(PosisionNumber_2)=0.001;
G_Edges_vi2_New(PosisionNumber_2)=0.001;
%Cap��PCI�ָ��ɷ��֮ǰ������
G_Edges_Cap_d1_New(PosisionNumber_2)=G.Edges.Cap_d1(PosisionNumber_3);
G_Edges_Cap_d2_New(PosisionNumber_2)=G.Edges.Cap_d2(PosisionNumber_3);
G_Edges_PCI1_New(PosisionNumber_2)=G.Edges.PCI1(PosisionNumber_3);
G_Edges_PCI2_New(PosisionNumber_2)=G.Edges.PCI2(PosisionNumber_3);


%% �������Ը��赽������ͼ��Edges������
G_Edges_New=table;
G_Edges_New.EndNodes=G_Edges_EndNodes_New;
% G_Edges_New.Weight=G_Edges_Weight_New;
G_Edges_New.Cap_d1=G_Edges_Cap_d1_New;
G_Edges_New.Cap_d2=G_Edges_Cap_d2_New;
G_Edges_New.Speed_d=G_Edges_Speed_d_New;
G_Edges_New.Distance=G_Edges_Distance_New;
G_Edges_New.PCI1=G_Edges_PCI1_New;
G_Edges_New.PCI2=G_Edges_PCI2_New;
% G_Edges_New.Num_lane=G_Edges_Num_lane_New;
G_Edges_New.X1i=G_Edges_X1i_New;
G_Edges_New.X2i=G_Edges_X2i_New;
G_Edges_New.vi1=G_Edges_vi1_New;
G_Edges_New.vi2=G_Edges_vi2_New;
% G_Edges_New.Ti=G_Edges_Ti_New;
%% �����µ�����ͼ
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% ���¾�����ͼ�ŵ�Ԫ��������
G1{1}=G_origin;
G1{2}=G_new;
G_save=G1{2};
end

function G_save=G_LaneProperty_renew_A(H_day2day_Property,G_origin,parameter_temp)
%��Ҫ��47��48��49��50����PCI�ͳ���CapҪ�ָ��ɽ�ͨ����֮ǰ��ֵ��18-19��19-18���������ػ�ʵ�ʳ����ϣ���������Ӱ��PCI
%function [G_new,r_parameter]=GraphChange_A(G,r_parameter)
% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
G1{1}=G_origin;%ԭʼ��G
G1{2}=H_day2day_Property;

% ȡ��ԭʼG�еĸ�������
G_Edges{1}=G1{1}.Edges;
G_Nodes{1}=G1{1}.Nodes;
G_Edges_EndNodes{1}=G_Edges{1}.EndNodes;

G_Edges{2}=G1{2}.Edges;
G_Nodes{2}=G1{2}.Nodes;
G_Edges_EndNodes{2}=G_Edges{2}.EndNodes;

%% �ҳ���Ҫ�޸Ľڵ��λ��
MaxNodeID=max(G_Nodes{2}.NODEID);
%��ȡ������򳵵��Ľڵ�[47��48]
for i=1:size(G_Edges_EndNodes{2},1)
    if G_Edges_EndNodes{2}(i,1)==MaxNodeID-3&&G_Edges_EndNodes{2}(i,2)==MaxNodeID-2
        PosisionNumber_1=i;
    end
end
%��ȡ����Ľڵ㣬����[49��50]
for i=1:size(G_Edges_EndNodes{2},1)
    if G_Edges_EndNodes{2}(i,1)==(MaxNodeID-1)&&G_Edges_EndNodes{2}(i,2)==MaxNodeID
        PosisionNumber_2=i;
    end
end

node1=parameter_temp(1);node2=parameter_temp(2);fd=parameter_temp(3);length=parameter_temp(4);
%��ȡ������򳵵��Ľڵ�[18��19]
for i=1:size(G_Edges_EndNodes{1},1)
    if G_Edges_EndNodes{1}(i,1)==node1&&G_Edges_EndNodes{1}(i,2)==node2
        PosisionNumber=i;
    end
end
%��ȡ����Ľڵ㣬����[19��18]
for i=1:size(G_Edges_EndNodes{1},1)
    if G_Edges_EndNodes{1}(i,1)==node2&&G_Edges_EndNodes{1}(i,2)==node1
        PosisionNumber_r=i;
    end
end

%% ��������ͼ{3}��Edges�������Իָ�Ϊ{1}�����ԣ�����{2}�����ԣ���Ҫ�ǰ��������ص���ʵ�峵����
%47-48����
G_Edges{2}.Cap_d1(PosisionNumber_1)=G_Edges{1}.Cap_d1(PosisionNumber);
G_Edges{2}.Cap_d2(PosisionNumber_1)=G_Edges{1}.Cap_d2(PosisionNumber);
G_Edges{2}.Speed_d(PosisionNumber_1)=G_Edges{2}.Speed_d(PosisionNumber);
G_Edges{2}.PCI1(PosisionNumber_1)=G_Edges{1}.PCI1(PosisionNumber);
G_Edges{2}.PCI2(PosisionNumber_1)=G_Edges{1}.PCI2(PosisionNumber);
G_Edges{2}.X1i(PosisionNumber_1)=0.001;
G_Edges{2}.X2i(PosisionNumber_1)=0.001;
G_Edges{2}.vi(PosisionNumber_1)=0.001;
% G_Edges{2}.Ti(PosisionNumber_1)=G_Edges{1}.Ti(PosisionNumber);%T�ڷ����ʱ�����¼���

%49-50����
% G_Edges{2}.Cap_d1(PosisionNumber_2)=G_Edges{1}.Cap_d1(PosisionNumber_r);
G_Edges{2}.Cap_d2(PosisionNumber_2)=G_Edges{1}.Cap_d2(PosisionNumber_r);
% G_Edges{3}.Speed_d(PosisionNumber_2)=G_Edges{2}.Speed_d(PosisionNumber_r);
% G_Edges{3}.PCI1(PosisionNumber_2)=G_Edges{1}.PCI1(PosisionNumber_r);
G_Edges{2}.PCI2(PosisionNumber_2)=G_Edges{1}.PCI2(PosisionNumber_r);
G_Edges{2}.X1i(PosisionNumber_2)=G_Edges{2}.X1i(PosisionNumber_2)+G_Edges{2}.X2i(PosisionNumber_2);
G_Edges{2}.X2i(PosisionNumber_2)=G_Edges{2}.X1i(PosisionNumber_1)+G_Edges{2}.X2i(PosisionNumber_1);
% G_Edges{2}.v1i(PosisionNumber_2)=G_Edges{2}.v1i(PosisionNumber_2)+G_Edges{2}.v1i(PosisionNumber_2);
% G_Edges{2}.v2i(PosisionNumber_2)=G_Edges{2}.v2i(PosisionNumber_1)+G_Edges{2}.v2i(PosisionNumber_1);
% G_Edges{2}.Ti(PosisionNumber_1)=G_Edges{1}.Ti(PosisionNumber);%T�ڷ����ʱ�����¼���

%% ���¾�����ͼ�ŵ�Ԫ��������
G_save=G1{2};
end

%����B�³������������������ط���
function G_save=G_LaneProperty_Change_B(G_MSA,G,node_2plus,node_2)
G_return_temp=G;
for i=1:size(G_return_temp.Edges.EndNodes,1)
    for j=1:size(G_MSA.Edges,1)
        if G_return_temp.Edges.EndNodes(i,1)==G_MSA.Edges.EndNodes(j,1)&&G_return_temp.Edges.EndNodes(i,2)==G_MSA.Edges.EndNodes(j,2)
%             G_return_temp.Edges.Weight(i,1)=G_after.Edges.Weight(j,1);
            G_return_temp.Edges.Cap_d1(i,1)=G_MSA.Edges.Cap_d1(j,1);
            G_return_temp.Edges.Cap_d2(i,1)=G_MSA.Edges.Cap_d2(j,1);
            G_return_temp.Edges.Speed_d(i,1)=G_MSA.Edges.Speed_d(j,1);
            G_return_temp.Edges.Distance(i,1)=G_MSA.Edges.Distance(j,1);
            G_return_temp.Edges.PCI1(i,1)=G_MSA.Edges.PCI1(j,1);
            G_return_temp.Edges.PCI2(i,1)=G_MSA.Edges.PCI2(j,1);
%             G_return_temp.Edges.Num_lane(i,1)=G_MSA.Edges.Num_lane(j,1);
            G_return_temp.Edges.X1i(i,1)=G_MSA.Edges.X1i(j,1);
            G_return_temp.Edges.X2i(i,1)=G_MSA.Edges.X2i(j,1);
            G_return_temp.Edges.vi1(i,1)=G_MSA.Edges.vi1(j,1);
            G_return_temp.Edges.vi2(i,1)=G_MSA.Edges.vi2(j,1);
        elseif G_return_temp.Edges.EndNodes(i,2)==node_2&&G_MSA.Edges.EndNodes(j,1)==G_return_temp.Edges.EndNodes(i,1)&&sum(G_MSA.Edges.EndNodes(j,2)==node_2plus)==1
%             G_return_temp.Edges.Weight(i,1)=G_after.Edges.Weight(j,1);
            G_return_temp.Edges.Cap_d1(i,1)=G_MSA.Edges.Cap_d1(j,1);
            G_return_temp.Edges.Cap_d2(i,1)=G_MSA.Edges.Cap_d2(j,1);
            G_return_temp.Edges.Speed_d(i,1)=G_MSA.Edges.Speed_d(j,1);
            G_return_temp.Edges.Distance(i,1)=G_MSA.Edges.Distance(j,1);
            G_return_temp.Edges.PCI1(i,1)=G_MSA.Edges.PCI1(j,1);
            G_return_temp.Edges.PCI2(i,1)=G_MSA.Edges.PCI2(j,1);
%             G_return_temp.Edges.Num_lane(i,1)=G_MSA.Edges.Num_lane(j,1);
            G_return_temp.Edges.X1i(i,1)=G_MSA.Edges.X1i(j,1);
            G_return_temp.Edges.X2i(i,1)=G_MSA.Edges.X2i(j,1);
            G_return_temp.Edges.vi1(i,1)=G_MSA.Edges.vi1(j,1);
            G_return_temp.Edges.vi2(i,1)=G_MSA.Edges.vi2(j,1);
        elseif G_return_temp.Edges.EndNodes(i,1)==node_2&&G_MSA.Edges.EndNodes(j,2)==G_return_temp.Edges.EndNodes(i,2)&&sum(G_MSA.Edges.EndNodes(j,1)==node_2plus)==1
%             G_return_temp.Edges.Weight(i,1)=G_after.Edges.Weight(j,1);
            G_return_temp.Edges.Cap_d1(i,1)=G_MSA.Edges.Cap_d1(j,1);
            G_return_temp.Edges.Cap_d2(i,1)=G_MSA.Edges.Cap_d2(j,1);
            G_return_temp.Edges.Speed_d(i,1)=G_MSA.Edges.Speed_d(j,1);
            G_return_temp.Edges.Distance(i,1)=G_MSA.Edges.Distance(j,1);
            G_return_temp.Edges.PCI1(i,1)=G_MSA.Edges.PCI1(j,1);
            G_return_temp.Edges.PCI2(i,1)=G_MSA.Edges.PCI2(j,1);
%             G_return_temp.Edges.Num_lane(i,1)=G_MSA.Edges.Num_lane(j,1);
            G_return_temp.Edges.X1i(i,1)=G_MSA.Edges.X1i(j,1);
            G_return_temp.Edges.X2i(i,1)=G_MSA.Edges.X2i(j,1);
            G_return_temp.Edges.vi1(i,1)=G_MSA.Edges.vi1(j,1);
            G_return_temp.Edges.vi2(i,1)=G_MSA.Edges.vi2(j,1);
        end
    end
end
G_save= G_return_temp;
end

%% G_Flow function
%������Ӹ������������ǰ·��·�λ���ͳһ�󣬰���ȡ�������������ں������(����UpdateFlow_LinkSeparator_A)
function [G_Flow,G_Property_Separator]=UpdateFlow_LinkSeparator(G_LaneProperty,paramater_plus,time_event)
%��ͬʱ��·��·�λ���ͳһ��������ȡ�������������ں������(Xi��PCI)
n=size(time_event,1);
paramater_r=cell(n-1,1);
G_Flow_Edges=cell(n,1);
G_Property_Separator=cell(n,1);
G_Flow_Nodes=cell(n,1);
G_Flow=cell(n,1);
for i=1:n
    if i==n
        G_Property_Separator{i}=G_LaneProperty{i};
    else
        paramater_temp=[];
        cn=0;
        for j=i+1:n
            if time_event(j,2)==1
                cn=cn+1;
                paramater_temp(cn,:)=paramater_plus{j};
            end
        end
        paramater_r{i}=paramater_temp;
        if isempty(paramater_r{i})==1
            G_Property_Separator{i}=G_LaneProperty{i};
        else
            G_Property_Separator{i}=UpdateFlow_LinkSeparator_A(G_LaneProperty{i},paramater_r{i});
            plot(G_Property_Separator{i})
        end
    end
    G_Flow_Edges{i}=table;
    G_Flow_Edges{i}=G_Property_Separator{i}.Edges(:,{'EndNodes','X1i','X2i','PCI1','PCI2'});%����ȡ�������������ں������
    G_Flow_Nodes{i}=G_Property_Separator{i}.Nodes;
    G_Flow{i}=digraph(G_Flow_Edges{i},G_Flow_Nodes{i});
end
end
%UpdateFlow_LinkSeparator��function���������ǰ·��·�λ���ͳһ
function G_save=UpdateFlow_LinkSeparator_A(G_LaneProperty,paramater_r)%����Ķ�
%�����������Ժ󣬶�ǰ���·�����зֶΣ�ʹ�ò�ͬʱ�ε�·������ֱ����ӵõ�·������

% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
m=size(paramater_r,1);
for ii=1:m
    %% 
    paramater_r_temp(1,:)=paramater_r(ii,:);
    paramater_r_temp(2,:)=[paramater_r(ii,2),paramater_r(ii,1),paramater_r(ii,3),paramater_r(ii,4),paramater_r(ii,5)+2];
    if ii==1
        G=G_LaneProperty;
    else
        G=G_new{ii-1};
    end
     %% 
    G_Edges=G.Edges;
    G_Nodes=G.Nodes;
    G_Edges_EndNodes=G_Edges.EndNodes;
    % G_Edges_Weight=G_Edges.Weight;
    G_Edges_Cap_d1=G_Edges.Cap_d1;
    G_Edges_Cap_d2=G_Edges.Cap_d2;
    G_Edges_Speed_d=G_Edges.Speed_d;
    G_Edges_Distance=G_Edges.Distance;
    G_Edges_PCI1=G_Edges.PCI1;
    G_Edges_PCI2=G_Edges.PCI2;
    % G_Edges_Num_lane=G_Edges.Num_lane;
    G_Edges_X1i=G_Edges.X1i;
    G_Edges_X2i=G_Edges.X2i;
    G_Edges_vi1=G_Edges.vi1;
    G_Edges_vi2=G_Edges.vi2;
    % G_Edges_Ti=G_Edges.Ti;
    G_Nodes_NODEID=G_Nodes.NODEID;
    G_Nodes_YoN=G_Nodes.YoN;
    G_Nodes_XDate=G_Nodes.XDate;  
    G_Nodes_YDate=G_Nodes.YDate;  
    %% ��������ͼ�ظ����Σ�1�α����Σ����������ڵ�
    G_Edges_Append_temp=cell(2,1);
     for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        % node1-node2����L
        for i=1:size(G_Nodes_NODEID)
            if G_Nodes_NODEID(i)==node1
                nodepossion1=i;
            end
        end
        for i=1:size(G_Nodes_NODEID)
            if G_Nodes_NODEID(i)==node2
                nodepossion2=i;
            end
        end
        XDate_node1=G_Nodes_XDate(nodepossion1);
        XDate_node2=G_Nodes_XDate(nodepossion2);
        YDate_node1=G_Nodes_YDate(nodepossion1);
        YDate_node2=G_Nodes_YDate(nodepossion2);
        a=[XDate_node1,YDate_node1];
        b=[XDate_node2,YDate_node2];
        L_node12=norm(b-a);
        %
        if jj==1
            fd=paramater_r_temp(jj,3);
        else
            fd=L_node12-paramater_r_temp(jj,3);
        end
        length=paramater_r_temp(jj,4);
        XDate_node3(jj)= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
        YDate_node3(jj)= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
        XDate_node4(jj)= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
        YDate_node4(jj)= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
        node5=paramater_r_temp(jj,5);
        for s=1:size(G_Edges_EndNodes,1) 
            if G_Edges_EndNodes(s,1)==node1&&G_Edges_EndNodes(s,2)==node2
                PosisionNumber=s;
            end
        end
        %% ��������ͼ��Edges���Բ����½ڵ㣬ʹ��һ�α������
        G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),node5;
            node5,node5+1;
            node5+1,G_Edges_EndNodes(PosisionNumber,2)];
        % G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
        %     G_Edges_Weight(PosisionNumber);
        %     G_Edges_Weight(PosisionNumber)];
        G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
            G_Edges_Cap_d1(PosisionNumber);
            G_Edges_Cap_d1(PosisionNumber)];
        G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
            G_Edges_Cap_d2(PosisionNumber);
            G_Edges_Cap_d2(PosisionNumber)];
        G_Edges_Speed_d_Append=[G_Edges_Speed_d(PosisionNumber);
            G_Edges_Speed_d(PosisionNumber);
            G_Edges_Speed_d(PosisionNumber)];
        G_Edges_Distance_Append=[fd;
            G_Edges_Distance(PosisionNumber)-fd-length;
            length];
        G_Edges_PCI1_Append=[G_Edges_PCI1(PosisionNumber);
            G_Edges_PCI1(PosisionNumber);
            G_Edges_PCI1(PosisionNumber)];
        G_Edges_PCI2_Append=[G_Edges_PCI2(PosisionNumber);
             G_Edges_PCI2(PosisionNumber);
            G_Edges_PCI2(PosisionNumber)];
        % G_Edges_Num_lane_Append=[G_Edges_Num_lane(PosisionNumber);
        %     G_Edges_Num_lane(PosisionNumber);
        %     G_Edges_Num_lane(PosisionNumber)];
        G_Edges_X1i_Append=[G_Edges_X1i(PosisionNumber);
            G_Edges_X1i(PosisionNumber);
            G_Edges_X1i(PosisionNumber)];
        G_Edges_X2i_Append=[G_Edges_X2i(PosisionNumber);
            G_Edges_X2i(PosisionNumber);
            G_Edges_X2i(PosisionNumber)];
        G_Edges_vi1_Append=[G_Edges_vi1(PosisionNumber);
            G_Edges_vi1(PosisionNumber);
            G_Edges_vi1(PosisionNumber)];
        G_Edges_vi2_Append=[G_Edges_vi2(PosisionNumber);
            G_Edges_vi2(PosisionNumber);
            G_Edges_vi2(PosisionNumber)];
        % G_Edges_Ti_Append=[G_Edges_Ti(PosisionNumber);
        %     G_Edges_Ti(PosisionNumber);
        %     G_Edges_Ti(PosisionNumber)];       
        %% �������ӵ�·���������һ���µ�table_append
        G_Edges_Append=table;
        G_Edges_Append.EndNodes=G_Edges_EndNodes_Append;
        % G_Edges_New.Weight=G_Edges_Weight_Append;
        G_Edges_Append.Cap_d1=G_Edges_Cap_d1_Append;
        G_Edges_Append.Cap_d2=G_Edges_Cap_d2_Append;
        G_Edges_Append.Speed_d=G_Edges_Speed_d_Append;
        G_Edges_Append.Distance=G_Edges_Distance_Append;
        G_Edges_Append.PCI1=G_Edges_PCI1_Append;
        G_Edges_Append.PCI2=G_Edges_PCI2_Append;
        % G_Edges_Append.Num_lane=G_Edges_Num_lane_Append;
        G_Edges_Append.X1i=G_Edges_X1i_Append;
        G_Edges_Append.X2i=G_Edges_X2i_Append;
        G_Edges_Append.vi1=G_Edges_vi1_Append;
        G_Edges_Append.vi2=G_Edges_vi2_Append;
        % G_Edges_Append.Ti=G_Edges_Ti_Append;
        G_Edges_Append_temp{jj}=G_Edges_Append;
     end
 %% ��������ͼ��Nodes������������½ڵ�
    G_Nodes_New=table;
    G_Nodes_NODEID_New=[G_Nodes_NODEID;node5;node5+1;node5+2;node5+3];    
    G_Nodes_YoN_New=[G_Nodes_YoN;0;0;0;0];%������YoN�趨Ϊ0
    G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node3(1);XDate_node4(1);XDate_node3(2);XDate_node4(2)];
    G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node3(1);YDate_node4(1);YDate_node3(2);YDate_node4(2)];
    G_Nodes_New.NODEID=G_Nodes_NODEID_New;
    G_Nodes_New.YoN=G_Nodes_YoN_New;
    G_Nodes_New.XDate=G_Nodes_XDate_New;
    G_Nodes_New.YDate=G_Nodes_YDate_New;
     %% ��ԭ�е�G_Edgesɾ������·��λ�õ�ԭʼ·�Σ��õ�G_Edges_0��18��19,19��18��
     G_Edges_0=G_Edges;
    for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        for s=1:size(G_Edges_0.EndNodes,1) 
            if G_Edges_0.EndNodes(s,1)==node1&&G_Edges_0.EndNodes(s,2)==node2
                p=s;
            end
        end
        G_Edges_0(p,:)=[];
    end
    %% �������������·�κϲ���������6��·��
    G_Edges_Add=[G_Edges_Append_temp{1};G_Edges_Append_temp{2}];
    %% ���������·��table
    G_Edges_New_temp=cell(m,1);
    G_Edges_New_temp{ii}=[G_Edges_0;G_Edges_Add];
    %% �����µ�����ͼ
    G_new=cell(m,1);
    G_new{ii}=digraph(G_Edges_New_temp{ii},G_Nodes_New);
    % plot(G_new);
end
G_save=G_new{m};
end
%����������ĩ���G_LaneProperty_T_v
function [G_Flow,G_FlowD_Separator]=UpdateFlow_LinkSeparator3(G_LaneProperty_T_v,paramater_plus,time_event)
%��ͬʱ��·��·�λ���ͳһ��������ȡ�������������ں������(Xi��PCI)
%���2����Cap_d������Ǩ��
n=size(time_event,1);
paramater_r=cell(n-1,1);
G_Flow_Edges=cell(n,1);
G_FlowD_Separator=cell(n,1);
G_Flow_Nodes=cell(n,1);
G_Flow=cell(n,1);
for i=1:n
    if i==n
        G_FlowD_Separator{i}=G_LaneProperty_T_v{i};
    else
        paramater_temp=[];
        cn=0;
        for j=i+1:n
            if time_event(j,2)==1
                cn=cn+1;
                paramater_temp(cn,:)=paramater_plus{j};
            end
        end
        paramater_r{i}=paramater_temp;
        if isempty(paramater_r{i})==1
            G_FlowD_Separator{i}=G_LaneProperty_T_v{i};
        else
            G_FlowD_temp=G_LaneProperty_T_v{i};
            G_FlowD_Separator{i}=UpdateFlow_LinkSeparator_A4_MSA(G_FlowD_temp,paramater_r{i});
%             plot(G_FlowD_Separator{i})
        end
    end
    G_Flow_Edges{i}=table;
    G_Flow_Edges{i}=G_FlowD_Separator{i}.Edges(:,{'EndNodes','X1i','X2i','PCI1','PCI2'});%����ȡ�������������ں������
    G_Flow_Nodes{i}=G_FlowD_Separator{i}.Nodes;
    G_Flow{i}=digraph(G_Flow_Edges{i},G_Flow_Nodes{i});
end
end
%UpdateFlow_LinkSeparator��function���������ǰ·��·�λ���ͳһ
function G_save=UpdateFlow_LinkSeparator_A4_MSA(G_LaneProperty,paramater_r)%eg38����Ķ�dertaX
%�����������Ժ󣬶�ǰ���·�����зֶΣ�ʹ�ò�ͬʱ�ε�·������ֱ����ӵõ�·������

% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
m=size(paramater_r,1);
for ii=1:m
    %% 
    paramater_r_temp(1,:)=paramater_r(ii,:);
    paramater_r_temp(2,:)=[paramater_r(ii,2),paramater_r(ii,1),paramater_r(ii,3),paramater_r(ii,4),paramater_r(ii,5)+2];
    if ii==1
        G=G_LaneProperty;
    else
        G=G_new{ii-1};
    end
     %% 
    G_Edges=G.Edges;
    G_Nodes=G.Nodes;
    G_Edges_EndNodes=G_Edges.EndNodes;
    G_Edges_Cap_d1=G_Edges.Cap_d1;
    G_Edges_Cap_d2=G_Edges.Cap_d2;
    G_Edges_Speed_d=G_Edges.Speed_d;
    G_Edges_Distance=G_Edges.Distance;
%     G_Edges_par_pavement=G_Edges.par_pavement;
    G_Edges_PCI1=G_Edges.PCI1;
    G_Edges_PCI2=G_Edges.PCI2;
    G_Edges_X1i=G_Edges.X1i;
    G_Edges_X2i=G_Edges.X2i;
%     G_Edges_N1i=G_Edges.N1i;
%     G_Edges_N2i=G_Edges.N2i;
    G_Edges_v1i=G_Edges.vi1;
    G_Edges_v2i=G_Edges.vi2;
    G_Nodes_NODEID=G_Nodes.NODEID;
    G_Nodes_YoN=G_Nodes.YoN;
    G_Nodes_XData=G_Nodes.XDate;  
    G_Nodes_YData=G_Nodes.YDate;  
    %% ��������ͼ�ظ����Σ�1�α����Σ����������ڵ�
    G_Edges_Append_temp=cell(2,1);
    dertaX=[-40,80];%eg38���꣬��ͼ������ƫ��
     for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        % node1-node2����L
        for i=1:size(G_Nodes_NODEID)
            if G_Nodes_NODEID(i)==node1
                nodepossion1=i;
            end
        end
        for i=1:size(G_Nodes_NODEID)
            if G_Nodes_NODEID(i)==node2
                nodepossion2=i;
            end
        end
        XData_node1=G_Nodes_XData(nodepossion1);
        XData_node2=G_Nodes_XData(nodepossion2);
        YData_node1=G_Nodes_YData(nodepossion1);
        YData_node2=G_Nodes_YData(nodepossion2);
        a=[XData_node1,YData_node1];
        b=[XData_node2,YData_node2];
        L_node12=norm(b-a);
        %
        if jj==1
            fd=paramater_r_temp(jj,3);
        else
            fd=L_node12-paramater_r_temp(jj,3)-paramater_r_temp(jj,4);
        end
        length=paramater_r_temp(jj,4);
        XData_node3(jj)= XData_node1+ (XData_node2 -XData_node1)/L_node12*fd;
        YData_node3(jj)= YData_node1 + (YData_node2 -YData_node1)/L_node12*fd;
        XData_node4(jj)= XData_node1 + (XData_node2 -XData_node1)/L_node12*(fd+length);
        YData_node4(jj)= YData_node1 + (YData_node2 -YData_node1)/L_node12*(fd+length);
        dertaY(jj)=dertaX(jj)*abs(XData_node1-XData_node2)/abs(YData_node1-YData_node2);
        
        %% ��������ͼ��Edges���Բ����½ڵ㣬ʹ��һ�α������
        node5=paramater_r_temp(jj,5);
        for s=1:size(G_Edges_EndNodes,1) 
            if G_Edges_EndNodes(s,1)==node1&&G_Edges_EndNodes(s,2)==node2
                PosisionNumber=s;
            end
        end
        G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),node5;
            node5,node5+1;
            node5+1,G_Edges_EndNodes(PosisionNumber,2)];
        G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber,:);
            G_Edges_Cap_d1(PosisionNumber,:);
            G_Edges_Cap_d1(PosisionNumber,:)];
        G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber,:);
            G_Edges_Cap_d2(PosisionNumber,:);
            G_Edges_Cap_d2(PosisionNumber,:)];
        G_Edges_Distance_Append=[G_Edges_Distance(PosisionNumber,:);
            G_Edges_Distance(PosisionNumber,:);
            G_Edges_Distance(PosisionNumber,:)];
        G_Edges_Speed_d_Append=[G_Edges_Speed_d(PosisionNumber,:);
            G_Edges_Speed_d(PosisionNumber,:);
            G_Edges_Speed_d(PosisionNumber,:)];
%         G_Edges_par_pavement_Append=[G_Edges_par_pavement(PosisionNumber,:);
%             G_Edges_par_pavement(PosisionNumber,:);
%             G_Edges_par_pavement(PosisionNumber,:)];
        G_Edges_PCI1_Append=[G_Edges_PCI1(PosisionNumber,:);
            G_Edges_PCI1(PosisionNumber,:);
            G_Edges_PCI1(PosisionNumber,:)];
        G_Edges_PCI2_Append=[G_Edges_PCI2(PosisionNumber,:);
            G_Edges_PCI2(PosisionNumber,:);
            G_Edges_PCI2(PosisionNumber,:)];
        G_Edges_X1i_Append=[G_Edges_X1i(PosisionNumber,:);
            G_Edges_X1i(PosisionNumber,:);
            G_Edges_X1i(PosisionNumber,:)];
        G_Edges_X2i_Append=[G_Edges_X2i(PosisionNumber,:);
            G_Edges_X2i(PosisionNumber,:);
            G_Edges_X2i(PosisionNumber,:)];
%         G_Edges_N1i_Append=[G_Edges_N1i(PosisionNumber,:);
%             G_Edges_N1i(PosisionNumber,:);
%             G_Edges_N1i(PosisionNumber,:)];
%         G_Edges_N2i_Append=[G_Edges_N2i(PosisionNumber,:);
%             G_Edges_N2i(PosisionNumber,:);
%             G_Edges_N2i(PosisionNumber,:)];
        G_Edges_v1i_Append=[G_Edges_v1i(PosisionNumber,:);
            G_Edges_v1i(PosisionNumber,:);
            G_Edges_v1i(PosisionNumber,:)];
        G_Edges_v2i_Append=[G_Edges_v2i(PosisionNumber,:);
            G_Edges_v2i(PosisionNumber,:);
            G_Edges_v2i(PosisionNumber,:)];
        %% �������ӵ�·���������һ���µ�table_append
        G_Edges_Append=table;
        G_Edges_Append.EndNodes=G_Edges_EndNodes_Append;
        G_Edges_Append.Cap_d1=G_Edges_Cap_d1_Append;
        G_Edges_Append.Cap_d2=G_Edges_Cap_d2_Append;
        G_Edges_Append.Speed_d=G_Edges_Speed_d_Append;
        G_Edges_Append.Distance=G_Edges_Distance_Append;
%         G_Edges_Append.par_pavement=G_Edges_par_pavement_Append;
        G_Edges_Append.PCI1=G_Edges_PCI1_Append;
        G_Edges_Append.PCI2=G_Edges_PCI2_Append;
        G_Edges_Append.X1i=G_Edges_X1i_Append;
        G_Edges_Append.X2i=G_Edges_X2i_Append;
%         G_Edges_Append.N1i=G_Edges_N1i_Append;
%         G_Edges_Append.N2i=G_Edges_N2i_Append;
        G_Edges_Append.vi1=G_Edges_v1i_Append;
        G_Edges_Append.vi2=G_Edges_v2i_Append;        
        G_Edges_Append_temp{jj}=G_Edges_Append;
     end
 %% ��������ͼ��Nodes�������2+2���½ڵ�
         %������ͼ
        G_Nodes_New=table;
        G_Nodes_NODEID_New=[G_Nodes_NODEID;node5;node5+1;node5+2;node5+3];    
        G_Nodes_YoN_New=[G_Nodes_YoN;0;0;0;0];%������YoN�趨Ϊ0
        G_Nodes_XData_New=[G_Nodes_XData;XData_node3(1)-dertaX(1);XData_node4(1)-dertaX(1);XData_node3(2)-dertaX(2);XData_node4(2)-dertaX(2)];%eg38����
        G_Nodes_YData_New=[G_Nodes_YData;YData_node3(1)-dertaY(1);YData_node4(1)-dertaY(1);YData_node3(2)-dertaY(2);YData_node4(2)-dertaY(2)];%eg38����
        G_Nodes_New.NODEID=G_Nodes_NODEID_New;
        G_Nodes_New.YoN=G_Nodes_YoN_New;
        G_Nodes_New.XDate=G_Nodes_XData_New;
        G_Nodes_New.YDate=G_Nodes_YData_New;
     %% ��ԭ�е�G_Edgesɾ������·��λ�õ�ԭʼ·�Σ��õ�G_Edges_0��18��19,19��18��
     G_Edges_0=G_Edges;
    for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        for s=1:size(G_Edges_0.EndNodes,1) 
            if G_Edges_0.EndNodes(s,1)==node1&&G_Edges_0.EndNodes(s,2)==node2
                p=s;
            end
        end
        G_Edges_0(p,:)=[];
    end
    %% �������������·�κϲ���������6��·��
    G_Edges_Add=[G_Edges_Append_temp{1};G_Edges_Append_temp{2}];
    %% ���������·��table
    G_Edges_New_temp=cell(m,1);
    G_Edges_New_temp{ii}=[G_Edges_0;G_Edges_Add];
    %% �����µ�����ͼ
    G_new=cell(m,1);
    G_new{ii}=digraph(G_Edges_New_temp{ii},G_Nodes_New);
    % plot(G_new);
end
G_save=G_new{m};
end
%�������
function G_FlowSum_derta=G_FlowSum_MSA(T,G_Flow,time_event)
%�������󣬷ֶμ���·��������
% ����PCI��������T��������Ӧ��·�������ܺͣ�eg.PCI��������Ϊ3m����ôFlowSum��������ҲΪ3m
%G_Property_Separator,time_eventΪ��Ӧ�����Ժ�ʱ���б�time_event[i-1,i)=G_Property_Separator{i}
time=[0;time_event(:,1)];
n=size(time,1);
m=fix(time(n,1)/T);
u=30*24;
%���������������������ͣ�
time_d=time;
for i=2:n
time_d(i)=time(i)-time(i-1);
end
for i=1:n-1
    X1i_area{i}=G_Flow{i}.Edges.X1i.*time_d(i+1)*u;
    X2i_area{i}=G_Flow{i}.Edges.X1i.*time_d(i+1)*u;
    if i==1
        X1i_areaSum{i}=X1i_area{i};
        X2i_areaSum{i}=X2i_area{i};
    else
        X1i_areaSum{i}=X1i_areaSum{i-1}+X1i_area{i};
        X2i_areaSum{i}=X2i_areaSum{i-1}+X2i_area{i};
    end
end
%�����ڼ���������
for j=1:m
    G_FlowSum{j}=G_Flow{1};
    t=T*j;
    %����ʱ��t��������time��λ��
    for ii=1:n
        if t>time(ii)&&t<=time(ii+1)
            possionnum_t=ii;
        end
    end
    %�������������Ĳ���
    X1i_fix(:,j)=G_Flow{possionnum_t}.Edges.X1i.*(t-time(possionnum_t))*u;
    X2i_fix(:,j)=G_Flow{possionnum_t}.Edges.X2i.*(t-time(possionnum_t))*u;
    if possionnum_t==1
        X1i_Sum{j}=X1i_fix(:,j);
        X2i_Sum{j}=X2i_fix(:,j);
    else
        X1i_Sum{j}=X1i_areaSum{possionnum_t-1}+X1i_fix(:,j);
        X2i_Sum{j}=X2i_areaSum{possionnum_t-1}+X2i_fix(:,j);
    end
    G_FlowSum{j}.Edges.X1i=X1i_Sum{j};
    G_FlowSum{j}.Edges.X2i=X2i_Sum{j};
end
%������������
G_FlowSum_derta{1}=G_FlowSum{1};
for j=2:m
    G_FlowSum_derta{j}=G_FlowSum{j};
    G_FlowSum_derta{j}.Edges.X1i=G_FlowSum{j}.Edges.X1i-G_FlowSum{j-1}.Edges.X1i;
    G_FlowSum_derta{j}.Edges.X2i=G_FlowSum{j}.Edges.X2i-G_FlowSum{j-1}.Edges.X2i;      
end
end
%% G_Flow function2
function G_FlowD_Separator=UpdateFlow_LinkSeparator2(G_FlowD,paramater_plus,time_event)
%��ͬʱ��·��·�λ���ͳһ��������ȡ�������������ں������(Xi��PCI)
n=size(time_event,1);
paramater_r=cell(n-1,1);
G_Flow_Edges=cell(n,1);
G_FlowD_Separator=cell(n,1);
G_Flow_Nodes=cell(n,1);
G_Flow=cell(n,1);
for i=1:n
    i
    if i==n
        G_FlowD_Separator{i}=G_FlowD{i};
    else
        paramater_temp=[];
        cn=0;
        for j=i+1:n
            if time_event(j,2)==1
                cn=cn+1;
                paramater_temp(cn,:)=paramater_plus{j};
            end
        end
        paramater_r{i}=paramater_temp;
        if isempty(paramater_r{i})==1
            G_FlowD_Separator{i}=G_FlowD{i};
        else
            G_FlowD_temp=G_FlowD{i};
            G_FlowD_Separator{i}=UpdateFlow_LinkSeparator_A2(G_FlowD_temp,paramater_r{i});
%             plot(G_FlowD_Separator{i})
        end
    end
%     G_Flow_Edges{i}=table;
%     G_Flow_Edges{i}=G_FlowD_Separator{i}.Edges(:,{'EndNodes','X1i','X2i','PCI1','PCI2'});%����ȡ�������������ں������
%     G_Flow_Nodes{i}=G_FlowD_Separator{i}.Nodes;
%     G_Flow{i}=digraph(G_Flow_Edges{i},G_Flow_Nodes{i});
end
end
function G_save=UpdateFlow_LinkSeparator_A2(G_FlowDi,paramater_r)
%�����������Ժ󣬶�ǰ���·�����зֶΣ�ʹ�ò�ͬʱ�ε�·������ֱ����ӵõ�·������

% �޸Ľڵ�ĺ��������G_save��һ��Ԫ�����飬G_save{1}��ԭʼ�ṹ��G_save{2}�Ǹ��º�Ľṹ��ʹ��ʱ������Ҫ����
%node1��node2�Ƕ�Ӧ·�����з����Եģ�node1��node2����Ϊǰ������
%fd,length�ֱ�Ϊ��������node1�ľ���͹���������
% node1=1;node2=2;r
n1=size(paramater_r,1);
for ii=1:n1
    paramater_r_temp(1,:)=paramater_r(ii,:);
    paramater_r_temp(2,:)=[paramater_r(ii,2),paramater_r(ii,1),paramater_r(ii,4),paramater_r(ii,3),paramater_r(ii,5)+2];
    if ii==1
        G=G_FlowDi;
    else
        G=G_new{ii-1};
    end
    %% 
    G_Edges=G.Edges;
    G_Nodes=G.Nodes;
    G_Edges_EndNodes=G_Edges.EndNodes;

%     G_Edges_Cap_d1=G_Edges.Cap_d1;
%     G_Edges_Cap_d2=G_Edges.Cap_d2;
%     G_Edges_Speed_d=G_Edges.Speed_d;
%     G_Edges_Distance=G_Edges.Distance;
%     G_Edges_PCI1=G_Edges.PCI1;
%     G_Edges_PCI2=G_Edges.PCI2;
%     G_Edges_Ti=G_Edges.Ti;

    G_Edges_X1i=G_Edges.X1i;
    G_Edges_X2i=G_Edges.X2i;
    G_Edges_vi1=G_Edges.vi1;
    G_Edges_vi2=G_Edges.vi2;
    G_Nodes_NODEID=G_Nodes.NODEID;
    G_Nodes_YoN=G_Nodes.YoN;
    
    %% ��������ͼ�ظ����Σ�1�α����Σ����������ڵ�
    G_Edges_Append_temp=cell(2,1);
     for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        fd=paramater_r_temp(jj,3);
        length=paramater_r_temp(jj,4);
        node5=paramater_r_temp(jj,5);
        for s=1:size(G_Edges_EndNodes,1) 
            if G_Edges_EndNodes(s,1)==node1&&G_Edges_EndNodes(s,2)==node2
                PosisionNumber=s;
            end
        end
        %% ��������ͼ��Nodes������������½ڵ�
        G_Nodes_New=table;
        G_Nodes_NODEID_New=[G_Nodes_NODEID;node5;node5+1;node5+2;node5+3];
        G_Nodes_YoN_New=[G_Nodes_YoN;0;0;0;0];%������YoN�趨Ϊ0
        G_Nodes_New.NODEID=G_Nodes_NODEID_New;
        G_Nodes_New.YoN=G_Nodes_YoN_New;
        %% ��������ͼ��Edges���Բ����½ڵ㣬ʹ��һ�α������
        G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),node5;
            node5,node5+1;
            node5+1,G_Edges_EndNodes(PosisionNumber,2)];
%         G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
%             G_Edges_Weight(PosisionNumber);
%             G_Edges_Weight(PosisionNumber)];
%         G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
%             G_Edges_Cap_d1(PosisionNumber);
%             G_Edges_Cap_d1(PosisionNumber)];
%         G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
%             G_Edges_Cap_d2(PosisionNumber);
%             G_Edges_Cap_d2(PosisionNumber)];
%         G_Edges_Speed_d_Append=[G_Edges_Speed_d(PosisionNumber);
%             G_Edges_Speed_d(PosisionNumber);
%             G_Edges_Speed_d(PosisionNumber)];
%         G_Edges_Distance_Append=[fd;
%             G_Edges_Distance(PosisionNumber)-fd-length;
%             length];
%          G_Edges_Ti_Append=[G_Edges_Ti(PosisionNumber);
%             G_Edges_Ti(PosisionNumber);
%             G_Edges_Ti(PosisionNumber)];
%         G_Edges_PCI1_Append=[G_Edges_PCI1(PosisionNumber);
%             G_Edges_PCI1(PosisionNumber);
%             G_Edges_PCI1(PosisionNumber)];
%         G_Edges_PCI2_Append=[G_Edges_PCI2(PosisionNumber);
%              G_Edges_PCI2(PosisionNumber);
%             G_Edges_PCI2(PosisionNumber)];

        G_Edges_X1i_Append=[G_Edges_X1i(PosisionNumber,:);
            G_Edges_X1i(PosisionNumber,:);
            G_Edges_X1i(PosisionNumber,:)];
        G_Edges_X2i_Append=[G_Edges_X2i(PosisionNumber,:);
            G_Edges_X2i(PosisionNumber,:);
            G_Edges_X2i(PosisionNumber,:)];
        G_Edges_vi1_Append=[G_Edges_vi1(PosisionNumber,:);
            G_Edges_vi1(PosisionNumber,:);
            G_Edges_vi1(PosisionNumber,:)];
        G_Edges_vi2_Append=[G_Edges_vi2(PosisionNumber,:);
            G_Edges_vi2(PosisionNumber,:);
            G_Edges_vi2(PosisionNumber,:)];

        
        %% �������ӵ�·���������һ���µ�table_append
        G_Edges_Append=table;
        G_Edges_Append.EndNodes=G_Edges_EndNodes_Append;
%         G_Edges_Append.Cap_d1=G_Edges_Cap_d1_Append;
%         G_Edges_Append.Cap_d2=G_Edges_Cap_d2_Append;
%         G_Edges_Append.Speed_d=G_Edges_Speed_d_Append;
%         G_Edges_Append.Distance=G_Edges_Distance_Append;
%         G_Edges_Append.PCI1=G_Edges_PCI1_Append;
%         G_Edges_Append.PCI2=G_Edges_PCI2_Append;
%         G_Edges_Append.Ti=G_Edges_Ti_Append;
        
        G_Edges_Append.X1i=G_Edges_X1i_Append;
        G_Edges_Append.X2i=G_Edges_X2i_Append;
        G_Edges_Append.vi1=G_Edges_vi1_Append;
        G_Edges_Append.vi2=G_Edges_vi2_Append;
        
        G_Edges_Append_temp{jj}=G_Edges_Append;
     end
     %% ��ԭ�е�G_Edgesɾ������·��λ�õ�ԭʼ·�Σ��õ�G_Edges_0��18��19,19��18��
     G_Edges_0=G_Edges;
    for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        for s=1:size(G_Edges_0.EndNodes,1) 
            if G_Edges_0.EndNodes(s,1)==node1&&G_Edges_0.EndNodes(s,2)==node2
                p=s;
            end
        end
        G_Edges_0(p,:)=[];
    end
    %% �������������·�κϲ���������6��·��
    G_Edges_Add=[G_Edges_Append_temp{1};G_Edges_Append_temp{2}];
    %% ���������·��table
    G_Edges_New_temp=cell(n1,1);
    G_Edges_New_temp{ii}=[G_Edges_0;G_Edges_Add];
    %% �����µ�����ͼ
    G_new=cell(n1,1);
    G_new{ii}=digraph(G_Edges_New_temp{ii},G_Nodes_New);
    % plot(G_new);
end
G_save=G_new{n1};
end

%% ��η���function
%�ؼ�1���г�ʱ����Ҫ��θ���
%�ؼ�2�����ζ�·���������dail����õ�Xi
%% ���շ����е��μ��㺯����·����,MSA���壬��Ҫ���ظ�����dail����
function G_save=day2day_mG(m,theta_input,G,OD_demand,erfa,beita,epsilon)
% G=H0;
% OD_demand=OD_input_p;
%% �������
EdNodes=G.Edges.EndNodes;
G_Nodes=G.Nodes;
TNodes=EdNodes(:,1);
ENodes=EdNodes(:,2);
G_vexnum=length(EdNodes);%����
Ci=G.Edges.Cap_d1+G.Edges.Cap_d2;
Di=G.Edges.Distance;
Xi_0=G.Edges.X1i+G.Edges.X2i;
Spi=G.Edges.Speed_d;
Wi0=Di./Spi.*3.6 ;
preca=zeros(G_vexnum,m);
ca_temp=zeros(G_vexnum,m);
Xi_temp=zeros(G_vexnum,m);
vi_temp=zeros(G_vexnum,m);
Xi_sum=zeros(G_vexnum,1);
%% ��μ��㣬���ж��Ƿ�ﵽ�ȶ���ʲôʱ��ﵽ�ȶ�
for ii=1:m
    if ii==1
        erfa_temp=0;
        G_temp=G;%����Xi����ӦWi
        preca_last=ones(G_vexnum,1);%��Ϊerfaȡ��0��Wi_lastȡ�κ������һ����������������ó���1����
        Xi_last=Xi_0;
        ca_last=cost_set(Wi0,Xi_last,Ci);
    else
        erfa_temp=erfa;
        G_temp=H{ii-1};%Ԥ�������Ҫ���ÿ��G����������޸�ΪH{ii}
        preca_last=preca(:,ii-1);
        ca_last=ca_temp(:,ii-1);
    end
    preca(:,ii)=ca_last.*(1-erfa_temp)+preca_last.*erfa_temp;%Ԥ���г�ʱ��
    [count1,count2]=dial_G(theta_input,G_temp,OD_demand,preca(:,ii));
    Y1=zeros(G_vexnum,1);%��ԭ������ͼ��ʽ
    for i=1:G_vexnum
        a=TNodes(i);
        b=ENodes(i);
        Y1(i,1)=round(count1(a,b));
    end
    %���Ǧ¼��㵥�η�������
    Xi_temp(:,ii)=Y1.*beita+Xi_last.*(1-beita);
    G_temp.Edges.X1i=round(Xi_temp(:,ii)./2);
    G_temp.Edges.X2i=round(Xi_temp(:,ii)./2);
    H{ii}=G_temp;
    %ʵ���г�ʱ��
    ca_temp(:,ii)=cost_set(Wi0,Xi_temp(:,ii),Ci);
    % ·��ƽ������
    vi_temp(:,ii)=Di./ca_temp(:,ii)*3.6;
end

    %% ����������״̬������X��Speed��������һ���ڳ�ʼ��
    G_temp2=G_temp;
    G_temp2.Edges.X1i=round(Xi_temp(:,m)./2);
    G_temp2.Edges.X2i=round(Xi_temp(:,m)./2);
    G_temp2.Edges.vi1=round(vi_temp(:,m));
    G_temp2.Edges.vi2=round(vi_temp(:,m)); 
    G1{1}=G_temp2;
    %% ����X��Speed�洢
    X1i_temp=round(Xi_temp./2);
    X2i_temp=round(Xi_temp./2);
    vi1_temp=round(vi_temp);
    vi2_temp=round(vi_temp);   
    EdgeTable_2 =  table(EdNodes,X1i_temp,X2i_temp,vi1_temp,vi2_temp,'VariableNames',{'EndNodes','X1i','X2i','vi1','vi2'});
    G1{2}=digraph(EdgeTable_2,G_Nodes);
     %% �����ϼƼ�ƽ�����ٴ����µ�����ͼ
    Xi_sum=sum(Xi_temp,2);%������
    vi_mu=sum(vi_temp,2)./m;%����ʱ��ƽ������
    G_temp.Edges.X1i=round(Xi_sum./2);
    G_temp.Edges.X2i=round(Xi_sum./2);
    G_temp.Edges.vi1=vi_mu;
    G_temp.Edges.vi2=vi_mu;
    % G.Edges.Ti=Ti;
    G1{3}=G_temp;
    G_save=G1;
     %% step4 �жϵ�ǰ·�������Ƿ�������������
     error= sum(abs(Xi_temp(:,m)-Xi_temp(:,m-1)))/sum(Xi_temp(:,m));
     if error<epsilon
         fprintf('�ȶ�,����%d,���error=%8.5f\n',m,error)
     else
         fprintf('���ȶ�,����%d,���error=%8.5f\n',m,error)
     end
 
end

%% MSA function
function [H] =MSA_G(theta_input,G,OD_demand,erfa,beita,epsilon)
%����ԭ����MSA_G�����ڳ�ʼ״̬����
%step1������dail���䣬dail�г��гɱ���t0���棬�����X1��
%step2����X���³��гɱ�����ΪȨ�أ�����0-1��������Y2��
%step3:�ü�Ȩ�����ķ������õ�X2,��X2ΪX1��Y2-X1�ļ�Ȩ�ͣ�
%step4���Ƚ�X2��X1���ж��Ƿ������������Ҫ�󣬲��������ظ�2.3.4step��
%�÷�����ֻ������һ��dail���㣬����ԭ��Ϊ��μ�Ȩ����
%erfa������ʱ����µ�ϰ����Ȩ��
% G=G0;
% OD_demand=OD_input_p;
error=0.3;%���ó�ʼ�������
k=1;%����������������Ϊ����Ȩ��
x_p=0.5;%����˫�����ڳ���������������
OD_Nodes=G.Nodes.NODEID(logical(G.Nodes.YoN));%ODʱ���������������������˽ṹ�ı�Ŷ�Ӧ��
G_ODnum=length(OD_Nodes);
topo_Nodes=G.Nodes.NODEID;%���˽ṹ���нڵ�ı��
EdNodes=G.Edges.EndNodes;
TNodes=EdNodes(:,1);
ENodes=EdNodes(:,2);
G_vexnum=length(EdNodes);%����
G_nodenum=length(topo_Nodes);%�ڵ���
% Num_lane= G.Edges.Num_lane;
Ci=G.Edges.Cap_d1+G.Edges.Cap_d2;
Di=G.Edges.Distance;
Spi=G.Edges.Speed_d;
Wi0 =Di./Spi.*3.6 ;%����speed_d�Լ�distance�й�,��Ҫ���£�����ԭ����Wi = G.Edges.Weight
% Xi=G.Edges.X1i+G.Edges.X2i;
%% ·��״̬��ʼ��
%dial����
%��ʼ��·������X1
[count1,count2]=dial_G(theta_input,G,OD_demand,Wi0);
X1=zeros(G_vexnum,1);%��ԭ������ͼ��ʽ
for i=1:G_vexnum
    a=TNodes(i);
    b=ENodes(i);
    X1(i,1)=round(count1(a,b));
end
%��ʼ���г�ʱ��

time_temp= cost_set(Wi0,X1,Ci);
while error>epsilon
    %% step1������·����ʻʱ��time
    time_temp=precost_update(Wi0,X1,Ci,time_temp,erfa);
    %% step2: ��step1�м������·��ʱ���OD��ͨ����ִ��һ��dail���䣬�õ�һ�鸽�ӽ�ͨ��Y1
    [count1,count2]=dial_G(theta_input,G,OD_demand,time_temp);
    Y1=zeros(G_vexnum,1);
    for i=1:G_vexnum
        a=TNodes(i);
        b=ENodes(i);
        Y1(i,1)=round(count1(a,b));
    end
    %% step3 �����·�ε�ǰ��ͨ��
    X2=X1+1/k.*(Y1-X1); %   0< phi<1 phi=0.%���ٵ�������������û�дﵽ������SUE״̬
%     X2=X1+beita.*(Y1-X1); %   ���Դﵽ���⣬���Ǽ������̫��
    %% step4 �жϵ�·�������Ƿ�������������
    error= sum(abs(X1-X2))/sum(X1);
    %��������
    X1=X2;
    k=k+1
end
 fprintf('��������%d,���error=%8.5f\n',k,error)

%% ����·��ƽ���ٶ�
T1= cost_set(Wi0,X1,Ci);
V1=Di./T1*3.6;
%% ��������ͼ���
% G.Edges.Weight=Wi;
G.Edges.X1i=round(X1./2);
G.Edges.X2i=round(X1./2);
G.Edges.vi1=V1;
G.Edges.vi2=V1;
% G.Edges.Ti=Ti;
H=G;
% count=count_positive+count_negative';
% plotHH(count)
end

%MSA��function��dail/time_update/construct2
function [count1 count2]=dial_G(theta_input,G,OD_demand,W)
%����㵽�յ����·
% OD_demand,��ͨ����ODʱ�����
OD_Nodes=G.Nodes.NODEID(logical(G.Nodes.YoN));%ODʱ���������������������˽ṹ�ı�Ŷ�Ӧ��
topo_Nodes=G.Nodes.NODEID;%���˽ṹ���нڵ�ı��
EdNodes=G.Edges.EndNodes;
TNodes=EdNodes(:,1);
ENodes=EdNodes(:,2);
G_nodenum=length(topo_Nodes);%�ڵ���
G_ODnum=length(OD_Nodes);%OD�ڵ���
% W = G.Edges.Distance./G.Edges.Speed_d.*3.6;
L=zeros(1,length(EdNodes));
% g=biograph(sparse(TNodes,ENodes,W));
% view(g);
%% ��������нڵ�����·����r�����нڵ����յ�����·����s
G_graph=G;
G_graph.Edges.Weight=W;
for i=1:G_nodenum
    for j=1:G_nodenum
        [path,dist] = shortestpath(G_graph,topo_Nodes(i),topo_Nodes(j));%bigraph�γɵ����·������
%         [dist,path]=graphshortestpath(linkweight,i,j);%biograph�γɵ����·������
        dis(topo_Nodes(i),topo_Nodes(j))=dist;
    end
end

%���巢��
%for i=1:9
Oi=cell(1,G_nodenum);
Di=cell(1,G_nodenum);
for i=1:G_nodenum
    k=topo_Nodes(i);
    m=find(TNodes==k);
    Oi{k}=[ENodes(m)];%�ڵ���
end
%�����յ�
for i=1:G_nodenum
    k=topo_Nodes(i);
    m=find(ENodes==k);
    Di{k}=[TNodes(m)];%�ڵ���
end

%%������
count_positive0=sparse(TNodes,ENodes,L);
count_positive=spear2full_sq(count_positive0);
for i=1:G_ODnum
    origin=OD_Nodes(i);
    for j=1:G_ODnum
        if i~=j
            destination=OD_Nodes(j);
            r=dis(origin,:);%�ڵ���
            s=dis(:,destination)';
            Q=OD_demand(i,j);%λ�ñ��
            Lsparse=likelihood_G(r,s,theta_input,W,EdNodes);
            Wi=weight2_G(r,Lsparse,Oi,Di,origin,EdNodes,destination);
            Xi=flow2_G(s,EdNodes,Oi,Di,Wi,destination,Q,origin);
            count_positive=count_positive+Xi;
        end
    end
end

%%������
count_negative0=sparse(TNodes,ENodes,L);
count_negative=spear2full_sq(count_negative0);
for i=1:G_ODnum
    origin=OD_Nodes(i);
    for j=1:G_ODnum
        if i~=j
        destination=OD_Nodes(j);
        r=dis(origin,:);
        s=dis(:,destination)';
        Q=OD_demand(i,j);%jΪOD�������ţ��ǽڵ��
        Lsparse=likelihood_G(r,s,theta_input,W,EdNodes);
        Wi=weight2_G(r,Lsparse,Oi,Di,origin,EdNodes,destination);
        Xi=flow2_G(s,EdNodes,Oi,Di,Wi,destination,Q,origin);
        count_negative=count_negative+Xi;
        end
    end
count1=full(count_positive);
count2=full(count_negative)';
end
end
%BPR����
function [time] = cost_set(W0,X,Cap)%���������ٶȵ�ʱ�����ڼ���ʵ����ʻʱ��
erfa=0.15;
berta=4;
time=W0.*(1+erfa.*(X./Cap).^berta);
time(isnan(time))=0;
end
function [time] = cost_set1(W0,X,Cap)%���������ٶȵ�ʱ�����ڼ���ʵ����ʻʱ��
erfa=0.15;
berta=4;
Cap(Cap-X<=0)=0.001;
time=W0.*(1+erfa.*(X./Cap).^berta);
time(isnan(time))=0;
end
 %����ϰ����ѡ��ı����г̸�֪ʱ��
 function  preca=precost_update(Wi0,X_last,Cap,preca_last,erfa)
 %������ϵ����ϰ����ϵ����erfa=0��˵����������ǰһ��ʵ��x��������ǰһ��Ԥ��Ľ�������֮ǰ����ƽ��״̬������Ҳ����ǰһ��ʵ�ʵ��г�ʱ�䡣
  ca_last = cost_set(Wi0,X_last,Cap);
  preca=ca_last.*(1-erfa)+preca_last.*erfa;
 end
%dail ��function:likehood/weight/flow
function Xi=flow2_G(s,EdNodes,Oi,Di,Wi,destination,Q,origin)
%Nodes��·��G·���ڽӹ�ϵ=G.Edges.EndNodes
[a,b]=sort(s);
Xi0=sparse(EdNodes(:,1),EdNodes(:,2),zeros(1,length(EdNodes)));
Xi=spear2full_sq(Xi0);
b=b(1:find(b==origin));
 for i=b
        %if Oi{i}==9 
        if i==destination
        %if Oi{i}==destination 
            Xi(Di{i},i)=Q*Wi(Di{i},i)/sum(Wi(Di{i},i));
            %Xi(Di{i},i)=Wi(Di{i},i)*sum(Xi(i,Oi{i}))/sum(Wi(Di{i},i));
        %end
        else 
        %if  Wi(Di{i},i)~=0
        %else
            Xi(Di{i},i)=Wi(Di{i},i)*sum(Xi(i,Oi{i}))/sum(Wi(Di{i},i));
        %end
        end
    %end
 end
Xi(isnan(Xi))=0;     
end
function Lsparse=likelihood_G(r,s,theta,W,EdNodes)
%WΪG�е�·���迹Ȩ��weight=G.Edges.Weight
%Nodes=G.Edges.EndNodes��·�����յ���
n=length(EdNodes);
L=zeros(1,n);

for k=1:n
    if r(EdNodes(k,1))<r(EdNodes(k,2))&&s(EdNodes(k,1))>s(EdNodes(k,2))
        L(k)=exp(theta*(r(EdNodes(k,2))-r(EdNodes(k,1))-W(k)));
    end
end
Lsparse=sparse(EdNodes(:,1),EdNodes(:,2),L);
end
function Wi=weight2_G(r,Lsparse,Oi,Di,origin,EdNodes,destination)
%r=r(origin,:);
[a,b]=sort(r);%��r����
Wi0=sparse(EdNodes(:,1),EdNodes(:,2),zeros(1,length(EdNodes)));
Wi=spear2full_sq(Wi0);
b=b(1:find(b==destination));
for i=b
    if i==origin
          Wi(i,Oi{i})=Lsparse(i,Oi{i});
    else 
          Wi(i,Oi{i})=Lsparse(i,Oi{i})*sum(Wi(Di{i},i));
    end
    
end

end 
%MSA����function
function possiablePaths = findPath(Graph, partialPath, destination, partialWeight)
% findPath����������������п��ܵĴ�partialPath������destination��·������Щ·���в�������·
% Graph: ·��ͼ���������0��ʾ���ڵ�֮��ֱ����ͨ������ֵ��Ϊ·��Ȩֵ
% partialPath: ������·�������partialPath��һ��������ʾ���������ʼ��
% destination: Ŀ��ڵ�
% partialWeight: partialPath��Ȩֵ����partialPathΪһ����ʱ��partialWeightΪ0
pathLength = length(partialPath);
lastNode = partialPath(pathLength); %�õ����һ���ڵ�
nextNodes = find(0<Graph(lastNode,:) & Graph(lastNode,:)<inf); %����Graphͼ�õ����һ���ڵ����һ���ڵ�
GLength = length(Graph);
possiablePaths = [];
if lastNode == destination
    % ���lastNode��Ŀ��ڵ���ȣ���˵��partialPath���Ǵ��������Ŀ��ڵ��·�������ֻ����һ����ֱ�ӷ���
    possiablePaths = partialPath;
    possiablePaths(GLength + 1) = partialWeight;
    return;
elseif length( find( partialPath == destination ) ) ~= 0
    return;
end
%nextNodes�е���һ������0,����Ϊ����nextNodes(i)ȥ�����Ƚ��丳ֵΪ0
for i=1:length(nextNodes)
    if destination == nextNodes(i)
        %���·��
        tmpPath = cat(2, partialPath, destination);      %���ӳ�һ��������·��
        tmpPath(GLength + 1) = partialWeight + Graph(lastNode, destination); %�ӳ����鳤����GLength+1, ���һ��Ԫ�����ڴ�Ÿ�·������·��
        possiablePaths( length(possiablePaths) + 1 , : ) = tmpPath;
        nextNodes(i) = 0;
    elseif length( find( partialPath == nextNodes(i) ) ) ~= 0
        nextNodes(i) = 0;
    end
end
nextNodes = nextNodes(nextNodes ~= 0); %��nextNodes��Ϊ0��ֵȥ������Ϊ��һ���ڵ�����Ѿ�����������������Ŀ��ڵ�
for i=1:length(nextNodes)
    tmpPath = cat(2, partialPath, nextNodes(i));
    tmpPsbPaths = findPath(Graph, tmpPath, destination, partialWeight + Graph(lastNode, nextNodes(i)));
    possiablePaths = cat(1, possiablePaths, tmpPsbPaths);
end
end
function [matrix_full] = spear2full_sq(matrix_sparse)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
a=matrix_sparse;
num=max(size(a,1),size(a,2));    
[m,n]=size(a);
if m>n
    s=m-n;
    for p=1:s
    for i=m:num
        for j=n+p:num
            a(i,j)=0;
        end
    end
    end
else
    s=n-m;
    for p=1:s
        for i=n:num 
            for j=m+p:num
                a(i,j)=0;
            end
        end
    end
end
matrix_full=a;
end

%% ���շ����е��μ��㺯��(˫·��)
function [ca,X,erfa,beita]=daytoday(erfa,m,du_s1,c2,t0,t1,thita,q,beita)
%������ϵ��erfa=1��˵����������ǰһ��ʵ��x
ca = cell(1,m,length(t0));%·��c1��Ϊnodenum
X = cell(1,m,length(t0));%·��c1��Ϊnodenum
for i=1:m
    if i==1
        erfa_temp=0;
        x_last=du_s1;
        ca_prlast=t1;
    else
        x_last=X{i-1};
        ca_prlast=ca{i-1};
        erfa_temp=erfa;
    end
        [ca_last]=time_1(x_last,c2,t0);
        [ca{i}]=preca(ca_last,ca_prlast,erfa_temp);
        [X{i}]=logit_x(q,ca{i},thita,x_last,beita);% ���η��� 
end
end
%daytoday��function��time_1��preca��logit_x
%��һ�ε�ʵ������Ԥ�Ȿ���г�ʱ��
 function [ca]=time_1(x,cap,t0)
 n=length(x);
 ca=zeros(1,n);
 for i=1:n
ca(1,i)=t0(1,i)*(1+0.15*(x(1,i)/cap(1,i))^4);
 end
 end
 %����ϰ����ѡ��ı����г̸�֪ʱ��
 function  [ca]=preca(ca_last,ca_prlast,erfa)
 %������ϵ����ϰ����ϵ����erfa=0��˵����������ǰһ��ʵ��x��������ǰһ��Ԥ��Ľ�������֮ǰ����ƽ��״̬������Ҳ����ǰһ��ʵ�ʵ��г�ʱ�䡣
  n=length(ca_last);
 ca=zeros(n,1);
 for i=1:n
ca(1,i)=ca_last(1,i)*(1-erfa)+ca_prlast(1,i)*erfa;
 end
 end
 %˫·��logit���η���
 function [X]=logit_x(Q,ca,thita,y,beita)
ca_mean=(ca(1,1)+ca(1,2))/2;
ca_sum=(exp(-thita*ca(1,1)/ ca_mean))+exp(-thita*ca(1,2)/ ca_mean);
ptrans1=exp(-thita*ca(1,1)/ ca_mean)/ca_sum;
ptrans2=exp(-thita*ca(1,2)/ ca_mean)/ca_sum;
X=[0,0];
X(1,1)=Q*ptrans1*beita+y(1,1)*(1-beita);
X(1,2)=Q*ptrans2*beita+y(1,2)*(1-beita);
 end

% ��ͼ
function fastplot(x,thita,t1,du_s1,ff1,ff2,m,figure_num,erfa,beita)
% ff1=ca2;
% ff2=X2;
% m=100;
% figure_num=x;
% erfa=0.2;
% beita=0.6;
n=m+40;
% x1=zeros(n,2);
% x2=zeros(n,2);
for j=-39:0
    i=40+j;
    x1(i,1)=t1(1,1);
    x1(i,2)=t1(1,2);
    x2(i,1)=du_s1(1,1);
    x2(i,2)=du_s1(1,2);
end
for j=1:m
    i=40+j;
    x1(i,1)=ff1{j}(1,1);
    x1(i,2)=ff1{j}(1,2);
    x2(i,1)=ff2{j}(1,1);
    x2(i,2)=ff2{j}(1,2);
end
figure
subplot(2,1,1)
plot(2,1,-39:m,x1(:,1),'r-*',-39:m,x1(:,2),'g-+');
str1=['��',num2str(x),'-',num2str(figure_num),'(��=',num2str(thita),')�η�����:�г�ʱ����α仯.','��=',num2str(erfa),',��=',num2str(beita)];
title(str1);
subplot(2,1,2)
plot(2,2,-39:m,x2(:,1),'r-*',-39:m,x2(:,2),'g-+');
str2=['��',num2str(x),'-',num2str(figure_num),'(��=',num2str(thita),')�η�����:������α仯.','��=',num2str(erfa),',��=',num2str(beita)];
title(str2);
end
%% ˫·��ƽ�����
function [s,t] = fsolve_flow(t01,t02,c1,c2,q,derta)
%UEƽ�������
syms x1 x2
t1=t01*(1+0.15*(x1/c1)^4)-derta;
t2=t02*(1+0.15*(x2/c2)^4);
% t1=t01*(1+(x1/c1))-10;
% t2=t02*(1+(x2/c2));
eq1=t1-t2;
eq2=x1+x2-q;
s=solve(eq1,eq2,x1,x2);
digits(6); %��������λ��
r=vpa([s.x1,s.x2]);
m=size(r,1);
j=1;
s=zeros(1,2);
for i=1:m
   if  r(i,1) >=0&&isreal(r(i,1))
       if r(i,2) >=0&&isreal(r(i,2))
           s(1,1)=r(i,1) ;
           s(1,2)=r(i,2) ;
       end
   end
end
t(1,1)=t01*(1+0.15*(s(1,1)/c1)^4);
t(1,2)=t02*(1+0.15*(s(1,2)/c2)^4);
end
function [s,t] = sue_fsolve_flow(t01,t02,c1,c2,q,sita)
%SUEƽ�������
syms x1 x2
t1=t01*(1+0.15*(x1/c1)^4);
t2=t02*(1+0.15*(x2/c2)^4);
% t1=t01*(1+(x1/c1));
% t2=t02*(1+(x2/c2));
eq1=sita*(t2-t1)-(log(x1)-log(x2));
eq2=x1+x2-q;
s=vpasolve(eq1,eq2,x1,x2);
% digits(6); %��������λ��
r=[s.x1,s.x2];
m=size(r,1);
j=1;
s=zeros(1,2);
for i=1:m
   if  r(i,1) >=0&&isreal(r(i,1))
       if r(i,2) >=0&&isreal(r(i,2))
           s(1,1)=r(i,1) ;
           s(1,2)=r(i,2) ;
       end
   end
end
t(1,1)=t01*(1+0.15*(s(1,1)/c1)^4);
t(1,2)=t02*(1+0.15*(s(1,2)/c2)^4);
end