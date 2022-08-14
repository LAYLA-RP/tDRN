%% MSA迭代次数k的倒数作为迭代权重，相比_k修改了图的显示，增加function G0_Plot(G_temp),修改function G_FlowPlot(G_temp)，修改%坐标改动
%增加od流量可视化，(G=G_append_OD(G0)用于改进G0)没用上，直接改node和line更方便
% 增加分车型
%% 路网基础数据输入
clear;clc
%%  step1： 路网基础数据录入
% GIS表格excel数据导入并以有向图G的格式存储;
[EdgeTable,NodeTable]=topotable();%拓扑结构基础数据导入；% load topo_basic  %%之前手动编的数据
% 生成基础有向路网
G0=digraph(EdgeTable,NodeTable);%node里面有name列的话，可以显示node的name% G=digraph(P);
% G0=G_append_OD(G0);%增加OD表示，没有必要，只要原始地图加上OD节点就可以
%%  step2：交通需求关于时段t的 OD矩阵录入.手动录入是必须保证原有表格node的排列顺序
[OD_input_p,OD_input_f,OD_struct] = OD_matr(2005,10,NodeTable,380,930);%缺乏数据,采用自动生成为了继续测试实验结果
% OD_input=OD_struct(1,1).volume;
% [mon_t]=daycal(year);
% for i=1:12
%     [OD_struct,OD_input_p,OD_pair] = OD_matr(year,mon_t(i,1));
% end
%% step4：养护日志相关数据录入，同上，完成表格并读取
%交互输入：施工交通组织方案编号及对应参数
%计划施工时间及对应参数：封道位置（桩号）/车道封闭方案参数
%输出：养护时序表→路网更新时序表
% 信号向量matrx_mon
% matrx_mon=timeseries([0:11],[0:11]);
%net=xlsread('maintain.xlsx','A2:C100');  
%对应时刻调用一次拓扑重构关键位置函数
% T=3;%PCI更新查询周期/月
%多次（无时段交叉）养护重分配
%t:维持的时长;programme:养护方案(0/1/2分别表示未进行养护/方案A养护/方案B养护);paremater:养护方案要提供的参数
[time_event,parameter]=xlsread('maintain.xlsx','A2:D100'); 
%% 初始状态
erfa1=0;%出行成本预测系数,该函数可以分析erfa对迭代效率是否能增加
theta_input=1.6;
beita1=1;epsilon = 0.0001;%设置初始误差及迭代次数
[H0{1}] =MSA_G(theta_input,G0,OD_input_p,erfa1,beita1,epsilon);
%初始状态图
% G0_Plot(G0)
% G_FlowPlot1(H0)
save G_basic0_2%存储变量方便调试
%% 单次路网更新重分配
% G_update=update_topu_meantime(time_event,parameter,G0);
% G_temp=MSA_G(theta_input,G_update,OD_input_p);%拓扑结构更新引起的交通重分配
% G_flow=flow_lane(time_event,G_temp);%time_event:t/programme(0表示没有养护方案/1、2分别表示两种方案)

%% 参数分析
clear;clc
load G_basic0_2.mat
x=1;%算例序号，用于自己分析时识别图
theta_input=1.6;
erfa=0.8;%出行成本预测系数,该函数可以分析erfa对迭代效率是否能增加
beita=0.8;epsilon = 0.0001;%设置初始误差及迭代次数
[G_LaneProperty,paramater_plus,H_MSA_Property]=G_dup(G0,time_event,parameter,theta_input,OD_input_p,erfa,beita,epsilon );
save G_basic_7

%% 路网历史流量数据统一划分路段并计算积累流量
clear;clc;
load G_basic_7.mat
[G_Flow,G_Property_Separator]=UpdateFlow_LinkSeparator3(G_LaneProperty,paramater_plus,time_event);%不同时段路网路段划分统一，并按需取出部分属性用于后面计算(Xi，PCI)
G_FlowSum=G_FlowSum_MSA(12,G_Flow,time_event);%根据PCI更新周期，计算相应的路段流量总和，eg.PCI更新周期为3m，那么FlowSum计算周期也为3m
save G_basic_8

%% 结果作图
clear;clc;
load G_basic_8.mat
%% 有向图可视化
G0_Plot(G0)
%原始结构静态分配结果
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
% G_FlowPlot2(G_LaneProperty)%内外侧车道同时显示
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
%% 路网数据更新
function [G_LaneProperty,paramater_plus,H_MSA_Property]=G_dup(G0,time_event,parameter,theta_input,OD_input_p,erfa1,beita1,epsilon)%路网拓扑更新-流量分配
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
        H_TopuUpdate{i}=G_temp;%拓扑变更的存储
        H_MSA_Property{i}=MSA_G(theta_input,H_TopuUpdate{i},OD_input_p,erfa1,beita1,epsilon);%拓扑结构更新引起的交通重分配，只用于计算，属性不是最后结果
        G_LaneProperty{i}=H_MSA_Property{i};%路网流量属性更新
    elseif time_event_temp(2)==1
        G_save_temp=GraphChange_A(G_temp,parameter_temp);%[18,19,300,500]方案A：路段局部封闭，中央分隔带打开分流到对向车;
        H_TopuUpdate{i}=G_save_temp{1};
        paramater_plus{i}=G_save_temp{2};
        H_MSA_Property{i}=MSA_G(theta_input,G_save_temp{1},OD_input_p,erfa1,beita1,epsilon);%拓扑结构更新引起的交通重分配
        G_LaneProperty{i}=G_LaneProperty_Change_A(H_MSA_Property{i},G_temp,parameter_temp);%[18,19,300,500]方案A
    elseif time_event_temp(2)==2
         parameter_temp_node=str2num(parameter{i,2});
         G_save_temp=GraphChange_B(G_temp,parameter_temp,parameter_temp_node(1,1));%[27,18;34,29],28方案B：匝道封闭%方案A：路段局部封闭，中央分隔带打开分流到对向车;
         H_TopuUpdate{i}=G_save_temp{2};
         H_MSA_Property{i}=MSA_G(theta_input,G_save_temp{2},OD_input_p,erfa1,beita1,epsilon);%拓扑结构更新引起的交通重分配
         G_LaneProperty{i}=G_LaneProperty_Change_B(H_MSA_Property{i},G_temp,G_save_temp{4},G_save_temp{3});%方案B
    end
end
end

%% 可视化function
%G0绘制，标出OD点
function G0_Plot(G_temp)
G=G_temp;
XData0=G.Nodes.XDate;
YData0=G.Nodes.YDate;
nn=length(G.Nodes.YoN);
j=0;
for i=1:nn
    if G.Nodes.YoN(i)==1
        j=j+1;
        hightlight_node(j)=i;%hightlight node提取
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
highlight(P,hightlight_node,'NodeColor','k')%突出显示
str0=['初始路网结构'];
title(str0)
P.XData=XData0';
P.YData=YData0';
% colorbar
end      

%车道pcu在路网的分布G_FlowPlot
function G_FlowPlot02(G_temp)%G_FlowPlot2基础上去坐标
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
            hightlight_node(jj)=j;%hightlight node提取
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
            str1=['第',num2str(i),'次更新路网外侧车道路段小时pcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['第',num2str(i),'次更新路网内侧车道路段小时pcu'];
        end
            title(str1)
%             P{j}.XData=XData0';
%             P{j}.YData=YData0';
            colorbar
    end
end
end
function G_FlowPlot1(G_temp)%页面上显示1个车道
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
            hightlight_node(jj)=j;%hightlight node提取
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
            str1=['第',num2str(i),'次更新路网外侧车道路段小时pcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['第',num2str(i),'次更新路网内侧车道路段小时pcu'];
        end
            title(str1)
            P{j}.XData=XData0';
            P{j}.YData=YData0';
            colorbar
    end
end
end
function G_FlowPlot01(G_temp)%G_FlowPlot1基础上去坐标
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
            hightlight_node(jj)=j;%hightlight node提取
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
            str1=['第',num2str(i),'次更新路网外侧车道路段小时pcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['第',num2str(i),'次更新路网内侧车道路段小时pcu'];
        end
            title(str1)
%             P{j}.XData=XData0';
%             P{j}.YData=YData0';
            colorbar
    end
end
end
function G_FlowPlot2(G_temp)%布局：页面上显示2个车道
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
            hightlight_node(jj)=j;%hightlight node提取
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
            str1=['第',num2str(i),'次更新路网外侧车道路段小时pcu'];
        else
            P{j}.EdgeLabel=G.Edges.X2i;
            str1=['第',num2str(i),'次更新路网内侧车道路段小时pcu'];
        end
            title(str1)
            P{j}.XData=XData0';
            P{j}.YData=YData0';
            colorbar
    end
end
end
%断面=车道和
function G_FlowPlot_T1(G_temp)%页面上显示n个车道和
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
            hightlight_node(jj)=j;%hightlight node提取
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
        str1=['第',num2str(i),'阶段路网小时pcu'];
        title(str1)
        P.XData=XData0';
        P.YData=YData0';
        colorbar
         colormap spring;  %colormap还可以设置为jet，hsv，hot，spring，summer，autumn，winter，gray，bone，copper，pink，lines。
        
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
        NODEID_od(j)=(i);%用于后面构造Edges
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
%有向图
G2=digraph(G_Edges,G_Nodes);
G1=digraph(G_Edges_New,G_Nodes_New);
end%增加OD表示

%车道流量逐次变化图
function EdgesPlot(G_FlowD,plot_Edges)
plot_vexnum=size(plot_Edges,1);%边数
s=fix(plot_vexnum/10);%边数太多，10次放一个图
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
            G_vexnum=length(G.Edges.EndNodes);%边数
            %查询位置
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

%% 初始function
function [EdgeTable,NodeTable] = topotable()
%
%   
%预处理line和node表格，形成topo需要的邻接关系表格EdgeTable和NodeTable 
line=readtable('line.xls');
% t0 = line.Shape_Length./line.Speed_d*(3600/1000);
% Cap=line.Cap_d.*cap_reduction(Si);
node=readtable('node.xls');
EdgeTable_a =  table([line.FNODE_,line.TNODE_],line.Cap_d1,line.Cap_d2,line.Speed_d,line.Shape_Length,line.PCI1,line.PCI2,'VariableNames',{'EndNodes','Cap_d1','Cap_d2','Speed_d','Distance','PCI1','PCI2'});
EdgeTable_b =  table([line.TNODE_,line.FNODE_],line.Cap_d1,line.Cap_d2,line.Speed_d,line.Shape_Length,line.PCI1,line.PCI2,'VariableNames',{'EndNodes','Cap_d1','Cap_d2','Speed_d','Distance','PCI1','PCI2'});
EdgeTable=[EdgeTable_a ;EdgeTable_b];%单向图数据变为双向图数据
NodeTable =  table(node.FID_node,node.ID,node.XDate,node.YDate,'VariableNames',{'NODEID','YoN','XDate','YDate'});

end
function[OD_input_p,OD_input_f, OD_struct] = OD_matr(year,mon,NodeTable,p,f)
%% 通过提取od点，形成od_pair表格，并提示完成od数据输入
k = find(NodeTable.YoN==1);
xlswrite('OD.xlsx',k','OD_input_f','B1');
xlswrite('OD.xlsx',k,'OD_input_f','	A2');
xlswrite('OD.xlsx',k','OD_input_p','B1');
xlswrite('OD.xlsx',k,'OD_input_p','	A2');
% disp('请按要求完成OD.xlsx中路段表格OD');
% disp('按任意键继续...');
% pause;
% OD_input_p=xlsread('OD.xlsx','OD_input_p','B2:ZZ10000');
% OD_input_f=xlsread('OD.xlsx','OD_input_f','B2:ZZ10000');

%% 为了完成实验，用自动生成随机od代替实际历史数据输入，并存入od_pair表格
t1=8;t2=16;%t1,t2分别为高峰期/平峰期的持续时长（/h）
OD_cell=cell(2,12);
A=cell(1,12);
B=cell(1,12);
OD_num=length(k);
OD_Node=zeros(OD_num,1);
for i = 1:OD_num
    OD_Node(i,1)=NodeTable.NODEID(k(i));
end
mon_i=[1,1,1.5,1.5,1,1,1,2,2,1,5,1,1];%不同月份，需求的变化
for i=1:12
    A{1,i}=round(rand(OD_num)+mon_i(i)*fix(mod(p,100)/10)*10+fix(p/100)*100);%随机生成OD代替io2OD().%提取十位数和百位数fix(mod(p,100)/10)+fix(p/100)*100)
    OD_cell{1,i}=A{1,i}-diag(diag(A{1,i}));%随机生成OD代替io2OD()平峰期(16)
    B{1,i}=round(rand(OD_num)+mon_i(i)*fix(mod(f,100)/10)*10+fix(f/100)*100);%随机生成OD代替io2OD()
    OD_cell{2,i}=B{1,i}-diag(diag(B{1,i}));%随机生成OD代替io2OD()高峰期(8)
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

%% 路网更新function
%% G_Topu function
%同一时段不同路段路网更新
function G_update=UpdateTopu_meantime(time_event,parameter,G0)
n=size(time_event,1);
G_temp=cell(n+1,1);
G_temp{1}=G0;
for i=1:n
     parameter_temp=str2num(parameter{i});
     if time_event(i,3)==1
         if time_event(i,2)==1
             G_temp{i+1}=GraphChange_A(G_temp{i},parameter_temp);%[18,19,300,500]方案A：路段局部封闭，中央分隔带打开分流到对向车;
         elseif time_event(i,2)==2
             G_temp{i+1}=GraphChange_B(G_temp{i},parameter_temp);%[27,18;34,29],28方案B：匝道封闭%方案A：路段局部封闭，中央分隔带打开分流到对向车;
         end
         elseif time_event(i,3)==0%方案AB需要修改为相应的反方案
             if time_event(i,2)==1
                 G_temp{i+1}=GraphChange_A(G_temp{i},parameter_temp);%[18,19,300,500]方案A：路段局部封闭，中央分隔带打开分流到对向车;
             elseif time_event(i,2)==2
                 G_temp{i+1}=GraphChange_B(G_temp{i},parameter_temp);%[27,18;34,29],28方案B：匝道封闭%方案A：路段局部封闭，中央分隔带打开分流到对向车;
             end
     end
end
G_update=G_temp{n+1};
end
%一个时段只有一个点进行养护,路网拓扑更新
%方案A
function G_save=GraphChange_A(G,parameter_temp)%相比main.k坐标改动+新增坐标平移100有利显示
% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
% node1=1;node2=2;r
G_origin=G;%原始的G
node1=parameter_temp(1);node2=parameter_temp(2);fd=parameter_temp(3);length=parameter_temp(4);
%% 取出原始G中的各个属性
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
%% 最大节点数，以便增加新节点
MaxNodeID=max(G_Nodes.NODEID);
%% 找出需要修改节点的位置
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node1&&G_Edges_EndNodes(i,2)==node2
        PosisionNumber=i;
    end
end
%对向车道对应edge位置
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node2&&G_Edges_EndNodes(i,2)==node1
        PosisionNumber_r=i;
    end
end
%% 把要修改的部分的前面和后面都取出来
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
%% 对新有向图的Nodes属性添加两个新节点，MaxNodeID+1,MaxNodeID+2
%目标路段节点坐标
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
%新增节点坐标
a=[XDate_node1,YDate_node1];
b=[XDate_node2,YDate_node2];
L_node12=norm(b-a);
XDate_node3= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
YDate_node3= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
XDate_node4= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
YDate_node4= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
%新增坐标平移
dertaX=100;
dertaY=dertaX*XDate_node1/YDate_node2;
%新有向图
G_Nodes_New=table;
G_Nodes_NODEID_New=[G_Nodes_NODEID;MaxNodeID+1;MaxNodeID+2];
G_Nodes_YoN_New=[G_Nodes_YoN;0;0];%新增的YoN设定为0
G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node3+dertaX;XDate_node4+dertaX];%新增XDate
G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node3+dertaY;YDate_node4+dertaY];%新增YDate
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% 对新有向图的Edges属性插入新节点，使得一段变成三段
G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),MaxNodeID+1;
    MaxNodeID+1,MaxNodeID+2;
    MaxNodeID+2,G_Edges_EndNodes(PosisionNumber,2)];
% G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber)];
G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
    G_Edges_Cap_d2(PosisionNumber_r);
    G_Edges_Cap_d1(PosisionNumber)];%考虑是否需要根据情境折减
G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
    0;
    G_Edges_Cap_d2(PosisionNumber)];%考虑是否需要根据情境折减
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
%     G_Edges_Num_lane(PosisionNumber)]*0.5;%车道减半
%% 将新增加的节点与之前没变化的拼接起来
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
%% 将新属性赋予到新有向图的Edges属性中
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
%% 构造新的有向图
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% 提取新增的节点
add_node1=max(G_new.Nodes.NODEID)-1;
r_parameter=[parameter_temp,add_node1];
%% 把新旧有向图放到元胞数组中
% G1{1}=G_origin;
G1{1}=G_new;
G1{2}= r_parameter;
G_save=G1;

end
function G_save=GraphChange_AA(G,parameter_temp)%相比A(代替A)，增加了对向分段（如果要使用，应该参照G_LaneProperty_renew_A修改G_LaneProperty_Change_A，因为19→已经在A_r中被删除
G_save_temp=GraphChange_A(G,parameter_temp);
G_new=GraphChange_A_r(G_save_temp{1},parameter_temp);
G_save={G_new,G_save_temp{2}};
end
function G_new=GraphChange_A_r(G,parameter_temp)%相比main.k坐标改动%在GraphChange_A基础上增加对向19→18分段
% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
% node1=1;node2=2;r
G_origin=G;%原始的G
node1=parameter_temp(2);node2=parameter_temp(1);fd=parameter_temp(4);length=parameter_temp(3);
%% 取出原始G中的各个属性
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
%% 最大节点数，以便增加新节点
MaxNodeID=max(G_Nodes.NODEID);
%% 找出需要修改节点的位置
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node1&&G_Edges_EndNodes(i,2)==node2
        PosisionNumber=i;
    end
end
%% 把要修改的部分的前面和后面都取出来
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
%% 对新有向图的Nodes属性添加两个新节点，MaxNodeID+1,MaxNodeID+2
%目标路段节点坐标
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
%新增节点坐标
a=[XDate_node1,YDate_node1];
b=[XDate_node2,YDate_node2];
L_node12=norm(b-a);
XDate_node3= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
YDate_node3= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
XDate_node4= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
YDate_node4= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
%新有向图
G_Nodes_New=table;
G_Nodes_NODEID_New=[G_Nodes_NODEID;MaxNodeID+1;MaxNodeID+2];
G_Nodes_YoN_New=[G_Nodes_YoN;0;0];%新增的YoN设定为0
G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node3;XDate_node4];%新增XDate
G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node3;YDate_node4];%新增YDate
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% 对新有向图的Edges属性插入新节点，使得一段变成三段
G_Edges_EndNodes_Append=[G_Edges_EndNodes(PosisionNumber,1),MaxNodeID+1;
    MaxNodeID+1,MaxNodeID+2;
    MaxNodeID+2,G_Edges_EndNodes(PosisionNumber,2)];
% G_Edges_Weight_Append=[G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber);
%     G_Edges_Weight(PosisionNumber)];
G_Edges_Cap_d1_Append=[G_Edges_Cap_d1(PosisionNumber);
    G_Edges_Cap_d2(PosisionNumber);
    G_Edges_Cap_d1(PosisionNumber)];%考虑是否需要根据情境折减
G_Edges_Cap_d2_Append=[G_Edges_Cap_d2(PosisionNumber);
    0;
    G_Edges_Cap_d2(PosisionNumber)];%考虑是否需要根据情境折减
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
%     G_Edges_Num_lane(PosisionNumber)]*0.5;%车道减半
%% 将新增加的节点与之前没变化的拼接起来
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
%% 将新属性赋予到新有向图的Edges属性中
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
%% 构造新的有向图
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% 提取新增的节点
add_node1=max(G_new.Nodes.NODEID)-1;
r_parameter=[parameter_temp,add_node1];
%% 把新旧有向图放到元胞数组中
% G1{1}=G_origin;
G1{1}=G_new;
G1{2}= r_parameter;
G_save=G1;

end

%方案B
function G_save=GraphChange_B(G,node1_3,node2)
%node1，node2和node3是对应路段是有方向性的，node1→node2→node3方向为前进方向
%node1_3，第一列是node1对应的数，第二列是node3对应的数，例如node1_3=[2,3;15,4]，即为封闭2→3,15→4两段
%fd,length分别为工作区至node1的距离和工作区长度
% node1=1;node2=2;r
G_origin=G;%原始的G
%% 取出原始G中的各个属性
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
%% 检索与node2相连的节点有多少个,即需要增加的节点数量node_num
cnt_node=[];
cnt_edges=[];
for i=1:size(G_Edges_EndNodes,1)
    node_nargin=G_Edges_EndNodes(i,:)-node2;
    if node_nargin(1)*node_nargin(2)==0%至少有一个为0，说明这一段中有node2
        cnt_node=[cnt_node;node_nargin(1)+node_nargin(2)];
        cnt_edges=[cnt_edges;G_Edges_EndNodes(i,:)];
    end
end
% node_num=size(unique(cnt_node),1);%不同元素的个数
node_num=size(cnt_node,1);%不同元素的个数
%% 最大节点数，以便增加新节点
MaxNodeID=max(G_Nodes.NODEID);
%% 新建node_num个节点
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
    G_Nodes_YoN_New=[G_Nodes_YoN_New;0];%新增的YoN设定为0
    G_Nodes_XDate_New=[G_Nodes_XDate_New;XDate];
    G_Nodes_YDate_New=[G_Nodes_YDate_New;YDate];
end
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% 对新有向图的Edges,把包含有node2和不包含node2的Edge分开
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
%% 把和node2相连的节点换成新节点，并和之前的合并
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

%% 增加拆分出来之后节点的相互连接，这些点之间的权重和距离设为0，node1→node2、node2→node3对应的不连接

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
%% 将三个部分放一起
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
%% 将新属性赋予到新有向图的Edges属性中
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
%% 构造新的有向图
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% 把新旧有向图放到元胞数组中
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
%方案A下车道流量及其他属性重分配
function G_save=G_LaneProperty_Change_A(G_MSA,G,parameter_temp)%坐标改动+新增坐标平移100
%function [G_new,r_parameter]=GraphChange_A(G,r_parameter)
% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
% node1=1;node2=2;r
G_origin=G_MSA;%原始的G
node1=parameter_temp(1);node2=parameter_temp(2);fd=parameter_temp(3);length=parameter_temp(4);
%% 取出原始G中的各个属性
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
%% 最大节点数，以便增加新节点
MaxNodeID=max(G_Nodes.NODEID);
%% 找出需要修改节点的位置
%提取封道对向车道的节点[19→18]
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==node2&&G_Edges_EndNodes(i,2)==node1
        PosisionNumber=i;
    end
end
%提取封道的节点，即G_MSA的最后两个节点[47→48]
for i=1:size(G_Edges_EndNodes,1)
    if G_Edges_EndNodes(i,1)==(MaxNodeID-1)&&G_Edges_EndNodes(i,2)==MaxNodeID
        PosisionNumber_1=i;
    end
end
%% 把要修改的部分的前面和后面都取出来
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
%% 对新有向图的Nodes属性添加两个新节点，MaxNodeID+1,MaxNodeID+2
%目标路段节点坐标
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
%新增节点坐标
a=[XDate_node1,YDate_node1];
b=[XDate_node2,YDate_node2];
L_node12=norm(b-a);
XDate_node3= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*fd;
YDate_node3= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*fd;
XDate_node4= XDate_node1 + (XDate_node2 -XDate_node1)/L_node12*(fd+length);
YDate_node4= YDate_node1 + (YDate_node2 -YDate_node1)/L_node12*(fd+length);
%新增坐标平移
dertaX=100;
dertaY=dertaX*XDate_node1/YDate_node2;
%新有向图
G_Nodes_New=table;
G_Nodes_NODEID_New=[G_Nodes_NODEID;MaxNodeID+1;MaxNodeID+2];
G_Nodes_YoN_New=[G_Nodes_YoN;0;0];%新增的YoN设定为0
G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node4-dertaX;XDate_node3-dertaX];%新增XDate
G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node4-dertaY;YDate_node3-dertaY];%新增YDate
G_Nodes_New.NODEID=G_Nodes_NODEID_New;
G_Nodes_New.YoN=G_Nodes_YoN_New;
G_Nodes_New.XDate=G_Nodes_XDate_New;
G_Nodes_New.YDate=G_Nodes_YDate_New;
%% 对新有向图的Edges属性插入新节点，使得一段变成三段
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
%     G_Edges_Num_lane(PosisionNumber)];%车道数不变
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

%% 将新增加的节点与之前没变化的拼接起来
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

%% 局部修正 
%提取封道前封道位置的PCI和Cap[18→19]
for i=1:size(G.Edges.EndNodes,1)
    if G.Edges.EndNodes(i,1)==node1&&G.Edges.EndNodes(i,2)==node2
        PosisionNumber_3=i;
    end
end
%封道处流量/速度修正为0[47→48]
for i=1:size(G_Edges_EndNodes_New,1)
    if G_Edges_EndNodes_New(i,1)==(MaxNodeID-1)&&G_Edges_EndNodes_New(i,2)==MaxNodeID
        PosisionNumber_2=i;
    end
end
G_Edges_X1i_New(PosisionNumber_2)=0.001;
G_Edges_X2i_New(PosisionNumber_2)=0.001;
G_Edges_vi1_New(PosisionNumber_2)=0.001;
G_Edges_vi2_New(PosisionNumber_2)=0.001;
%Cap和PCI恢复成封道之前的数据
G_Edges_Cap_d1_New(PosisionNumber_2)=G.Edges.Cap_d1(PosisionNumber_3);
G_Edges_Cap_d2_New(PosisionNumber_2)=G.Edges.Cap_d2(PosisionNumber_3);
G_Edges_PCI1_New(PosisionNumber_2)=G.Edges.PCI1(PosisionNumber_3);
G_Edges_PCI2_New(PosisionNumber_2)=G.Edges.PCI2(PosisionNumber_3);


%% 将新属性赋予到新有向图的Edges属性中
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
%% 构造新的有向图
G_new=digraph(G_Edges_New,G_Nodes_New);
% plot(G_new);
%% 把新旧有向图放到元胞数组中
G1{1}=G_origin;
G1{2}=G_new;
G_save=G1{2};
end

function G_save=G_LaneProperty_renew_A(H_day2day_Property,G_origin,parameter_temp)
%主要是47→48和49→50处的PCI和车道Cap要恢复成交通分配之前的值即18-19和19-18；流量加载回实际车道上，便于流量影响PCI
%function [G_new,r_parameter]=GraphChange_A(G,r_parameter)
% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
% node1=1;node2=2;r
G1{1}=G_origin;%原始的G
G1{2}=H_day2day_Property;

% 取出原始G中的各个属性
G_Edges{1}=G1{1}.Edges;
G_Nodes{1}=G1{1}.Nodes;
G_Edges_EndNodes{1}=G_Edges{1}.EndNodes;

G_Edges{2}=G1{2}.Edges;
G_Nodes{2}=G1{2}.Nodes;
G_Edges_EndNodes{2}=G_Edges{2}.EndNodes;

%% 找出需要修改节点的位置
MaxNodeID=max(G_Nodes{2}.NODEID);
%提取封道对向车道的节点[47→48]
for i=1:size(G_Edges_EndNodes{2},1)
    if G_Edges_EndNodes{2}(i,1)==MaxNodeID-3&&G_Edges_EndNodes{2}(i,2)==MaxNodeID-2
        PosisionNumber_1=i;
    end
end
%提取封道的节点，对向[49→50]
for i=1:size(G_Edges_EndNodes{2},1)
    if G_Edges_EndNodes{2}(i,1)==(MaxNodeID-1)&&G_Edges_EndNodes{2}(i,2)==MaxNodeID
        PosisionNumber_2=i;
    end
end

node1=parameter_temp(1);node2=parameter_temp(2);fd=parameter_temp(3);length=parameter_temp(4);
%提取封道对向车道的节点[18→19]
for i=1:size(G_Edges_EndNodes{1},1)
    if G_Edges_EndNodes{1}(i,1)==node1&&G_Edges_EndNodes{1}(i,2)==node2
        PosisionNumber=i;
    end
end
%提取封道的节点，对向[19→18]
for i=1:size(G_Edges_EndNodes{1},1)
    if G_Edges_EndNodes{1}(i,1)==node2&&G_Edges_EndNodes{1}(i,2)==node1
        PosisionNumber_r=i;
    end
end

%% 对新有向图{3}的Edges部分属性恢复为{1}的属性，除了{2}的属性，主要是把流量加载到了实体车道上
%47-48部分
G_Edges{2}.Cap_d1(PosisionNumber_1)=G_Edges{1}.Cap_d1(PosisionNumber);
G_Edges{2}.Cap_d2(PosisionNumber_1)=G_Edges{1}.Cap_d2(PosisionNumber);
G_Edges{2}.Speed_d(PosisionNumber_1)=G_Edges{2}.Speed_d(PosisionNumber);
G_Edges{2}.PCI1(PosisionNumber_1)=G_Edges{1}.PCI1(PosisionNumber);
G_Edges{2}.PCI2(PosisionNumber_1)=G_Edges{1}.PCI2(PosisionNumber);
G_Edges{2}.X1i(PosisionNumber_1)=0.001;
G_Edges{2}.X2i(PosisionNumber_1)=0.001;
G_Edges{2}.vi(PosisionNumber_1)=0.001;
% G_Edges{2}.Ti(PosisionNumber_1)=G_Edges{1}.Ti(PosisionNumber);%T在分配的时候重新计算

%49-50部分
% G_Edges{2}.Cap_d1(PosisionNumber_2)=G_Edges{1}.Cap_d1(PosisionNumber_r);
G_Edges{2}.Cap_d2(PosisionNumber_2)=G_Edges{1}.Cap_d2(PosisionNumber_r);
% G_Edges{3}.Speed_d(PosisionNumber_2)=G_Edges{2}.Speed_d(PosisionNumber_r);
% G_Edges{3}.PCI1(PosisionNumber_2)=G_Edges{1}.PCI1(PosisionNumber_r);
G_Edges{2}.PCI2(PosisionNumber_2)=G_Edges{1}.PCI2(PosisionNumber_r);
G_Edges{2}.X1i(PosisionNumber_2)=G_Edges{2}.X1i(PosisionNumber_2)+G_Edges{2}.X2i(PosisionNumber_2);
G_Edges{2}.X2i(PosisionNumber_2)=G_Edges{2}.X1i(PosisionNumber_1)+G_Edges{2}.X2i(PosisionNumber_1);
% G_Edges{2}.v1i(PosisionNumber_2)=G_Edges{2}.v1i(PosisionNumber_2)+G_Edges{2}.v1i(PosisionNumber_2);
% G_Edges{2}.v2i(PosisionNumber_2)=G_Edges{2}.v2i(PosisionNumber_1)+G_Edges{2}.v2i(PosisionNumber_1);
% G_Edges{2}.Ti(PosisionNumber_1)=G_Edges{1}.Ti(PosisionNumber);%T在分配的时候重新计算

%% 把新旧有向图放到元胞数组中
G_save=G1{2};
end

%方案B下车道流量及其他属性重分配
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
%流量相加辅助：流量相加前路网路段划分统一后，按需取出部分属性用于后面计算(包含UpdateFlow_LinkSeparator_A)
function [G_Flow,G_Property_Separator]=UpdateFlow_LinkSeparator(G_LaneProperty,paramater_plus,time_event)
%不同时段路网路段划分统一，并按需取出部分属性用于后面计算(Xi，PCI)
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
    G_Flow_Edges{i}=G_Property_Separator{i}.Edges(:,{'EndNodes','X1i','X2i','PCI1','PCI2'});%按需取出部分属性用于后面计算
    G_Flow_Nodes{i}=G_Property_Separator{i}.Nodes;
    G_Flow{i}=digraph(G_Flow_Edges{i},G_Flow_Nodes{i});
end
end
%UpdateFlow_LinkSeparator中function：流量相加前路网路段划分统一
function G_save=UpdateFlow_LinkSeparator_A(G_LaneProperty,paramater_r)%坐标改动
%流量计算完以后，对前面的路网进行分段，使得不同时段的路网可以直接相加得到路段流量

% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
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
    %% 对新有向图重复两次：1段变三段，增加两个节点
    G_Edges_Append_temp=cell(2,1);
     for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        % node1-node2距离L
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
        %% 对新有向图的Edges属性插入新节点，使得一段变成三段
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
        %% 将新增加的路段属性组成一张新的table_append
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
 %% 对新有向图的Nodes属性添加两个新节点
    G_Nodes_New=table;
    G_Nodes_NODEID_New=[G_Nodes_NODEID;node5;node5+1;node5+2;node5+3];    
    G_Nodes_YoN_New=[G_Nodes_YoN;0;0;0;0];%新增的YoN设定为0
    G_Nodes_XDate_New=[G_Nodes_XDate;XDate_node3(1);XDate_node4(1);XDate_node3(2);XDate_node4(2)];
    G_Nodes_YDate_New=[G_Nodes_YDate;YDate_node3(1);YDate_node4(1);YDate_node3(2);YDate_node4(2)];
    G_Nodes_New.NODEID=G_Nodes_NODEID_New;
    G_Nodes_New.YoN=G_Nodes_YoN_New;
    G_Nodes_New.XDate=G_Nodes_XDate_New;
    G_Nodes_New.YDate=G_Nodes_YDate_New;
     %% 将原有的G_Edges删除新增路段位置的原始路段，得到G_Edges_0（18→19,19→18）
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
    %% 两个方向的新增路段合并，共增加6个路段
    G_Edges_Add=[G_Edges_Append_temp{1};G_Edges_Append_temp{2}];
    %% 组成完整的路段table
    G_Edges_New_temp=cell(m,1);
    G_Edges_New_temp{ii}=[G_Edges_0;G_Edges_Add];
    %% 构造新的有向图
    G_new=cell(m,1);
    G_new{ii}=digraph(G_Edges_New_temp{ii},G_Nodes_New);
    % plot(G_new);
end
G_save=G_new{m};
end
%适用于周期末结果G_LaneProperty_T_v
function [G_Flow,G_FlowD_Separator]=UpdateFlow_LinkSeparator3(G_LaneProperty_T_v,paramater_plus,time_event)
%不同时段路网路段划分统一，并按需取出部分属性用于后面计算(Xi，PCI)
%相比2多了Cap_d等属性迁移
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
    G_Flow_Edges{i}=G_FlowD_Separator{i}.Edges(:,{'EndNodes','X1i','X2i','PCI1','PCI2'});%按需取出部分属性用于后面计算
    G_Flow_Nodes{i}=G_FlowD_Separator{i}.Nodes;
    G_Flow{i}=digraph(G_Flow_Edges{i},G_Flow_Nodes{i});
end
end
%UpdateFlow_LinkSeparator中function：流量相加前路网路段划分统一
function G_save=UpdateFlow_LinkSeparator_A4_MSA(G_LaneProperty,paramater_r)%eg38坐标改动dertaX
%流量计算完以后，对前面的路网进行分段，使得不同时段的路网可以直接相加得到路段流量

% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
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
    %% 对新有向图重复两次：1段变三段，增加两个节点
    G_Edges_Append_temp=cell(2,1);
    dertaX=[-40,80];%eg38坐标，出图美化，偏移
     for jj=1:2
        node1=paramater_r_temp(jj,1);
        node2=paramater_r_temp(jj,2);
        % node1-node2距离L
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
        
        %% 对新有向图的Edges属性插入新节点，使得一段变成三段
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
        %% 将新增加的路段属性组成一张新的table_append
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
 %% 对新有向图的Nodes属性添加2+2个新节点
         %新有向图
        G_Nodes_New=table;
        G_Nodes_NODEID_New=[G_Nodes_NODEID;node5;node5+1;node5+2;node5+3];    
        G_Nodes_YoN_New=[G_Nodes_YoN;0;0;0;0];%新增的YoN设定为0
        G_Nodes_XData_New=[G_Nodes_XData;XData_node3(1)-dertaX(1);XData_node4(1)-dertaX(1);XData_node3(2)-dertaX(2);XData_node4(2)-dertaX(2)];%eg38坐标
        G_Nodes_YData_New=[G_Nodes_YData;YData_node3(1)-dertaY(1);YData_node4(1)-dertaY(1);YData_node3(2)-dertaY(2);YData_node4(2)-dertaY(2)];%eg38坐标
        G_Nodes_New.NODEID=G_Nodes_NODEID_New;
        G_Nodes_New.YoN=G_Nodes_YoN_New;
        G_Nodes_New.XDate=G_Nodes_XData_New;
        G_Nodes_New.YDate=G_Nodes_YData_New;
     %% 将原有的G_Edges删除新增路段位置的原始路段，得到G_Edges_0（18→19,19→18）
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
    %% 两个方向的新增路段合并，共增加6个路段
    G_Edges_Add=[G_Edges_Append_temp{1};G_Edges_Append_temp{2}];
    %% 组成完整的路段table
    G_Edges_New_temp=cell(m,1);
    G_Edges_New_temp{ii}=[G_Edges_0;G_Edges_Add];
    %% 构造新的有向图
    G_new=cell(m,1);
    G_new{ii}=digraph(G_Edges_New_temp{ii},G_Nodes_New);
    % plot(G_new);
end
G_save=G_new{m};
end
%流量相加
function G_FlowSum_derta=G_FlowSum_MSA(T,G_Flow,time_event)
%根据需求，分段计算路段流量和
% 根据PCI更新周期T，计算相应的路段流量总和，eg.PCI更新周期为3m，那么FlowSum计算周期也为3m
%G_Property_Separator,time_event为对应的属性和时段列表time_event[i-1,i)=G_Property_Separator{i}
time=[0;time_event(:,1)];
n=size(time,1);
m=fix(time(n,1)/T);
u=30*24;
%计算整块矩形面积（流量和）
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
%分周期计算流量和
for j=1:m
    G_FlowSum{j}=G_Flow{1};
    t=T*j;
    %求算时刻t所在序列time的位置
    for ii=1:n
        if t>time(ii)&&t<=time(ii+1)
            possionnum_t=ii;
        end
    end
    %计算非整分以外的部分
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
%计算周期增长
G_FlowSum_derta{1}=G_FlowSum{1};
for j=2:m
    G_FlowSum_derta{j}=G_FlowSum{j};
    G_FlowSum_derta{j}.Edges.X1i=G_FlowSum{j}.Edges.X1i-G_FlowSum{j-1}.Edges.X1i;
    G_FlowSum_derta{j}.Edges.X2i=G_FlowSum{j}.Edges.X2i-G_FlowSum{j-1}.Edges.X2i;      
end
end
%% G_Flow function2
function G_FlowD_Separator=UpdateFlow_LinkSeparator2(G_FlowD,paramater_plus,time_event)
%不同时段路网路段划分统一，并按需取出部分属性用于后面计算(Xi，PCI)
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
%     G_Flow_Edges{i}=G_FlowD_Separator{i}.Edges(:,{'EndNodes','X1i','X2i','PCI1','PCI2'});%按需取出部分属性用于后面计算
%     G_Flow_Nodes{i}=G_FlowD_Separator{i}.Nodes;
%     G_Flow{i}=digraph(G_Flow_Edges{i},G_Flow_Nodes{i});
end
end
function G_save=UpdateFlow_LinkSeparator_A2(G_FlowDi,paramater_r)
%流量计算完以后，对前面的路网进行分段，使得不同时段的路网可以直接相加得到路段流量

% 修改节点的函数，输出G_save是一个元胞数组，G_save{1}是原始结构，G_save{2}是更新后的结构，使用时根据需要调用
%node1和node2是对应路段是有方向性的，node1→node2方向为前进方向
%fd,length分别为工作区至node1的距离和工作区长度
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
    
    %% 对新有向图重复两次：1段变三段，增加两个节点
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
        %% 对新有向图的Nodes属性添加两个新节点
        G_Nodes_New=table;
        G_Nodes_NODEID_New=[G_Nodes_NODEID;node5;node5+1;node5+2;node5+3];
        G_Nodes_YoN_New=[G_Nodes_YoN;0;0;0;0];%新增的YoN设定为0
        G_Nodes_New.NODEID=G_Nodes_NODEID_New;
        G_Nodes_New.YoN=G_Nodes_YoN_New;
        %% 对新有向图的Edges属性插入新节点，使得一段变成三段
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

        
        %% 将新增加的路段属性组成一张新的table_append
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
     %% 将原有的G_Edges删除新增路段位置的原始路段，得到G_Edges_0（18→19,19→18）
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
    %% 两个方向的新增路段合并，共增加6个路段
    G_Edges_Add=[G_Edges_Append_temp{1};G_Edges_Append_temp{2}];
    %% 组成完整的路段table
    G_Edges_New_temp=cell(n1,1);
    G_Edges_New_temp{ii}=[G_Edges_0;G_Edges_Add];
    %% 构造新的有向图
    G_new=cell(n1,1);
    G_new{ii}=digraph(G_Edges_New_temp{ii},G_Nodes_New);
    % plot(G_new);
end
G_save=G_new{n1};
end

%% 逐次仿真function
%关键1：行程时间需要逐次更新
%关键2：单次多路径分配采用dail分配得到Xi
%% 逐日仿真中单次计算函数（路网）,MSA变体，主要是重复进行dail计算
function G_save=day2day_mG(m,theta_input,G,OD_demand,erfa,beita,epsilon)
% G=H0;
% OD_demand=OD_input_p;
%% 网络参数
EdNodes=G.Edges.EndNodes;
G_Nodes=G.Nodes;
TNodes=EdNodes(:,1);
ENodes=EdNodes(:,2);
G_vexnum=length(EdNodes);%边数
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
%% 逐次计算，并判断是否达到稳定，什么时候达到稳定
for ii=1:m
    if ii==1
        erfa_temp=0;
        G_temp=G;%包含Xi及对应Wi
        preca_last=ones(G_vexnum,1);%因为erfa取了0，Wi_last取任何数结果一样，所以这里可以用常数1代替
        Xi_last=Xi_0;
        ca_last=cost_set(Wi0,Xi_last,Ci);
    else
        erfa_temp=erfa;
        G_temp=H{ii-1};%预留，如果要输出每日G，这里可以修改为H{ii}
        preca_last=preca(:,ii-1);
        ca_last=ca_temp(:,ii-1);
    end
    preca(:,ii)=ca_last.*(1-erfa_temp)+preca_last.*erfa_temp;%预测行程时间
    [count1,count2]=dial_G(theta_input,G_temp,OD_demand,preca(:,ii));
    Y1=zeros(G_vexnum,1);%还原至有向图格式
    for i=1:G_vexnum
        a=TNodes(i);
        b=ENodes(i);
        Y1(i,1)=round(count1(a,b));
    end
    %考虑β计算单次分配流量
    Xi_temp(:,ii)=Y1.*beita+Xi_last.*(1-beita);
    G_temp.Edges.X1i=round(Xi_temp(:,ii)./2);
    G_temp.Edges.X2i=round(Xi_temp(:,ii)./2);
    H{ii}=G_temp;
    %实际行程时间
    ca_temp(:,ii)=cost_set(Wi0,Xi_temp(:,ii),Ci);
    % 路段平均车速
    vi_temp(:,ii)=Di./ca_temp(:,ii)*3.6;
end

    %% 本周期最终状态，包括X和Speed，用于下一周期初始化
    G_temp2=G_temp;
    G_temp2.Edges.X1i=round(Xi_temp(:,m)./2);
    G_temp2.Edges.X2i=round(Xi_temp(:,m)./2);
    G_temp2.Edges.vi1=round(vi_temp(:,m));
    G_temp2.Edges.vi2=round(vi_temp(:,m)); 
    G1{1}=G_temp2;
    %% 逐日X及Speed存储
    X1i_temp=round(Xi_temp./2);
    X2i_temp=round(Xi_temp./2);
    vi1_temp=round(vi_temp);
    vi2_temp=round(vi_temp);   
    EdgeTable_2 =  table(EdNodes,X1i_temp,X2i_temp,vi1_temp,vi2_temp,'VariableNames',{'EndNodes','X1i','X2i','vi1','vi2'});
    G1{2}=digraph(EdgeTable_2,G_Nodes);
     %% 流量合计及平均车速存入新的有向图
    Xi_sum=sum(Xi_temp,2);%流量和
    vi_mu=sum(vi_temp,2)./m;%整个时段平均车速
    G_temp.Edges.X1i=round(Xi_sum./2);
    G_temp.Edges.X2i=round(Xi_sum./2);
    G_temp.Edges.vi1=vi_mu;
    G_temp.Edges.vi2=vi_mu;
    % G.Edges.Ti=Ti;
    G1{3}=G_temp;
    G_save=G1;
     %% step4 判断当前路段流量是否满足收敛条件
     error= sum(abs(Xi_temp(:,m)-Xi_temp(:,m-1)))/sum(Xi_temp(:,m));
     if error<epsilon
         fprintf('稳定,周期%d,误差error=%8.5f\n',m,error)
     else
         fprintf('不稳定,周期%d,误差error=%8.5f\n',m,error)
     end
 
end

%% MSA function
function [H] =MSA_G(theta_input,G,OD_demand,erfa,beita,epsilon)
%就是原来的MSA_G，用于初始状态计算
%step1：采用dail分配，dail中出行成本用t0代替，计算出X1；
%step2：用X更新出行成本，作为权重，采用0-1分配计算出Y2；
%step3:用加权迭代的方法，得到X2,即X2为X1和Y2-X1的加权和；
%step4：比较X2和X1，判断是否满足迭代精度要求，不满足则重复2.3.4step。
%该方法中只采用了一次dail计算，迭代原理为逐次加权分配
%erfa，出行时间更新的习惯性权重
% G=G0;
% OD_demand=OD_input_p;
error=0.3;%设置初始迭代误差
k=1;%迭代次数，倒数作为迭代权重
x_p=0.5;%设置双车道内车道的流量比例。
OD_Nodes=G.Nodes.NODEID(logical(G.Nodes.YoN));%OD时变矩阵的行列索引，与拓扑结构的编号对应；
G_ODnum=length(OD_Nodes);
topo_Nodes=G.Nodes.NODEID;%拓扑结构所有节点的编号
EdNodes=G.Edges.EndNodes;
TNodes=EdNodes(:,1);
ENodes=EdNodes(:,2);
G_vexnum=length(EdNodes);%边数
G_nodenum=length(topo_Nodes);%节点数
% Num_lane= G.Edges.Num_lane;
Ci=G.Edges.Cap_d1+G.Edges.Cap_d2;
Di=G.Edges.Distance;
Spi=G.Edges.Speed_d;
Wi0 =Di./Spi.*3.6 ;%仅与speed_d以及distance有关,需要更新，区别原来的Wi = G.Edges.Weight
% Xi=G.Edges.X1i+G.Edges.X2i;
%% 路网状态初始化
%dial分配
%初始化路段流量X1
[count1,count2]=dial_G(theta_input,G,OD_demand,Wi0);
X1=zeros(G_vexnum,1);%还原至有向图格式
for i=1:G_vexnum
    a=TNodes(i);
    b=ENodes(i);
    X1(i,1)=round(count1(a,b));
end
%初始化行程时间

time_temp= cost_set(Wi0,X1,Ci);
while error>epsilon
    %% step1：更新路段行驶时间time
    time_temp=precost_update(Wi0,X1,Ci,time_temp,erfa);
    %% step2: 按step1中计算出的路段时间和OD交通量，执行一次dail分配，得到一组附加交通量Y1
    [count1,count2]=dial_G(theta_input,G,OD_demand,time_temp);
    Y1=zeros(G_vexnum,1);
    for i=1:G_vexnum
        a=TNodes(i);
        b=ENodes(i);
        Y1(i,1)=round(count1(a,b));
    end
    %% step3 计算各路段当前交通量
    X2=X1+1/k.*(Y1-X1); %   0< phi<1 phi=0.%减少迭代次数，但是没有达到真正的SUE状态
%     X2=X1+beita.*(Y1-X1); %   可以达到均衡，但是计算次数太大
    %% step4 判断当路段流量是否满足收敛条件
    error= sum(abs(X1-X2))/sum(X1);
    %继续迭代
    X1=X2;
    k=k+1
end
 fprintf('迭代次数%d,误差error=%8.5f\n',k,error)

%% 计算路段平均速度
T1= cost_set(Wi0,X1,Ci);
V1=Di./T1*3.6;
%% 存入有向图表格
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

%MSA中function：dail/time_update/construct2
function [count1 count2]=dial_G(theta_input,G,OD_demand,W)
%求起点到终点最短路
% OD_demand,交通需求OD时变矩阵；
OD_Nodes=G.Nodes.NODEID(logical(G.Nodes.YoN));%OD时变矩阵的行列索引，与拓扑结构的编号对应；
topo_Nodes=G.Nodes.NODEID;%拓扑结构所有节点的编号
EdNodes=G.Edges.EndNodes;
TNodes=EdNodes(:,1);
ENodes=EdNodes(:,2);
G_nodenum=length(topo_Nodes);%节点数
G_ODnum=length(OD_Nodes);%OD节点数
% W = G.Edges.Distance./G.Edges.Speed_d.*3.6;
L=zeros(1,length(EdNodes));
% g=biograph(sparse(TNodes,ENodes,W));
% view(g);
%% 起点至所有节点的最短路集合r，所有节点至终点的最短路集合s
G_graph=G;
G_graph.Edges.Weight=W;
for i=1:G_nodenum
    for j=1:G_nodenum
        [path,dist] = shortestpath(G_graph,topo_Nodes(i),topo_Nodes(j));%bigraph形成的最短路工具箱
%         [dist,path]=graphshortestpath(linkweight,i,j);%biograph形成的最短路工具箱
        dis(topo_Nodes(i),topo_Nodes(j))=dist;
    end
end

%定义发点
%for i=1:9
Oi=cell(1,G_nodenum);
Di=cell(1,G_nodenum);
for i=1:G_nodenum
    k=topo_Nodes(i);
    m=find(TNodes==k);
    Oi{k}=[ENodes(m)];%节点标记
end
%定义收点
for i=1:G_nodenum
    k=topo_Nodes(i);
    m=find(ENodes==k);
    Di{k}=[TNodes(m)];%节点标记
end

%%正向流
count_positive0=sparse(TNodes,ENodes,L);
count_positive=spear2full_sq(count_positive0);
for i=1:G_ODnum
    origin=OD_Nodes(i);
    for j=1:G_ODnum
        if i~=j
            destination=OD_Nodes(j);
            r=dis(origin,:);%节点标记
            s=dis(:,destination)';
            Q=OD_demand(i,j);%位置标记
            Lsparse=likelihood_G(r,s,theta_input,W,EdNodes);
            Wi=weight2_G(r,Lsparse,Oi,Di,origin,EdNodes,destination);
            Xi=flow2_G(s,EdNodes,Oi,Di,Wi,destination,Q,origin);
            count_positive=count_positive+Xi;
        end
    end
end

%%反向流
count_negative0=sparse(TNodes,ENodes,L);
count_negative=spear2full_sq(count_negative0);
for i=1:G_ODnum
    origin=OD_Nodes(i);
    for j=1:G_ODnum
        if i~=j
        destination=OD_Nodes(j);
        r=dis(origin,:);
        s=dis(:,destination)';
        Q=OD_demand(i,j);%j为OD需求的序号，非节点号
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
%BPR函数
function [time] = cost_set(W0,X,Cap)%计算运行速度的时候，用于计算实际行驶时间
erfa=0.15;
berta=4;
time=W0.*(1+erfa.*(X./Cap).^berta);
time(isnan(time))=0;
end
function [time] = cost_set1(W0,X,Cap)%计算运行速度的时候，用于计算实际行驶时间
erfa=0.15;
berta=4;
Cap(Cap-X<=0)=0.001;
time=W0.*(1+erfa.*(X./Cap).^berta);
time(isnan(time))=0;
end
 %考虑习惯性选择的本次行程感知时间
 function  preca=precost_update(Wi0,X_last,Cap,preca_last,erfa)
 %经验性系数（习惯性系数）erfa=0，说明单独依靠前一次实际x，不考虑前一次预测的结果，如果之前处于平衡状态，这里也等于前一次实际的行程时间。
  ca_last = cost_set(Wi0,X_last,Cap);
  preca=ca_last.*(1-erfa)+preca_last.*erfa;
 end
%dail 中function:likehood/weight/flow
function Xi=flow2_G(s,EdNodes,Oi,Di,Wi,destination,Q,origin)
%Nodes，路网G路段邻接关系=G.Edges.EndNodes
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
%W为G中的路段阻抗权重weight=G.Edges.Weight
%Nodes=G.Edges.EndNodes，路段起终点编号
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
[a,b]=sort(r);%对r排序
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
%MSA辅助function
function possiablePaths = findPath(Graph, partialPath, destination, partialWeight)
% findPath按深度优先搜索所有可能的从partialPath出发到destination的路径，这些路径中不包含环路
% Graph: 路网图，非无穷或0表示两节点之间直接连通，矩阵值就为路网权值
% partialPath: 出发的路径，如果partialPath就一个数，表示这个就是起始点
% destination: 目标节点
% partialWeight: partialPath的权值，当partialPath为一个数时，partialWeight为0
pathLength = length(partialPath);
lastNode = partialPath(pathLength); %得到最后一个节点
nextNodes = find(0<Graph(lastNode,:) & Graph(lastNode,:)<inf); %根据Graph图得到最后一个节点的下一个节点
GLength = length(Graph);
possiablePaths = [];
if lastNode == destination
    % 如果lastNode与目标节点相等，则说明partialPath就是从其出发到目标节点的路径，结果只有这一个，直接返回
    possiablePaths = partialPath;
    possiablePaths(GLength + 1) = partialWeight;
    return;
elseif length( find( partialPath == destination ) ) ~= 0
    return;
end
%nextNodes中的数一定大于0,所以为了让nextNodes(i)去掉，先将其赋值为0
for i=1:length(nextNodes)
    if destination == nextNodes(i)
        %输出路径
        tmpPath = cat(2, partialPath, destination);      %串接成一条完整的路径
        tmpPath(GLength + 1) = partialWeight + Graph(lastNode, destination); %延长数组长度至GLength+1, 最后一个元素用于存放该路径的总路阻
        possiablePaths( length(possiablePaths) + 1 , : ) = tmpPath;
        nextNodes(i) = 0;
    elseif length( find( partialPath == nextNodes(i) ) ) ~= 0
        nextNodes(i) = 0;
    end
end
nextNodes = nextNodes(nextNodes ~= 0); %将nextNodes中为0的值去掉，因为下一个节点可能已经遍历过或者它就是目标节点
for i=1:length(nextNodes)
    tmpPath = cat(2, partialPath, nextNodes(i));
    tmpPsbPaths = findPath(Graph, tmpPath, destination, partialWeight + Graph(lastNode, nextNodes(i)));
    possiablePaths = cat(1, possiablePaths, tmpPsbPaths);
end
end
function [matrix_full] = spear2full_sq(matrix_sparse)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
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

%% 逐日仿真中单次计算函数(双路径)
function [ca,X,erfa,beita]=daytoday(erfa,m,du_s1,c2,t0,t1,thita,q,beita)
%经验性系数erfa=1，说明单独依靠前一次实际x
ca = cell(1,m,length(t0));%路网c1改为nodenum
X = cell(1,m,length(t0));%路网c1改为nodenum
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
        [X{i}]=logit_x(q,ca{i},thita,x_last,beita);% 单次分配 
end
end
%daytoday中function：time_1，preca，logit_x
%上一次的实际流量预测本次行程时间
 function [ca]=time_1(x,cap,t0)
 n=length(x);
 ca=zeros(1,n);
 for i=1:n
ca(1,i)=t0(1,i)*(1+0.15*(x(1,i)/cap(1,i))^4);
 end
 end
 %考虑习惯性选择的本次行程感知时间
 function  [ca]=preca(ca_last,ca_prlast,erfa)
 %经验性系数（习惯性系数）erfa=0，说明单独依靠前一次实际x，不考虑前一次预测的结果，如果之前处于平衡状态，这里也等于前一次实际的行程时间。
  n=length(ca_last);
 ca=zeros(n,1);
 for i=1:n
ca(1,i)=ca_last(1,i)*(1-erfa)+ca_prlast(1,i)*erfa;
 end
 end
 %双路径logit单次分配
 function [X]=logit_x(Q,ca,thita,y,beita)
ca_mean=(ca(1,1)+ca(1,2))/2;
ca_sum=(exp(-thita*ca(1,1)/ ca_mean))+exp(-thita*ca(1,2)/ ca_mean);
ptrans1=exp(-thita*ca(1,1)/ ca_mean)/ca_sum;
ptrans2=exp(-thita*ca(1,2)/ ca_mean)/ca_sum;
X=[0,0];
X(1,1)=Q*ptrans1*beita+y(1,1)*(1-beita);
X(1,2)=Q*ptrans2*beita+y(1,2)*(1-beita);
 end

% 画图
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
str1=['第',num2str(x),'-',num2str(figure_num),'(θ=',num2str(thita),')次仿真结果:行程时间逐次变化.','α=',num2str(erfa),',β=',num2str(beita)];
title(str1);
subplot(2,1,2)
plot(2,2,-39:m,x2(:,1),'r-*',-39:m,x2(:,2),'g-+');
str2=['第',num2str(x),'-',num2str(figure_num),'(θ=',num2str(thita),')次仿真结果:流量逐次变化.','α=',num2str(erfa),',β=',num2str(beita)];
title(str2);
end
%% 双路径平衡求解
function [s,t] = fsolve_flow(t01,t02,c1,c2,q,derta)
%UE平衡解析解
syms x1 x2
t1=t01*(1+0.15*(x1/c1)^4)-derta;
t2=t02*(1+0.15*(x2/c2)^4);
% t1=t01*(1+(x1/c1))-10;
% t2=t02*(1+(x2/c2));
eq1=t1-t2;
eq2=x1+x2-q;
s=solve(eq1,eq2,x1,x2);
digits(6); %设置运算位数
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
%SUE平衡解析解
syms x1 x2
t1=t01*(1+0.15*(x1/c1)^4);
t2=t02*(1+0.15*(x2/c2)^4);
% t1=t01*(1+(x1/c1));
% t2=t02*(1+(x2/c2));
eq1=sita*(t2-t1)-(log(x1)-log(x2));
eq2=x1+x2-q;
s=vpasolve(eq1,eq2,x1,x2);
% digits(6); %设置运算位数
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