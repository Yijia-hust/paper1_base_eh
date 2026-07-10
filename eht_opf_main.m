clear;
clc;
yalmip('clear')
eht_opf_init;  % 系统拓扑等参数

%% ================== 定义未知量 (整合) ==================
% --- 电网侧变量 ---
PG=sdpvar(m,t,'full');
PL=sdpvar(k,t,'full');
theta=sdpvar(p,t,'full');
F_tu=sdpvar(m,t,'full');
ru=sdpvar(m,t,'full');
rd=sdpvar(m,t,'full');
DW=sdpvar(Nw,t,'full');
Dd=sdpvar(p,t,'full');

% --- 氢能侧变量  ---
% HGS 变量
PH=sdpvar(Nh,t,'full');
P_p2h=sdpvar(Nh,t,'full');
P_ep=sdpvar(Nh,t,'full');
Q_gen=sdpvar(Nh,t,'full');
Q_in=sdpvar(Nh,t,'full');
Q_out=sdpvar(Nh,t,'full');
Q=sdpvar(Nh,t,'full');
Q_hgsdt=sdpvar(Nh,t,'full');
Q_DT=sdpvar(1,t,'full');  % 假设只有1个HGS
SOC=sdpvar(Nh,t,'full');
% HLC 变量
Q_in_hlc  = sdpvar(Nc, t, 'full');
Q_hlc     = sdpvar(Nc, t, 'full');
HLC_load_shed = sdpvar(Nc,t,'full');
% 运输事件态变量
I = binvar(N_dt,N_HN,N_HN,'full');
t_arr = sdpvar(N_dt,N_HN,'full');
t_return = sdpvar(N_dt,1,'full');
Q_ser_hlc = sdpvar(N_dt,Nc,'full');
T_ser_hlc = sdpvar(N_dt,Nc,'full');
Q_dt_event = sdpvar(N_dt,1,'full');
Qv_arr = sdpvar(N_dt, N_HN, 'full');
Qv_dep = sdpvar(N_dt, N_HN, 'full');
U = sdpvar(N_dt, N_HN, 'full');
v=binvar(N_dt,Nc,'full');  % 与x_HLC仅有是否精确到时段t的区别
% 运输时段态与耦合变量
Q_dt = sdpvar(N_dt,t,'full');
x_HGS = binvar(N_dt, t, 'full');
x_HLC = binvar(N_dt, Nc, t, 'full');
q_ser_dit = sdpvar(N_dt, Nc, t, 'full');
% 松弛变量
s1=sdpvar(Nh,1,'full');
s2=sdpvar(Nh,1,'full');
% 增加出车成本
v_is_used = binvar(N_dt, 1, 'full'); % 车辆是否被使用的0-1变量

%% ================== 定义约束 (整合) ==================
Constraints = []; % 初始化约束集合

% ================== 第一部分: 电网与HGS经济调度 ==================
% 发电机实际输出功率限制
Constraints=[Constraints;repmat(PGmin,1,t)<=PG<=repmat(PGmax,1,t)];
% 节点功率平衡
Constraints=[Constraints;KG*PG+KW*(PW-DW)==KL*PL+KD*(bus_day-Dd)+KH*PH];
% 线路实际传输功率限制
Constraints=[Constraints;repmat(-PLmax,1,t)<=PL<=repmat(PLmax,1,t)];
% 直流潮流方程约束
Constraints=[Constraints;PL.*repmat(xl,1,t)==baseMVA*Ktheta*theta];
% 定义参考节点
Constraints=[Constraints;theta(1,:)==0];
% 节点电压相角约束
Constraints=[Constraints;-pi<=theta<=pi];
% 发电机成本线性化
Constraints=[Constraints;F_tu>=0];
for j=1:K
    Constraints=[Constraints;F_tu>=repmat(A(:,j),1,t).*PG+repmat(B(:,j),1,t)];
end
% 爬坡约束
if t~=1
    Constraints=[Constraints;repmat(-RD,1,t-1)*delta_t<=PG(:,2:t)-PG(:,1:t-1)<=repmat(RU,1,t-1)*delta_t];
end
% 上旋转备用约束
Constraints=[Constraints;0<=ru<=repmat(RU,1,t)];
Constraints=[Constraints;ru<=repmat(PGmax,1,t)-PG];
Constraints=[Constraints;sum(ru,1)>=repmat(SRU,1,t)];
% 下旋转备用约束
Constraints=[Constraints;0<=rd<=repmat(RD,1,t)];
Constraints=[Constraints;rd<=PG-repmat(PGmin,1,t)];
Constraints=[Constraints;sum(rd,1)>=repmat(SRD,1,t)];
% 弃风量限制
Constraints=[Constraints;0<=DW<=PW];
% 切负荷量约束
Constraints=[Constraints;0<=Dd<=bus_day];

% --- HGS 内部逻辑 (已验证) ---
% HGS耗电与产氢关系
Constraints=[Constraints;P_p2h*eta_p2h*delta_t==Q_gen];
Constraints=[Constraints;P_ep*delta_t==Q_gen*eta_h2c];
Constraints=[Constraints;PH==P_p2h+P_ep];
Constraints=[Constraints;PHmin<=PH<=PHmax];
% HGS氢气平衡
Constraints=[Constraints;Q_gen==Q_in+Q_hgsdt];
Constraints=[Constraints;Q_DT==Q_out+Q_hgsdt];
Constraints=[Constraints;Q_DT(1,:)==sum(Q_dt,1)];
Constraints=[Constraints;0<=Q_hgsdt];
Constraints=[Constraints;Q(:,1)==Q_init+Q_in(:,1)-Q_out(:,1)];
for i=2:t
    Constraints=[Constraints;Q(:,i)==Q(:,i-1)+Q_in(:,i)-Q_out(:,i)];
end
Constraints=[Constraints;SOC==Q./QN];
Constraints=[Constraints;Qmin<=Q<=Qmax];
% --- HGS储氢罐调度周期末状态约束 ---
Constraints=[Constraints; SOC(:,t) == SOC_init+s1-s2];
Constraints=[Constraints;s1>=0];
Constraints=[Constraints;s2>=0];
% Constraints=[Constraints; SOC(:,t) == SOC_init];
Constraints=[Constraints;0<=Q_in<=Qin_max];
Constraints=[Constraints;0<=Q_out<=Qout_max];

% 20260711增加
Constraints = [Constraints; 0 <= Q_gen <= Qgen_max];  % 限制HGS每小时最大制氢量
Constraints = [Constraints; 0 <= Q_DT  <= QDT_hour_max];  % 限制车辆在HGS处每小时最大装载量


% ================== 第二、三、四部分: 运输与HLC逻辑 (核心) ==================
% --- HLC 内部逻辑 ---time-based
for hh=1:Nc   % 总量守恒（假设每辆车只对每个hlc加一次
    for i=1:t
        Constraints = [Constraints; Q_in_hlc(hh,i) == sum(q_ser_dit(:,hh,i),1)];
    end
end
Constraints = [Constraints; Q_hlc(:,1) == Q_hlc_init + Q_in_hlc(:,1) - (HLC_load(:,1)-HLC_load_shed(:,1))];
for i=2:t
    Constraints = [Constraints; Q_hlc(:,i) == Q_hlc(:,i-1) + Q_in_hlc(:,i) - (HLC_load(:,i)-HLC_load_shed(:,i))];
end
Constraints = [Constraints; 0 <= HLC_load_shed <= HLC_load];
Constraints = [Constraints; Q_hlc >= 0];
Constraints = [Constraints; Q_hlc <= repmat(Q_hlc_cap, 1, t)];
for i=1:t
  Constraints = [Constraints; 0 <= Q_in_hlc(:,i) <= Qin_hlc_max]; % HLC站每小时接收氢气的上限
end

% --- 运输与耦合逻辑 ---
for d=1:N_dt
    for n=1:N_HN
        Constraints = [Constraints; I(d,n,n) == 0];  % 禁止自环（文中没有）
    end
end
for d = 1:N_dt
    for m = 1:N_HN
        Constraints = [Constraints; sum(I(d,m,:)) <= 1];  % 6a，节点的出弧限制在1条。这样可允许车辆访问多个HLC
    end
end
for d = 1:N_dt
    for n = 1:N_HN
        in_flow = sum(I(d, :, n));
        out_flow = sum(I(d, n, :));
        Constraints = [Constraints; in_flow == out_flow];  % 6b
    end
end
for d=1:N_dt
%     Constraints = [Constraints; sum(I(d,id_hgs,:)) <= 1];  % 允许车辆不出任务
%     Constraints = [Constraints; sum(I(d,:,id_hgs)) <= 1];  % 因为inflow=outflow，所以只要出发就必返回
    Constraints = [Constraints; sum(I(d,id_hgs,:)) == v_is_used(d)]; % 离开HGS的次数等于“被使用”状态
    Constraints = [Constraints; sum(I(d,:,id_hgs)) == v_is_used(d)]; % 返回HGS的次数也等于“被使用”状态
end
for d=1:N_dt
    for j=2:N_HN
        for k=2:N_HN
            if j~=k
                Constraints = [Constraints; U(d,j) - U(d,k) + N_HN*I(d,j,k) <= N_HN-1];
            end
        end
    end
    Constraints = [Constraints; 0 <= U(d,2:end) <= N_HN-1];
end
% Constraints = [Constraints; T_ser_hlc == 0];
Constraints = [Constraints; T_ser_hlc == 0.001*Q_ser_hlc];
for d=1:N_dt  % 如果放松模型可以注释掉大M约束的>=边
    for j=1:N_HN  % 出发点 j
        for k=1:N_HN  % 到达点 k
            if j~=k
                if Adj(j,k)==0
                    % 不相连的边，一定不存在路径
                    Constraints = [Constraints; I(d,j,k) == 0];
                    continue; % 跳过下面所有关于这条边的约束
                end

                % 对于相连的边，保留原来逻辑
                % 情况 1: 从 HGS (j=id_hgs) 出发
                if j == id_hgs
                    expression = t_arr(d,j) + T_tra(j,k) - t_arr(d,k);
                    Constraints = [Constraints; expression <= M_time * (1 - I(d,j,k))];
                    Constraints = [Constraints; expression >= -M_time * (1 - I(d,j,k))];
                
                % 情况 2: 从 HLC (j~=id_hgs) 出发
                else
                    hlc_idx = j - 1; % HLC 的索引
                    % 情况 2a: 从 HLC 返回 HGS (k=id_hgs)
                    if k == id_hgs
                        expression = t_arr(d,j) + T_ser_hlc(d,hlc_idx) + T_tra(j,k) - t_return(d);
                        Constraints = [Constraints; expression <= M_time * (1 - I(d,j,k))];
                        Constraints = [Constraints; expression >= -M_time * (1 - I(d,j,k))];
                    % 情况 2b: 从 HLC 前往另一个 HLC
                    else
                        expression = t_arr(d,j) + T_ser_hlc(d,hlc_idx) + T_tra(j,k) - t_arr(d,k);
                        Constraints = [Constraints; expression <= M_time * (1 - I(d,j,k))];
                        Constraints = [Constraints; expression >= -M_time * (1 - I(d,j,k))];
                    end
                end
            end
        end
    end
end

Constraints = [Constraints; 0 <= t_arr <= t];
Constraints = [Constraints; 0 <= t_return <= t];
for d=1:N_dt
    Constraints = [Constraints; Qv_dep(d,id_hgs) == Qv_arr(d,id_hgs) + Q_dt_event(d) ];
end
for d=1:N_dt
    for m=1:Nc  % 这里m和n是同一个点
        n = id_hlc(m);
        Constraints = [Constraints; Qv_dep(d,n) == Qv_arr(d,n) - Q_ser_hlc(d,m)];
    end
end
for d=1:N_dt
    Constraints = [Constraints; 0 <= Qv_arr(d,:) <= Q_dt_cap; 0 <= Qv_dep(d,:) <= Q_dt_cap];  % 文献里加了0-1约束，只有车子经过这个点时生效
    Constraints = [Constraints; Qv_arr(d,id_hgs) == 0]; %目前你每辆车只跑一趟（一次出发一次返回，返回HGS时氢量必须为0），恰好合理。
    % 但如果将来扩展到多次往返，车返回后还没卸完氢，这个约束就会强制它把氢留在路上，导致infeasible或者错误解
end
for d=1:N_dt
    for j=1:N_HN
        for k=1:N_HN
            if j~=k
                if Adj(j,k)==0  % 不相连边，不存在路径
                    Constraints = [Constraints; I(d,j,k) == 0];
                else
                    expression = Qv_dep(d,j) - Qv_arr(d,k);  % 路上无损耗
                    Constraints = [Constraints; expression <= M_h2*(1-I(d,j,k)); expression >= -M_h2*(1-I(d,j,k))];
                end
            end
        end
    end
end
epsilon = 1e-4;
for d=1:N_dt
%     Constraints = [Constraints; sum(x_HGS(d,:)) <= 1]; % 每辆车最多装载1次
    for i=1:t
        Constraints = [Constraints; t_arr(d,id_hgs) >= i - M_time * (1 - x_HGS(d,i))];
        Constraints = [Constraints; t_arr(d,id_hgs) <= (i + 1 - epsilon) + M_time * (1 - x_HGS(d,i))];
    end
    Constraints = [Constraints; sum(Q_dt(d,:)) == Q_dt_event(d)];  % 增加约束
end
for d=1:N_dt
    for i=1:t
        Constraints = [Constraints; Q_dt_min * x_HGS(d,i) <= Q_dt(d,i) <= Q_dt_cap * x_HGS(d,i)];
    end
end
for d=1:N_dt
    for m=1:Nc
        n = id_hlc(m);
%         Constraints = [Constraints; sum(x_HLC(d,m,:)) <= 1];
        for i=1:t
            Constraints = [Constraints; t_arr(d,n) >= (i - 1 + epsilon) - M_time * (1 - x_HLC(d,m,i))];
            Constraints = [Constraints; t_arr(d,n) <= i + M_time * (1 - x_HLC(d,m,i))];
            Constraints = [Constraints; Q_ser_min * x_HLC(d,m,i) <= q_ser_dit(d,m,i) <= Q_ser_max * x_HLC(d,m,i)];
        end
        Constraints = [Constraints; sum(q_ser_dit(d,m,:),3) == Q_ser_hlc(d,m)];  % 增加约束
    end
end


% 改回HGS处帽子约束
for d=1:N_dt
    for i=1:t
        Constraints = [Constraints; Q_dt(d,i) <= Q_dt_event(d,1) + M_h2*(1 - x_HGS(d,i))];
        Constraints = [Constraints; Q_dt(d,i) >= Q_dt_event(d,1) - M_h2*(1 - x_HGS(d,i))];
    end
end
% 改回HLC处帽子约束
for d=1:N_dt
    for m = 1:Nc
        for i=1:t
            Constraints = [Constraints; q_ser_dit(d,m,i) <= Q_ser_hlc(d,m) + M_h2*(1 - x_HLC(d,m,i))];
            Constraints = [Constraints; q_ser_dit(d,m,i) >= Q_ser_hlc(d,m) - M_h2*(1 - x_HLC(d,m,i))];
        end
    end
end
% 添加sum(x_HGS)<=sum(I)--约束11c
for d = 1:N_dt
    % Constraints = [Constraints; sum(x_HGS(d,:)) <= sum(I(d,id_hgs,:))];
    Constraints = [Constraints;sum(x_HGS(d,:)) == v_is_used(d);];  % 每辆出车车辆必须且只能有一个装氢时段，未出车车辆不能产生装氢事件
end
% 添加z=sum(x_HLC)<=sum(I)--约束12c
for d = 1:N_dt
    for m=1:Nc
        n = id_hlc(m); % 获取HLC对应的节点编号
        Constraints = [Constraints; v(d,m) == sum(x_HLC(d,m,:),3)];
        Constraints = [Constraints; sum(x_HLC(d,m,:),3) <= sum(I(d,n,:))];
    end
end

% 补充时间窗约束
for d=1:N_dt
  Constraints = [Constraints; tw_min(d,id_hgs) <= t_arr(d,id_hgs) <= tw_max(d,id_hgs)];
  for m=1:Nc
    n = id_hlc(m);
    Constraints = [Constraints; tw_min(d,n) <= t_arr(d,n) <= tw_max(d,n)];
  end
end

% 一天结束后，HLC库存保持在一定的水平，能满足第二天前几个小时的氢负荷。
Constraints = [Constraints; Q_hlc(:,t) >= Q_hlc_init];

%% ================== 求解 (整合版) ==================
tic
% 最终的经济目标函数
% p 是母线数量, 根据 bus 矩阵的行数自动获取
high_penalty_cost = 1000;  % 设置一个较高的惩罚成本, e.g., 1000 $/MWh
c_Dd = zeros(1, p); % 初始化一个全零的 1xp 向量

% 找到所有Pd>0的负荷母线，并为它们赋予惩罚成本
% load_bus_indices = find(bus(:, 3) > 0);
% c_Dd(load_bus_indices) = high_penalty_cost;
% 将更重要的负荷节点的成本设置得更高
c_Dd(1)  = 1000; % 节点 1 的切负荷成本
c_Dd(2)  = 1000;
c_Dd(3)  = 1200; % 例如，可以把节点3的负荷看得更重要一些
c_Dd(4)  = 900;
c_Dd(5)  = 900;
c_Dd(6)  = 1000;
c_Dd(7)  = 1100;
c_Dd(8)  = 1000;
c_Dd(9)  = 1100;
c_Dd(10) = 1200;
c_Dd(13) = 1500; % 例如，节点13的负荷非常重要
c_Dd(14) = 1000;
c_Dd(15) = 1200;
c_Dd(16) = 900;
c_Dd(18) = 1500;
c_Dd(19) = 1100;
c_Dd(20) = 1000;
c_Dd(21) = 1000;
c_Dd(22) = 1000;
c_Dd(23) = 1200;


b_tra = 5; % 每单位路途的成本系数，根据[41]Table I，成本为 5 $/KM
travel_cost = 0;
for d=1:N_dt
    for m=1:N_HN
        for n=1:N_HN
            if Adj(m,n) == 1   % 只对相连边计成本
                travel_cost = travel_cost + I(d,m,n) * b_tra * x(m,n);
            end
        end
    end
end
penalty_h2_shed = 1000;
vv=500; % 每出1辆车是500元-----不要这个+vv*sum(v_is_used)了
objective = sum(F_tu(:)) + [50,48]*sum(DW,2) + c_Dd*sum(Dd,2)+travel_cost+penalty_h2_shed*sum(HLC_load_shed(:))+5000*(sum(s1)+sum(s2));

% options=sdpsettings('solver','mosek', 'verbose', 2);
% options = sdpsettings('solver', 'mosek', 'verbose', 2, 'mosek.MSK_DPAR_MIO_TOL_REL_GAP', 0.02);  % 写0.02约20秒解出，写0.013约2分钟解出，但结果完全不同且都有错
options = sdpsettings('solver', 'gurobi', 'verbose', 2, 'gurobi.MIPGap', 0.001); 
% options = sdpsettings('solver', 'gurobi', 'verbose', 2);  % 不允许gap


% % ======================= NaN自动诊断脚本 =======================
% tic
% fprintf('\n\n--- 开始逐个检查约束是否存在NaN问题 ---\n');
% 
% F_all = Constraints; % 备份所有约束
% problem_found = false;
% 
% for i = 1:length(F_all)
%     fprintf('正在测试第 %d / %d 个约束...\n', i, length(F_all));
%     F_test = F_all(i); % 只取出当前这一个约束
%     
%     try
%         % 我们用一个虚拟的目标函数来测试这个约束能否被编译
%         % 'verbose',0 让它保持安静，我们只关心是否报错
%         % 'solver','' 让YALMIP只做编译，不调用求解器
%         options_test = sdpsettings('verbose',0,'solver','');
%         
%         % 关键一步：尝试编译这个单独的约束
%         sol = optimize(F_test, [], options_test);
%         
%         % 检查YALMIP的内部错误标志
%         if sol.problem ~= 0 && sol.problem ~= 1 && sol.problem ~= 4
%              % 如果错误不是“成功”、“无解”或“找不到求解器”，那就是编译错误
%              error('YALMIP编译时返回错误代码: %d. 信息: %s', sol.problem, sol.info);
%         end
%         
%         fprintf('--- 第 %d 个约束... 通过!\n', i);
%         
%     catch ME
%         fprintf('\n!!! >>> 错误源头找到! 问题出在第 %d 个约束! <<< !!!\n', i);
%         fprintf('这个有问题的约束是: \n');
%         disp(F_test); % 打印出有问题的约束
%         fprintf('YALMIP在编译它时产生的错误信息是: %s\n', ME.message);
%         problem_found = true;
%         break; % 找到错误就停止
%     end
% end
% 
% if ~problem_found
%     fprintf('\n--- 所有约束单独检查均通过。问题可能由约束之间的相互作用引起。---\n');
% end
% 
% fprintf('--- 诊断结束 ---\n');
% toc
% % ======================= 脚本结束 =======================

results=optimize(Constraints,objective,options);
if results.problem==0
    disp('求解成功')
else
    disp('求解失败')
    disp(results.info)
end

toc

results

% --- 结果展示 ---
cost=double(objective);
PG_show=double(PG);
PL_show=double(PL);
PH_show=double(PH);
DW_show = double(DW);
Dd_show = double(Dd);
Q_gen_show=double(Q_gen);
SOC_show = double(SOC);
Q_hlc_show = double(Q_hlc);
HLC_load_shed_show = double(HLC_load_shed);
I_show = double(I);  % I的第一个维度是车辆编号，可以按I_show(1,:,: )查看
t_arr_show = double(t_arr);
t_return_show = double(t_return);
Q_dt_event_show = double(Q_dt_event);
Q_ser_hlc_show = double(Q_ser_hlc);
v_is_used_show = double(v_is_used);

Q_DT_show = double(Q_DT);
T_ser_hlc_show = double(T_ser_hlc);
x_HGS_show = double(x_HGS);
x_HLC_show = double(x_HLC);
T_tra_show = double(T_tra);
% 检查每个HLC被派的车辆数量，不为0则说明有车到达并服务
hlc_visited = squeeze(any(I_show(:, :, 2:end), [1 2]));
disp(hlc_visited)

% 对每辆车，把未访问节点的t_arr置为0
for d = 1:N_dt
    for n = 1:N_HN
        % 检查节点n是否被车d访问
        if sum(I_show(d,:,n)) < 0.5 && n ~= id_hgs
            t_arr_show(d,n) = 0;
        end
    end
end
disp('到达时间 t_arr（未访问节点显示0）:')
disp(t_arr_show)


% 检查结果1：在命令行运行以下代码，可以看每辆车走了哪条路线
for d = 1:N_dt
    if v_is_used_show(d) < 0.5
        fprintf('车%d: 未出车\n', d); continue;
    end
    current = id_hgs;
    fprintf('车%d: HGS', d);
    while true
        for next = 1:N_HN
            if I_show(d, current, next) > 0.5
                if next == id_hgs
                    fprintf('->HGS\n');
                else
                    fprintf('->HLC%d', next-1);
                end
                current = next;
                break;
            end
        end
        if current == id_hgs, break; end
    end
end

% 检查结果2：看时间和氢量
disp('到达时间 t_arr:')
disp(t_arr_show)

disp('装载量 Q_dt_event:')
disp(Q_dt_event_show)

disp('卸货量 Q_ser_hlc:')
disp(Q_ser_hlc_show)

disp('是否出车:')
disp(v_is_used_show)

% 检查结果3：看HLC库存
disp('HLC库存 Q_hlc:')
disp(Q_hlc_show)

disp('切氢负荷:')
disp(HLC_load_shed_show)

% 查看电网侧结果
% 弃风量
fprintf('各时段弃风量(MW):\n')
disp(DW_show)
fprintf('总弃风量(MWh): %.2f\n', sum(DW_show(:)))

% 切负荷
fprintf('总切负荷(MWh): %.2f\n', sum(Dd_show(:)))

% 发电量
fprintf('各机组总发电量(MWh):\n')
disp(sum(PG_show, 2))

% HGS耗电
fprintf('HGS各时段耗电(MW):\n')
disp(PH_show)

% SOC
fprintf('储氢罐SOC:\n')
disp(SOC_show)


fprintf('总缺氢量：%.4f kg\n',sum(HLC_load_shed_show(:)));
fprintf('期末SOC偏差 s1：%.6f\n',value(s1));
fprintf('期末SOC偏差 s2：%.6f\n',value(s2));
fprintf('车辆总装氢量：%.4f kg\n',sum(Q_dt_event_show));

% 想检查是否存在“同一时段同时充氢和放氢”
Q_in_show = double(Q_in);
Q_out_show = double(Q_out);
tol = 1e-6;
sim_idx = find(Q_in_show > tol & Q_out_show > tol);

if isempty(sim_idx)
    disp('未出现同时充氢和放氢');
else
    fprintf('同时充放氢的时段：');
    fprintf('%d ',sim_idx);
    fprintf('\n最大同时充放量：%.4f kg\n', ...
        max(min(Q_in_show(sim_idx),Q_out_show(sim_idx))));
end

% 约束残差，若结果不低于约 -1e-5 就可用
residual = check(Constraints);
fprintf('最小约束裕度：%.3e\n',min(residual));

% 判断HLC期末库存是否回到初始库存
tol_hlc = 1e-4;
hlc_terminal_diff = Q_hlc_show(:,end) - Q_hlc_init;

if max(abs(hlc_terminal_diff)) <= tol_hlc
    fprintf('HLC期末库存已回到初始库存。\n');
else
    fprintf('HLC期末库存未回到初始库存；总增加量：%.4f kg，最大单站偏差：%.4f kg。\n', ...
        sum(hlc_terminal_diff), max(abs(hlc_terminal_diff)));
end

% 方便比较不同gap
fprintf('目标函数值：%.4f\n',cost);