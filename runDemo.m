function runDemo()
    %% — 参数 —
    p.g = 9.81;
    p.R = 0.028575;                  % 球半径 (m)
    p.m = 0.17;                       % 质量 (kg)
    p.I = 2/5*p.m*p.R^2;              % 实心球转动惯量

    p.Lx = 2.54; p.Ly = 1.27;         % 9英尺桌面近似 (m)

    p.mu_k = 0.2; p.mu_r = 0.01; p.c_spin=0.02;% 台布滑动摩擦系数，台布滚动摩擦系数，自旋衰减（空气阻力）
    p.e_ball = 0.93; p.mu_ball = 0.2; %球恢复系数，球滑动摩擦系数
    p.e_cushion = 0.85; p.mu_cushion = 0.2; %库边恢复系数，库边动摩擦系数
    p.e_tip = 0.9; p.mu_tip = 0.5; %杆恢复系数，杆动摩擦系数

    p.dt = 0.0001;p.dt_store=0.02; p.T=5.0; p.v_cut=0.0;

    %% — 初始状态 —
    balls = newBallArray(2, p);
    balls(1).pos = [-0.3,0.0,p.R]; balls(1).vel=[0.0 0.0 0.0]; balls(1).omega=[0.0 0.0 0.0];
    balls(1).Rmat=eye(3); balls(1).struck=false;
    balls(2).pos=[0.5,0.01,p.R]; balls(2).Rmat=eye(3);
    sim_stop=false;

    %% — 击球设定 —
    cue.time=0.5;
    cue.point=balls(1).pos;
    cue.dir_tip=normalize([1.00 0.00 -1.0]);
    cue.v_tip=0.5;
    cue.offset=[0.005 0.0];%半球面到整个yz的仿射变换。[offset_y,offset_z]

    %% — 仿真数据存储初始化 —
    nFrames = ceil(p.T/p.dt_store);
    nBalls = numel(balls);
    positions = zeros(nFrames, nBalls, 3);
    Rmats = zeros(3,3,nBalls,nFrames); % 3x3旋转矩阵，球编号，帧编号
    for i = 1:nBalls
        Rmats(:,:,i,1) = eye(3); % 初始旋转矩阵
    end

    f_forces = zeros(nFrames, nBalls, 3); % 用于存摩擦力

    %% — 主循环（存储位置，不绘图)
    t = 0; frame_idx = 1; t_next_store = 0;
    
    while t < p.T && ~sim_stop
        % 击球事件
        if ~balls(1).struck && t >= cue.time
            balls(1)=applyCueStrike(balls(1),cue,p);
            balls(1).struck=true;
        end
    
        % 推进一步
        balls=stepOnce(balls,p);
    
        t = t + p.dt;
    
        % 存储指定帧
        if t >= t_next_store - p.dt/2
            for i=1:nBalls
                positions(frame_idx,i,:) = balls(i).pos;
                Rmats(:,:,i,frame_idx) = balls(i).Rmat;
                f_forces(frame_idx,i,:) = balls(i).f_friction;
            end
            frame_idx = frame_idx + 1;
            t_next_store = t_next_store + p.dt_store;
        end
    
        % 停止判据
        if t >= cue.time && all(arrayfun(@(b) norm(b.vel(1:2))<p.v_cut && norm(b.omega)<0.2, balls))
            sim_stop = true; 
        end
    end
    disp('Simulation finished.')


    %% — 拖动回放阶段 —
    replaySimulation(positions,Rmats,f_forces, nBalls, p,cue);
end

function replaySimulation(positions, Rmats, f_forces, nBalls, p, cue)
    nFrames = size(positions,1);
    f = figure('Name','Billiards3D Replay','Color','w');
    ax = axes('Parent',f); hold(ax,'on'); axis(ax,'equal'); grid(ax,'on');
    xlabel('x'); ylabel('y'); zlabel('z'); view(ax,30,20);
    xlim(ax,[-p.Lx/2 p.Lx/2]); ylim(ax,[-p.Ly/2 p.Ly/2]); zlim(ax,[0 0.15]);
    drawTable(ax,p);

    cueLength = 1.5; % 杆长度 m
    hCue = plot3(ax, [0 0], [0 0], [0 0], 'Color','b','LineWidth',2);

    hBalls = gobjects(1,nBalls); 
    hMarks = gobjects(1,nBalls);
    hTraj  = gobjects(1,nBalls);
    hFriction = gobjects(1,nBalls); % 摩擦力箭头

    for i = 1:nBalls
        color = 'k'; if i==1,color='w'; elseif i==2,color='r'; end
        [hBalls(i), hMarks(i)] = drawBall(ax, squeeze(positions(1,i,:))', p.R, color);
        hTraj(i) = plot3(ax, positions(1,i,1), positions(1,i,2), positions(1,i,3), ...
                         '-', 'Color', get(hBalls(i),'FaceColor'), 'LineWidth',1.2, 'Visible','off');
        % 初始化摩擦力箭头
        pos = squeeze(positions(1,i,:));
        hFriction(i) = quiver3(pos(1), pos(2), pos(3), 0, 0, 0, 'r', 'LineWidth',1.5,'MaxHeadSize',0.5);
    end

    timeText = text(ax, -p.Lx/2, p.Ly/2, 0.15, sprintf('t = %.2f s',0), ...
                    'FontSize', 12, 'Color', 'k', 'FontWeight','bold');

    slider = uicontrol('Style','slider','Min',1,'Max',nFrames,'Value',1, ...
        'SliderStep',[1/(nFrames-1), 0.1],'Position',[50 20 400 20], ...
        'Callback', @(src,~) redraw(round(src.Value)));

    isPlaying = false;
    btn = uicontrol('Style','togglebutton','String','Play','Position',[460 20 60 20], ...
        'Callback', @(src,~) togglePlay(src));

    chk = uicontrol('Style','checkbox','String','Show Trajectory','Position',[530 20 120 20], ...
        'Value',0,'Callback',@(src,~) toggleTrajectory(src));

    % ✅ 新增：摩擦力箭头开关
    chkFric = uicontrol('Style','checkbox','String','Show Friction','Position',[660 20 120 20], ...
        'Value',1,'Callback',@(src,~) toggleFriction(src));

    tmr = timer('ExecutionMode','fixedRate','Period',0.02, 'TimerFcn', @(~,~) nextFrame());

    currentFrame = 1;
    redraw(currentFrame);

    function redraw(frameIdx)
        currentFrame = frameIdx;
        set(slider,'Value',currentFrame);

        for i = 1:nBalls
            pos = squeeze(positions(frameIdx,i,:));
            % 更新球体
            set(hBalls(i),'XData',p.R*hBalls(i).UserData.Xs+pos(1), ...
                          'YData',p.R*hBalls(i).UserData.Ys+pos(2), ...
                          'ZData',p.R*hBalls(i).UserData.Zs+pos(3));

            % 更新标记点
            R = Rmats(:,:,i,frameIdx);  
            markR = hMarks(i).UserData.markR;
            local = [p.R; 0; 0];  
            world = pos(:) + R*local;
            mx = hMarks(i).UserData.mx; my = hMarks(i).UserData.my; mz = hMarks(i).UserData.mz;
            set(hMarks(i), 'XData', markR*mx + world(1), ...
                           'YData', markR*my + world(2), ...
                           'ZData', markR*mz + world(3));
            % 更新轨迹
            if strcmp(get(hTraj(i),'Visible'),'on')
                set(hTraj(i),'XData',positions(1:frameIdx,i,1), ...
                             'YData',positions(1:frameIdx,i,2), ...
                             'ZData',positions(1:frameIdx,i,3));
            end
            % 更新摩擦力箭头（缩放5倍显示）
            Fvec = squeeze(f_forces(frameIdx,i,:)) / 0.5;
            set(hFriction(i),'XData', pos(1),'YData',pos(2),'ZData',pos(3)-p.R,...
                             'UData', Fvec(1),'VData', Fvec(2),'WData', Fvec(3));
        end

        % 更新杆子（只显示母球未被击打前一小段时间）
        b1pos = squeeze(positions(1,1,:))';

        if frameIdx* p.dt_store < cue.time  % 前0.2s显示杆
            % 杆方向：从杆头到击球点
            cueEnd = b1pos + offset2contact(cue.offset,p.R);
            cueStart = cueEnd - cue.dir_tip*cueLength; % 杆尾
            set(hCue, 'XData',[cueStart(1), cueEnd(1)], ...
                      'YData',[cueStart(2), cueEnd(2)], ...
                      'ZData',[cueStart(3), cueEnd(3)]);
            set(hCue,'Visible','on');
        else
            set(hCue,'Visible','off');
        end

        simTime = (frameIdx-1) * p.dt_store;
        set(timeText, 'String', sprintf('t = %.2f s', simTime));
        drawnow limitrate nocallbacks;
    end

    function togglePlay(src)
        if src.Value
            src.String = 'Pause'; isPlaying = true; start(tmr);
        else
            src.String = 'Play';  isPlaying = false; stop(tmr);
        end
    end

    function toggleTrajectory(src)
        if src.Value, set(hTraj,'Visible','on'); else set(hTraj,'Visible','off'); end
    end

    function toggleFriction(src)
        if src.Value
            set(hFriction,'Visible','on');
        else
            set(hFriction,'Visible','off');
        end
    end

    function nextFrame()
        if currentFrame < nFrames
            redraw(currentFrame + 1);
        else
            stop(tmr); isPlaying = false; btn.Value = 0; btn.String = 'Play';
        end
    end
end

%% ====== 工具函数 ======

function balls = newBallArray(N, p)
    balls = repmat(struct(...
        'pos',[0 0 p.R], ...
        'vel',[0 0 0], ...
        'omega',[0 0 0], ...
        'color',[0.6 0.6 0.6], ...
        'hBall',[], ...
        'hMark',[], ...
        'Rmat',eye(3), ...
        'R',p.R, ...
        'struck',false, ...
        'f_friction',[0 0 0] ... % 新增字段
    ),1,N);
end

function [hBall,hMark] = drawBall(ax,pos,R,color)
    [Xs,Ys,Zs]=sphere(20);
    hBall=surf(ax,R*Xs+pos(1),R*Ys+pos(2),R*Zs+pos(3),'FaceColor',color,'EdgeColor','none','FaceLighting','gouraud');

    % 标记点（球前方小黑点）
    [mx,my,mz]=sphere(6);
    markPos=pos+[R,0,0];
    markR = 0.08*R;  % 半径
    hMark = surf(ax, markR*mx+markPos(1), markR*my+markPos(2), markR*mz+markPos(3), ...
             'FaceColor', 'k', ...          % 黑色填充
             'EdgeColor', 'k', ...          % 白色边框
             'LineWidth', 1.2, ...          % 边框更粗
             'FaceLighting', 'gouraud');    % 光照更柔和;

    % 保存网格用于更新
    hBall.UserData.Xs=Xs; hBall.UserData.Ys=Ys; hBall.UserData.Zs=Zs;
    hMark.UserData.mx=mx; hMark.UserData.my=my; hMark.UserData.mz=mz;
    hMark.UserData.markR=markR;
end

function Rnew = updateBallRotation(R, omega, dt)
    % omega: 角速度向量 [wx wy wz]
    theta = norm(omega)*dt;
    if theta < 1e-12
        Rnew = R;
        return;
    end
    axis = omega / norm(omega);
    Rdelta = axisAngle2Rotm(axis, theta);
    Rnew = Rdelta * R; % 左乘表示增量旋转
end

function R = axisAngle2Rotm(axis,theta)
    axis=axis/norm(axis);
    c=cos(theta); s=sin(theta); t=1-c;
    x=axis(1); y=axis(2); z=axis(3);
    R=[t*x^2+c t*x*y-s*z t*x*z+s*y;
       t*x*y+s*z t*y^2+c t*y*z-s*x;
       t*x*z-s*y t*y*z+s*x t*z^2+c];
end

%% ===== 其他原有物理函数保持不变 =====
% drawTable, applyCueStrike, stepOnce, collideCushion, collideInfiniteMass, collideBalls, normalize
% 你原来的逻辑保持，只要调用 updateBallGraphics 更新网格即可

function drawTable(ax, p)
    % 台面（绿色矩形）
    [X,Y] = meshgrid(linspace(-p.Lx/2,p.Lx/2,2), linspace(-p.Ly/2,p.Ly/2,2));
    Z = 0*X;
    s = surf(ax, X, Y, Z, 'FaceColor',[0.05 0.4 0.2], 'FaceAlpha',0.9, 'EdgeColor','none'); %#ok<NASGU>
    % 库边
    plot3(ax, [-p.Lx/2 p.Lx/2 p.Lx/2 -p.Lx/2 -p.Lx/2], ...
              [-p.Ly/2 -p.Ly/2 p.Ly/2 p.Ly/2 -p.Ly/2], ...
              [0 0 0 0 0], 'w-', 'LineWidth',1.2);
end

function b = applyCueStrike(b, cue, p)
    % 击球：在母球上施加瞬时冲量（法向 + 切向），更新 (v, omega)
    % contact：相对球心的击打点位置（单位 m），限定在球面内
    % 直接用杆头击球方向的反向作为法向量
    contact = offset2contact(cue.offset,p.R);

    n = normalize(contact);  

    v_tip_vec = cue.v_tip * normalize(cue.dir_tip);
    % 接触点速度差：球面点 vs 杆头
    g = (b.vel + cross(b.omega, p.R*n)) - v_tip_vec;
    gn = dot(g,n);

    % 法向冲量（杆头等效无限质量）
    Jn = -(1+p.e_tip)*p.m*gn;

    % 切向：理想粘着所需冲量
    gt = g - gn*n; gt_norm = norm(gt);
    if gt_norm>1e-9
        t_hat = gt/gt_norm;
        Jt_star = -(2*p.m/7)*dot(gt,t_hat);
        Jt = Jt_star;
        if abs(Jt) > p.mu_tip*abs(Jn)
            Jt = -p.mu_tip*abs(Jn)*sign(dot(gt,t_hat));
        end
        J = Jn*n + Jt*t_hat;
    else
        J = Jn*n;
        Jt = 0; %#ok<NASGU>
    end
    
    J_linear = J;
    J_linear(3) = 0;% 归并约束力冲量

    % 更新
    b.vel = b.vel + J_linear/p.m;
    b.omega = b.omega + (cross(p.R*n, J))/p.I;
end

function balls = stepOnce(balls, p)
    dt = p.dt;
    N = numel(balls);

    for i=1:N
        b = balls(i);
        n = [0 0 1]; r = -p.R*n;  
        v_contact = b.vel + cross(b.omega, r);
        gt = v_contact; % 切向接触速度
        gt_norm = norm(gt);
        
        % 初始化摩擦力
        f_friction = [0 0 0];
        
        if gt_norm > 1e-6
            % 滑动摩擦
            Ft = -p.mu_k*p.m*p.g * gt;
            a = Ft/p.m;
            alpha = cross(r, Ft)/p.I;
            b.vel = b.vel + a*dt;
            b.omega = b.omega + alpha*dt;
            f_friction = Ft; % 保存摩擦力
        else
            % 近似纯滚
            vxy = b.vel; vxy(3)=0; vxy_n = norm(vxy);
            if vxy_n>1e-8
                a = -p.mu_r*p.m*p.g * vxy / (p.m*vxy_n);
                b.vel = b.vel + a*dt;
                f_friction = -p.mu_r*p.m*p.g * (vxy/vxy_n);
            end
            b.omega = b.omega - p.c_spin*[0 0 b.omega(3)]*dt;
        end
        
        % 重力约束
        b.pos = b.pos + b.vel*dt;
        b.pos(3) = p.R; b.vel(3)=0; 
        b.Rmat = updateBallRotation(b.Rmat, b.omega, p.dt);

        % 保存摩擦力到结构体
        b.f_friction = f_friction;
        balls(i) = b;
    end

    % 碰撞处理：库边
    for i=1:N
        balls(i) = collideCushion(balls(i), p);
    end
    % 两球对碰
    for i=1:N-1
        for j=i+1:N
            [balls(i), balls(j)] = collideBalls(balls(i), balls(j), p);
        end
    end
end


function b = collideCushion(b, p)
    % 与四条边的碰撞：按瞬时冲量处理
    n = [0 0 1]; r = -p.R*n;  %#ok<NASGU>
    % X 方向边
    if b.pos(1) >  p.Lx/2 - p.R
        nx = [-1 0 0]; b = collideInfiniteMass(b, nx, p);
        b.pos(1) =  p.Lx/2 - p.R;
    elseif b.pos(1) < -p.Lx/2 + p.R
        nx = [1 0 0];  b = collideInfiniteMass(b, nx, p);
        b.pos(1) = -p.Lx/2 + p.R;
    end
    % Y 方向边
    if b.pos(2) >  p.Ly/2 - p.R
        ny = [0 -1 0]; b = collideInfiniteMass(b, ny, p);
        b.pos(2) =  p.Ly/2 - p.R;
    elseif b.pos(2) < -p.Ly/2 + p.R
        ny = [0 1 0];  b = collideInfiniteMass(b, ny, p);
        b.pos(2) = -p.Ly/2 + p.R;
    end
end

function b = collideInfiniteMass(b, n, p)
    % 与静止刚性平面的碰撞（库边），法向 n 指向台内
    r = -p.R*n;
    v_contact = b.vel + cross(b.omega, r);
    gn = dot(v_contact, n); gt = v_contact - gn*n; gt_n = norm(gt);

    % 法向冲量（平面等效无限大质量）
    Jn = -(1+p.e_cushion)*p.m*gn;

    if gt_n > 1e-9
        t_hat = gt/gt_n;
        Jt_star = -(2*p.m/7)*dot(gt,t_hat);
        Jt = Jt_star;
        if abs(Jt) > p.mu_cushion*abs(Jn)
            Jt = -p.mu_cushion*abs(Jn)*sign(dot(gt,t_hat));
        end
        J = Jn*n + Jt*t_hat;
    else
        J = Jn*n;
    end

    b.vel = b.vel + J/p.m;
    b.omega = b.omega + cross(r, J)/p.I;
end

function [b1, b2] = collideBalls(b1, b2, p)
    % 两球相交则按冲量求解
    d = b2.pos - b1.pos; dxy = d; dxy(3)=0; dist = norm(dxy);
    if dist >= 2*p.R - 1e-6
        return; % 未接触
    end
    if dist < 1e-6
        % 极端重合：微调
        dxy = [1,0,0]; dist = 1; 
    else
        dxy = dxy/dist;
    end
    n = dxy;                     % 法向：1->2
    r1 = p.R*n; r2 = -p.R*n;

    % 接触点相对速度
    v1c = b1.vel + cross(b1.omega, r1);
    v2c = b2.vel + cross(b2.omega, r2);
    g = v1c - v2c;
    gn = dot(g,n); gt = g - gn*n; gt_n = norm(gt);

    % 法向冲量（等质量）
    Jn = -(1+p.e_ball)*(gn) / (1/p.m + 1/p.m); % = -(1+e) m/2 * gn

    if gt_n>1e-9
        t_hat = gt/gt_n;
        Jt_star = -(dot(gt,t_hat)) / ( (1/p.m+1/p.m) + (p.R^2/p.I + p.R^2/p.I) );
        % 对实心球化简 = - (m/7) * dot(gt,t_hat)
        Jt = Jt_star;
        if abs(Jt) > p.mu_ball*abs(Jn)
            Jt = -p.mu_ball*abs(Jn)*sign(dot(gt,t_hat));
        end
        J = Jn*n + Jt*t_hat;
    else
        J = Jn*n;
    end

    % 施加冲量
    b1.vel = b1.vel +  J/p.m;  b2.vel = b2.vel - J/p.m;
    b1.omega = b1.omega + cross(r1, J)/p.I;
    b2.omega = b2.omega + cross(r2, -J)/p.I;

    % 穿透修正：把两球推回切触
    overlap = 2*p.R - dist;
    if overlap>0
        corr = 0.5*overlap*n;
        b1.pos = b1.pos - corr; b2.pos = b2.pos + corr; 
    end
end

% ========================= 小工具 =========================
function v = normalize(v)
    n = norm(v);
    if n<1e-12, return; end
    v = v/n;
end

function contact = offset2contact(offset,R)

    if numel(offset) == 2
        y = offset(1); z = offset(2);
    elseif numel(offset) == 3
        y = offset(2); z = offset(3);
    else
        error('offset 必须是2或3维向量');
    end
    
    x = -R;
    alpha = R / sqrt(x^2+y^2+z^2);
    
    % 输出向量
    contact = alpha * [x, y, z];
end
