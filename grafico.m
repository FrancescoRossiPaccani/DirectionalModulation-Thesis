%% Antenna and Users Parameters
ulaSize = 3;          % Number of antennas
lambda = 0.19;        % Fictitious wavelength
spacing = lambda/2;

% Array position (aligned along y-axis)
antennaLocation = zeros(3,ulaSize);
antennaLocation(2,:) = (0:ulaSize-1) * spacing;

% Intended receiver
theta_ref = deg2rad(0);
phi_ref   = deg2rad(-20);

% Eavesdropper 1
theta_int = deg2rad(0);
phi_int   = deg2rad(25);

% Eavesdropper 2
theta_int2 = deg2rad(0);
phi_int2   = deg2rad(-50);

%% Convert to Cartesian coordinates
r = 10; % Fictitious distance from the transmitter
[x_user, y_user, z_user]     = sph2cart(phi_ref, theta_ref, r);
[x_int1, y_int1, z_int1]     = sph2cart(phi_int, theta_int, r);
[x_int2, y_int2, z_int2]     = sph2cart(phi_int2, theta_int2, r);

%% 3D Plot
figure('Color','w'); hold on; grid on; box on;

% Antenna array
plot3(antennaLocation(1,:), antennaLocation(2,:), antennaLocation(3,:), ...
      'ks','MarkerSize',12,'MarkerFaceColor','k');

% Users
plot3(x_user, y_user, z_user, 'go','MarkerSize',12,'MarkerFaceColor','g');
plot3(x_int1, y_int1, z_int1, 'ro','MarkerSize',12,'MarkerFaceColor','r');
plot3(x_int2, y_int2, z_int2, 'mo','MarkerSize',12,'MarkerFaceColor',[1 0.5 0]);

% Lines from array to users
line([0 x_user],[0 y_user],[0 z_user],'Color','g','LineStyle','--','LineWidth',2,'HandleVisibility','off');
line([0 x_int1],[0 y_int1],[0 z_int1],'Color','r','LineStyle','--','LineWidth',2,'HandleVisibility','off');
line([0 x_int2],[0 y_int2],[0 z_int2],'Color',[1 0.5 0],'LineStyle','--','LineWidth',2,'HandleVisibility','off');

% Labels using LaTeX
%text(x_user+.5, y_user + deg2rad(5), z_user+0.5, '$Intended\ Receiver$', 'Interpreter','latex', 'FontSize',14, 'Color','g','FontWeight','bold');
%text(x_int1+0.5, y_int1, z_int1+0.5, '$Eavesdropper\ 1$', 'Interpreter','latex', 'FontSize',14, 'Color','r','FontWeight','bold');
%text(x_int2+0.5, y_int2, z_int2+0.5, '$Eavesdropper\ 2$', 'Interpreter','latex', 'FontSize',14, 'Color',[1 0.5 0],'FontWeight','bold');

% Optional: show beam direction as a cone
coneLength = r; 
coneRadius = 0.5; 
[Xc,Yc,Zc] = cylinder([0 coneRadius]);
Zc = Zc*coneLength;
h = surf(Xc,Yc,Zc,'FaceAlpha',0.2,'EdgeColor','none','FaceColor','g');
v = [x_user, y_user, z_user]; % direction vector
theta = acos(v(3)/norm(v));
phi   = atan2(v(2),v(1));
rotate(h,[0 1 0], rad2deg(theta), [0 0 0]);
rotate(h,[0 0 1], rad2deg(phi), [0 0 0]);

% Axis labels with LaTeX
xlabel('$X$ [m]','Interpreter','latex','FontSize',14); 
ylabel('$Y$ [m]','Interpreter','latex','FontSize',14); 
zlabel('$Z$ [m]','Interpreter','latex','FontSize',14); 

% Title with LaTeX
%title('Antenna Array and Users Positions in 3D Space','Interpreter','latex','FontSize',16);

% Legend with LaTeX
legend({'Antenna Array','$\text{Intended Receiver}$','$\text{Eavesdropper 1}$','$\text{Eavesdropper 2}$'},...
    'Location','northeast','Interpreter','latex', 'FontSize', 12);

% 3D view
view(15,15);      % Azimuth and elevation
axis equal;
grid on;
rotate3d on;      % Enable interactive rotation
