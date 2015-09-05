%Here I make 3D sets of 2D plots to help visualize data

%sets defaults for the plots to look good
width = 5.1;     % Width in inches
height = 3;    % Height in inches
alw = 1;       % AxesLineWidth
fsz = 14;      % Fontsize
fna = 'Helvetica'; %Fontname
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize
interp = 'tex';


%Get Data from fig files.
fig10 = load('3D_DOS_10^7_Jx_100_h_0_2.fig','-mat');
fig14 = load('3D_DOS_10^7_Jx_143_h_0_2.fig','-mat');
fig03 = load('3D_DOS_10^7_Jx_30_h_0_2.fig','-mat');

[d1,d2] = fig10.hgS_070000.children.children.properties;
Ev10  = d1.XData;
dos10 = d1.YData;

[d1,d2] = fig14.hgS_070000.children.children.properties;
Ev14  = d1.XData;
dos14 = d1.YData;

[d1,d2] = fig03.hgS_070000.children.children.properties;
Ev03  = d1.XData;
dos03 = d1.YData;

%hgS_070000
%hgM_070000

%set up a vector of the same length with constant value
xv = Ev10*0+1;

%Plot the 2D plots in a 3D environment
hh = figure;

plot3(1.43*xv,Ev14,dos14,xv,Ev10,dos10,xv*.3,Ev03,dos03)

xs = .3;
x2= 1;
xb = 1.43;
yb = 17;
zb = 1.5;

axis([xs xb 0 yb 0 zb])
daspect([1.3 4 .8])
az=80;el=34;
view(az,el)


c = .2;
w =1.2;

line([xs,xs],[0,0],[0,zb],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([xs,xs],[0,yb],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

line([x2,x2],[0,0],[0,zb],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([x2,x2],[0,yb],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

% line([x3,x3],[0,0],[0,zb],'LineWidth',w,'Color',[c c c],'LineStyle','-')
% line([x3,x3],[0,yb],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

line([xb,xb],[0,0],[0,zb],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([xb,xb],[0,yb],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

grid on
set(gca,'XGrid','off');
set(gca,'YGrid','off');
%set(gca,'ZGrid','off');

set(gca,'XTick',[.3,1,1.43])
%set(gca,'YTick',[.1,.2])
%set(gca,'ZTick',[.1,.2])


title('DOS')
ylabel('\omega/J_z');
zlabel('DOS');
legend({'J_x = 1.43 J_z','J_x = 1.0 J_z','J_x=0.3 J_z'}, 'Location', 'NorthEast');
hold off;
filename = ['3D_DOSs_h_0'];
saveas(hh,filename)
print(hh, '-dpng', filename);
print(hh, '-depsc', filename);

