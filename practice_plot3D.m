%set(fh, ‘color’, ‘white’);

vals = (1:100)/100;
nums = vals*0 +1;
plot3(nums*0,vals,sin(vals).*(1-vals).*vals,nums*.2 ,vals,5*(vals-.5).^2.*(1-vals).*vals,nums*.4 ,...
    vals,5*(vals-.5).^2.*(1-vals).*vals,nums*.6 ,vals,(.5-(vals-.6).^2).*(1-vals).*vals)
%plot3(nums*0,vals,sin(vals),nums*.2 ,vals,(vals-.5).^2)
axis([0 .6 0 1 0 .2])
daspect([1.8 1 .5])
az=-110;el=12;
view(az,el)

%line([0,0],[1,1],[0,.2],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')
%line([0,0],[0,1],[0,0],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')


% line([.2,.2],[1,1],[0,.2],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')
% line([.2,.2],[0,1],[0,0],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')
% 
% line([.4,.4],[1,1],[0,.2],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')
% line([.4,.4],[0,1],[0,0],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')


%line([.6,.6],[1,1],[0,.2],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')
%line([.6,.6],[0,1],[0,0],'LineWidth',1.5,'Color',[.5 .5 .5],'LineStyle','--')
c = .2;
w =1.2;

line([0,0],[0,0],[0,.2],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([0,0],[0,1],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

line([.2,.2],[0,0],[0,.2],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([.2,.2],[0,1],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

line([.4,.4],[0,0],[0,.2],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([.4,.4],[0,1],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

line([.6,.6],[0,0],[0,.2],'LineWidth',w,'Color',[c c c],'LineStyle','-')
line([.6,.6],[0,1],[0,0],'LineWidth',w,'Color',[c c c],'LineStyle','-')

grid on
set(gca,'XGrid','off');
set(gca,'YGrid','off');
%set(gca,'ZGrid','off');

set(gca,'ZTick',[.1,.2])

