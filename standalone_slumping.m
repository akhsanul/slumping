th_depth=1;
mcr=0.2;

%%

x=1:1:100;
z=1:1:100;

z(1:41)=-6:0.05:-4;
z(41:53)=-4:0.5:2;
z(53:100)=2;
z=z+8;
%% define the permafrost line
nelayer=1:1:100;
nelayer=nelayer*0;
nelayer=nelayer+th_depth;
pline=z-nelayer;
z0=z;


%% find the gradient
mark=1:1:100;
mark=mark*0+1; % this is initial dummy variable to check the grid points for tigger conditions

tic % the function tic -tac is used to estimate the efficiency of the coding
f = figure;
f.Position = [100 100 800 200];
for i=1:1:1000000 % smaller number of iteration gives almost the same result
slope(1)=0; % first grid point can not have a slope, so zero is used as initial condition
for j=100:-1:2
    
    slope(j)=(z(j)-z(j-1))/(x(j)-x(j-1));
    if mark(j)==0
        slope(j)=mcr;
    end
end




% find all the indexes where the slope is greater than the critical slope
index=slope>mcr; % which index has critical slopes

% using find
index=1*index;
result=max(find(index==1)); % this will find the first index point where slumping starts
test=size(result);
if test(2)==0
    break
elseif test(2)==100
    print('simulation reached end of profile')
end

% now find the relevant points
% point B is the grid point where slumping starts
A.x=x(result+1);
B.x=x(result);
C.x=x(result-1);
D.x=x(result-2);



A.y=z(result+1);
B.y=z(result);
C.y=z(result-1);
D.y=z(result-2);


B2.y=B.y-(mcr-0.05)*(A.x-B.x);

% we need to check if the B2.y is lower
% than the permafrost line

if B2.y<pline(result)
    B2.y=pline(result);
    mark(result)=0;
else
    B2.y=B2.y;
end


% we need to find C2.y

%area before fail

% Area1=0.5*(A.y+B.y)*(A.x-B.x)+0.5*(C.y+B.y)*(B.x-C.x)+0.5*(C.y+D.y)*(C.x-D.x);
% 
% Area2=0.5*(A.y+B2.y)*(A.x-B.x);
% 
% syms y
% eqn=0.5*(B2.y+y)*(B.x-C.x)+0.5*(y+D.y)*(C.x-D.x)+Area2==Area1;
% 
% solx = solve(eqn, y);

A0=2*(0.5*(A.y+B.y)*(A.x-B.x)+0.5*(C.y+B.y)*(B.x-C.x)+0.5*(C.y+D.y)*(C.x-D.x));
B0=(A.y+B2.y)*(A.x-B.x);
rhs1=A0-B0;
rhs2=rhs1-(B.x-C.x)*B2.y-(C.x-D.x)*D.y;
y=rhs2/((B.x-C.x)+(C.x-D.x)); % y is the vertical position of point C

% now update the two points, B2 is the new position of B
z(result)=B2.y; % result is the index of point B 
z(result-1)=y;


changed_profile(i,:)=z;
% scatter(A.x,A.y)
% hold on,
% f = figure;
f.Position = [100 100 800 200];
plot(x,pline,'--k', 'linewidth',2)
hold on;
plot(x,z0,'r--', 'linewidth',1)
hold on
plot(x,z, 'b', 'linewidth',1)
ylim([1,12])
pause (.01)
hold off
end
toc
%%
% scatter(A.x,A.y)
% hold on,
choices=[1,100,400,900];
for k=choices;
f = figure;
f.Position = [100 100 800 200];
plot(x,pline,'-k', 'linewidth',1.2)
hold on;

p=area(x,pline);
p.FaceColor = 'blue';
p.FaceAlpha=.3;
p.EdgeColor = 'none';
% p.alpha=0.7;
hold on;
plot(x,z0,'r--', 'linewidth',0.7)
hold on
% plot(x,changed_profile(k,:), 'b', 'linewidth',0.2)
% hold on
area(x,changed_profile(k,:),'FaceColor','k','FaceAlpha',0.2)
ylim([1,12])
xlim([0, 80])
legend('permafrost table', 'permafrost','initial profile', 'thawed layer', 'Location','southeast', fontsize=14)
dim = [.15 .5 .35 .12];
str =append(sprintf('m_{cr}= %.2f,', mcr),sprintf(' initial thawing depth= %.2f m', th_depth)) ;
annotation('textbox',dim,'String',str,FontSize=14)
grid on
xlabel('(m)')
ylabel('(m)')
saveas(gcf,append(num2str(k),'.png'))
% pause (.01)
% hold off
end

