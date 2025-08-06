function [] = Func10_DesignproblemMesh(fem,~,~)

%% Dimension
X0 = 45;                
X1 = 80;                
X2 = 140;                
X3 = 80;               
X4 = 45;                

Y0 = 65;                
Y1 = 5;                
Y2 = 75;                 
Y3 = 20;                
Y4 = 55;                
Y5 = 20;                
Y6 = 10;                

% Generate the meshes and nodes number
dx = 5;                 % X direction mesh step
dy = 5;                 % Y direction mesh step

% Node number information
x0=round(X0/dx+1); 
x1=round(x0+X1/dx); 
x2=round(x1+X2/dx);
x3=round(x2+X3/dx); 
x4=round(x3+X4/dx);

NodeX=x4;               % Nodes in X direction

y0=round(Y0/dy+1);
y1=round(y0+Y1/dy);
y2=round(y1+Y2/dy); 
y3=round(y2+Y3/dy);
y4=round(y3+Y4/dy);
y5=round(y4+Y5/dy);
y6=round(y5+Y6/dy);

NodeY=y6;               % Nodes in Y direction
%%
% Find x coordinates values for material critical nodes

for j=1:NodeX
    if(j<=x0)        Xaxis(j)=(j-1)*dx;end
    if(j>x0)&&(j<=x1) Xaxis(j)=(j-x0)*dx+X0;end
    if(j>x1)&&(j<=x2) Xaxis(j)=(j-x1)*dx+X0+X1; end
    if(j>x2)&&(j<=x3) Xaxis(j)=(j-x2)*dx+X0+X1+X2; end
    if(j>x3)&&(j<=x4) Xaxis(j)=(j-x3)*dx+X0+X1+X2+X3; end
end

for i=1:NodeY
    if(i<=y0)        Yaxis(i)=(i-1)*dy; end
    if(i>y0)&&(i<=y1) Yaxis(i)=(i-y0)*dy+Y0; end
    if(i>y1)&&(i<=y2) Yaxis(i)=(i-y1)*dy+Y0+Y1; end
    if(i>y2)&&(i<=y3) Yaxis(i)=(i-y2)*dy+Y0+Y1+Y2; end
    if(i>y3)&&(i<=y4) Yaxis(i)=(i-y3)*dy+Y0+Y1+Y2+Y3; end
    if(i>y4)&&(i<=y5) Yaxis(i)=(i-y4)*dy+Y0+Y1+Y2+Y3+Y4; end
    if(i>y5)&&(i<=y6) Yaxis(i)=(i-y5)*dy+Y0+Y1+Y2+Y3+Y4+Y5; end
end
%%
% Plot model
    PP(1,:)=[0, 0];
    PP(2,:)=[X0+X1+X2+X3+X4,0];
    PP(3,:)=[0, Y0];
    PP(4,:)=[X0+X1+X2+X3+X4, Y0];
    PP(5,:)=[X0, Y0+Y1];
    PP(6,:)=[X0+X1, Y0+Y1];
    PP(7,:)=[X0+X1+X2, Y0+Y1];
    PP(8,:)=[X0+X1+X2+X3, Y0+Y1];
    PP(9,:)=[X0+X1, Y0+Y1+Y2];
    PP(10,:)=[X0+X1+X2, Y0+Y1+Y2];
    PP(11,:)=[X0+X1, Y0+Y1+Y2+Y3];
    PP(12,:)=[X0+X1+X2, Y0+Y1+Y2+Y3];
    PP(13,:)=[X0, Y0+Y1+Y2+Y3+Y4];
    PP(14,:)=[X0+X1, Y0+Y1+Y2+Y3+Y4]; 
    PP(15,:)=[X0+X1+X2, Y0+Y1+Y2+Y3+Y4];
    PP(16,:)=[X0+X1+X2+X3, Y0+Y1+Y2+Y3+Y4];
    PP(17,:)=[X0+X1, Y0+Y1+Y2+Y3+Y4+Y5];
    PP(18,:)=[X0+X1+X2, Y0+Y1+Y2+Y3+Y4+Y5];    
    PP(19,:)=[0, Y0+Y1+Y2+Y3+Y4+Y5+Y6];
    PP(20,:)=[X0+X1+X2+X3+X4, Y0+Y1+Y2+Y3+Y4+Y5+Y6];

% Plot iron
    IronX1=[PP(1,1) PP(2,1) PP(4,1) PP(3,1)];
    IronY1=[PP(1,2) PP(2,2) PP(4,2) PP(3,2)];
    IronX2=[PP(5,1) PP(6,1) PP(11,1) PP(12,1),PP(7,1),PP(8,1), PP(16,1), PP(15,1), PP(14,1), PP(13,1)];
    IronY2=[PP(5,2) PP(6,2) PP(11,2),PP(12,2),PP(7,2),PP(8,2), PP(16,2), PP(15,2), PP(14,2), PP(13,2)];

% Plot Coil1
    Coil1X1=[PP(9,1) PP(10,1) PP(12,1) PP(11,1)];
    Coil1Y1=[PP(9,2) PP(10,2) PP(12,2) PP(11,2)];

% Plot Coil2
    Coil2X1=[PP(14,1) PP(15,1) PP(18,1) PP(17,1)];
    Coil2Y1=[PP(14,2) PP(15,2) PP(18,2) PP(17,2)];

% Plot Meshes
    figure(9);

    fill(IronX1,IronY1,"y"); hold on;
    fill(IronX2,IronY2,"y"); hold on;
    fill(Coil1X1,Coil1Y1,"r"); hold on;
    fill(Coil2X1,Coil2Y1,"r"); hold on;
    axis equal; axis([0 390/2 0 250]);

    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Design Setting Mesh');
    grid on;


    for e=1:fem.ne
        i = fem.IX(e,1);        j = fem.IX(e,2);        k = fem.IX(e,3);

        xi = fem.X(i,1);  yi = fem.X(i,2);
        xj = fem.X(j,1);  yj = fem.X(j,2);
        xk = fem.X(k,1);  yk = fem.X(k,2);

    line([xi,xj],[yi,yj]); line([xj,xk],[yj,yk]); line([xk,xi], [yk,yi]);
    end

end
