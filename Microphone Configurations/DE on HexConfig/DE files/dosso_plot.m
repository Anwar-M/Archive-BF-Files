% program for visualisation of "dosso's function

clear all
close all

% calculate dosso function
x = -3:0.05:3;
for m1 = 1:length(x)
      z1(m1) = dosso([x(m1) 0     0     0     0     0]); 
      z2(m1) = dosso([0     x(m1) 0     0     0     0]); 
      z3(m1) = dosso([0     0     x(m1) 0     0     0]); 
      z4(m1) = dosso([0     0     0     x(m1) 0     0]); 
      z5(m1) = dosso([0     0     0     0     x(m1) 0]); 
      z6(m1) = dosso([0     0     0     0     0     x(m1)]); 
end

% plotting
figure(1)
subplot(231)
plot(x,z1,'k')
xaxis([-3 3])
h = get(gca,'Children');
set(h,'LineWidth',1)
h = get(gcf,'Children');
set(h,'FontSize',16)
set(h,'LineWidth',1)
xlabel('m_1')
subplot(232)
plot(x,z2,'k')
xaxis([-3 3])
h = get(gca,'Children');
set(h,'LineWidth',1)
h = get(gcf,'Children');
set(h,'FontSize',16)
set(h,'LineWidth',1)
xlabel('m_2')
subplot(233)
plot(x,z3,'k')
xaxis([-3 3])
h = get(gca,'Children');
set(h,'LineWidth',1)
h = get(gcf,'Children');
set(h,'FontSize',16)
set(h,'LineWidth',1)
xlabel('m_3')
subplot(234)
plot(x,z4,'k')
xaxis([-3 3])
h = get(gca,'Children');
set(h,'LineWidth',1)
h = get(gcf,'Children');
set(h,'FontSize',16)
set(h,'LineWidth',1)
xlabel('m_4')
ylabel('E_1')
subplot(235)
plot(x,z5,'k')
xaxis([-3 3])
h = get(gca,'Children');
set(h,'LineWidth',1)
h = get(gcf,'Children');
set(h,'FontSize',16)
set(h,'LineWidth',1)
xlabel('m_5')
subplot(236)
plot(x,z6,'k')
xaxis([-3 3])
h = get(gca,'Children');
set(h,'LineWidth',1)
h = get(gcf,'Children');
set(h,'FontSize',16)
set(h,'LineWidth',1)
xlabel('m_6')

clear z1 z2 z3 z4 z5 z6

% calculate dosso function
x = -3:0.01:3;
y = -3:0.01:3;
for m1 = 1:length(x)
   for m2 = 1:length(y)   
      z1(m1,m2) = dosso([x(m1) y(m2) 0     0     0     0]); 
      z2(m1,m2) = dosso([0     0     x(m1) y(m2) 0     0]); 
      z3(m1,m2) = dosso([0     0     0     0     x(m1) y(m2)]); 
  end
end

N = 14;

figure(2)
pcolor(x,y,z1')
shading interp
colormap(gray(N))
caxis([0 1])
colorbar
hold on
v = [0:1/(N):1];
contour(x,y,z1',v,'k')
h = get(gca,'Children'); 
set(h,'LineWidth',1)
h = get(gcf,'Children'); 
set(h(1),'Fontsize',[16])
set(h(2),'Fontsize',[16])
xlabel('m_1')
ylabel('m_2')
hold off

figure(3)
pcolor(x,y,z2')
shading interp
colormap(gray(N))
caxis([0 1])
colorbar
hold on
contour(x,y,z2',v,'k')
h = get(gca,'Children'); 
set(h,'LineWidth',1)
h = get(gcf,'Children'); 
set(h(1),'Fontsize',[16])
set(h(2),'Fontsize',[16])
xlabel('m_3')
ylabel('m_4')
hold off

figure(4)
pcolor(x,y,z3')
shading interp
colormap(gray(N))
caxis([0 1])
colorbar
hold on
contour(x,y,z3',v,'k')
h = get(gca,'Children'); 
set(h,'LineWidth',1)
h = get(gcf,'Children'); 
set(h(1),'Fontsize',[16])
set(h(2),'Fontsize',[16])
xlabel('m_5')
ylabel('m_6')
hold off

% end of program