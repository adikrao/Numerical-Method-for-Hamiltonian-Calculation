%%KP MODEL WITHOUT EXTERNAL FIELD


%% Global Variables determining System constants
syms x

global V0 L nb b xr a
V0=100;   %Magnitude of Periodic Potential Barriar 
N=150;    % No. of eigenvalues used for formulation
a=1;      % Lattice constant
nb=20;    %No. of barriars
L=nb*a;   %Length of Inf. Potential Well
b=1/16;   %Width of Potential Barrier
E_up=80;  %Limit of energy for plotting
xr=-a/2 + a*(linspace(1,nb,nb));  %Centres of Periodic Potential Barriar


%% Generation of H-matrix
H=zeros(N,N);

E0=((pi^2)/L^2)*linspace(1,N,N).^2; %Energies of N-States of Inf. Potential Well


k=1;
while k<N+1
    j=1;
    while j<N/2+1
       H(k,j)=V0*h(k,j);
       H(j,k)=V0*h(k,j);
       j=j+1;
    end
    k=k+1;
    
end

H=diag(E0)+H;

%% Eigenvalues and Eigenvectors of the energy levels

[V,D]=eig(H);
T1=eig(H);

index=find(T1>0,1);

% Eliminating negative potentials since E>0
T=T1(index:end);
Vnew=V(:,index:end);  %c_n vector for each energy level

%% PLOTS

%Energy levels vs k (Energy Band Gap is easily visible) 
figure (1)
num=floor(sqrt(E_up*L^2/pi^2));
scatter((pi/L)*linspace(1,num,num),T(1:num),'DisplayName','KP Model')
hold on
scatter((pi/L)*linspace(1,num,num),E0(1:num),'DisplayName','Inf. Well')
xlabel("k")
ylabel("E_{n}")

%Analytical KP Model
% P=V0*(b/2)*(a-b);
% alp=linspace(2.113,pi,100);
% alp2=linspace(4.546,2*pi,100);
% 
% E=alp.^2;
% E2=alp2.^2;
% k=acos(P*(sin(alp)./alp) +cos(alp));
% k2=2*pi-acos(P*(sin(alp2)./alp2) +cos(alp2));
% plot(k,E)
% plot(k2,E2)


%KP model
E=linspace(0,E_up,1000);
alp=sqrt(E);
beta=sqrt(V0-E);
fun=((beta.^2-alp.^2).*sin(alp*(a-b)).*sinh(beta*b)./(2*alp.*beta)) + cos(alp*(a-b)).*cosh(beta*b);
counter=0;
ind=cell(4,1);
count=cell(4,1);
n=1;
for i=1:1000
    if (abs(fun(i))>1)
        if counter>0
        ind{n}=i;  
        count{n}=counter;
        n=n+1;
        end
        E(i)=0;
        fun(i)=0;
        counter=0;
    else
        counter=counter+1;
        
    end
    
    
end
if counter>0
ind{n}=i;  
count{n}=counter;
end        

k=(1/a)*acos(fun);
plot(k(ind{2}-count{2}:ind{2}-1),E(ind{2}-count{2}:ind{2}-1),'r','DisplayName','Analytical KP Band-1')
plot(2*pi-k(ind{3}-count{3}:ind{3}-1),E(ind{3}-count{3}:ind{3}-1),'r','DisplayName','Analytical KP Band-2')
plot(2*pi+k(ind{4}-count{4}+1:ind{4}-1),E(ind{4}-count{4}+1:ind{4}-1),'r','DisplayName','Analytical KP Band-3')
title("E-k Diagram")
legend('Location','northwest')
hold off



i=1;
n=[1 4 7 13];
figure (2)

while i<5
    %Formation of wavefunction via Matrix Solution of Schrodinger's Equation
     %Energy Level to be plotted
si=0;   
q=1;
while q<N+1
    si=si-Vnew(q,n(i))*((2/L)^0.5)*sin(q*pi*x/L);  %Addition of Harmonics
    q=q+1;
end

%Wavefunction Plot
subplot(2,2,i)
fplot(si,[0,L])
hold on
fplot(0.4*sin(n(i)*pi*x/L),[0,L],'--')
axis padded
title("\psi_{"+num2str(n(i))+"}(x)")
xlabel("x")
i=i+1;
grid on
hold off
end


%% Intermediary Functions
function x=h(n,m)
    global xr b nb
    i=1;
    x=0;
    while i<nb+1
        x=x+(F(n,m,xr(i)+b/2)-F(n,m,xr(i)-b/2));
        i=i+1;
    end
    

end


function f=F(n,m,x)
global L
if n==m
    f=(x/L) - (sin(2*pi*n*x/L)/2*pi*n);
else
    f=(sin((m-n)*pi*x/L)/(pi*(m-n)))  -  (sin((m+n)*pi*x/L)/(pi*(m+n)));
       
end

end




