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
e=[0 0.1 0.8 2];
xr=-a/2 + a*(linspace(1,nb,nb));  %Centres of Periodic Potential Barriar


%% Generation of H-matrix
H=zeros(4,N,N);

E0=((pi^2)/L^2)*linspace(1,N,N).^2; %Energies of N-States of Inf. Potential Well


m=1;
while m<N+1
    j=1;
    while j<N/2+1
       t1=V0*h1(m,j);
       l=1;
       t2=h2(m,j);
       while l<5
            
            H(l,m,j)=(m==j)*E0(m)+t1+t2*e(l);
            H(l,j,m)=(m==j)*E0(m)+t1+t2*e(l);
            l=l+1;          
       end
       j=j+1;
    end
    
    m=m+1;
    
end



%% Eigenvalues and Eigenvectors of the energy levels
V=zeros(4,N,N);
D=zeros(4,N,N);
T=zeros(4,N);
m=zeros(4,1);
c=['r' 'b' 'k' 'm'];
l=1;
while l<5
    [V(l,:,:),D(l,:,:)]=eig(reshape(H(l,:,:),[N,N]));
    
    T(l,:)=eig(reshape(H(l,:,:),[N,N]));
    m(l)=find(T(l,:)>0.5,1);    
    l=l+1;
end

%% PLOTS
l=1;

figure(1)
scatter((1/L)*linspace(1,40,40),E0(1:40),'DisplayName','Inf. Pot')

while l<5
    
    %Energy levels vs k (Energy Band Gap is easily visible) 
    figure(1)
    scatter((pi/L)*linspace(1,40,40),T(l,m(l):m(l)+39),'o',c(l),'filled','DisplayName',strcat('\epsilon = ',num2str(e(l))))
    hold on
    legend ('Location','northwest')
    xlabel("k")
    ylabel("E_{n}")

    %Formation of wavefunction via Matrix Solution of Schrodinger's Equation
    n=1;     %Energy Level to be plotted
    si=0;    
    k=1;
    
    while k<N+1
        si=si+V(l,k,n-1+m(l))*((2/L)^0.5)*sin(k*pi*x/L);  %Addition of Harmonics
        k=k+1;
    end



    %Wavefunction Plot
    figure (2)
    fplot(si,[0,L],c(l),'DisplayName',strcat('\epsilon = ',num2str(e(l))))
    hold on
    axis padded
    grid on
    title("\psi_{"+num2str(n)+"}(x)")
    xlabel("x") 
    legend
    
    l=l+1;
end

%% Intermediary Functions
function y=h1(n,m)
    global xr b nb
    i=1;
    y=0;
    while i<nb+1
        y=y+(F(n,m,xr(i)+b/2)-F(n,m,xr(i)-b/2));
        i=i+1;
    end
    

end

function y=h2(m,n)
    global L
    if m==n
      y=L/2;
    elseif mod(m+n,2)==0
        y=0;
    else
        y=-8*m*n*L/((pi^2)*(m^2-n^2)^2);         
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




