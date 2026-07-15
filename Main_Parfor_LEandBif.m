clc; clear; close all;

%% ================= SYSTEM PARAMETERS =================
SIGMA = 16;
RHO   = 45.92;

q1=0.99; q2=0.99; q3=0.99;

h = 1e-3;

N = 20;            % Number of integration steps between GSR operations
h_norm = N*h;      % Physical duration of one GSR interval

tn = 100;
T = 0:h:tn;
n = length(T);

%% ================= FRACTIONAL COEFFICIENTS =================
cp1=1; cp2=1; cp3=1;
c1=zeros(1,n); c2=zeros(1,n); c3=zeros(1,n);
for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1; cp1=c1(j);
    c2(j)=(1-(1+q2)/j)*cp2; cp2=c2(j);
    c3(j)=(1-(1+q3)/j)*cp3; cp3=c3(j);
end

%% ================= PARAMETER SCAN =================
b_min = 1;
b_max = 10;
b_step = 0.1;
BETA_vals = b_min:b_step:b_max;
K = length(BETA_vals);

L1 = nan(1,K); L2 = nan(1,K); L3 = nan(1,K);
bif_cell = cell(K,1);

%% ================= PARALLEL =================
if isempty(gcp('nocreate'))
    parpool;
end

parfor kpar = 1:K
iter_tic = tic;
    BETA = BETA_vals(kpar);

    %% Initial conditions
    x=zeros(1,n); y=zeros(1,n); z=zeros(1,n);
    x(1)=20; y(1)=10; z(1)=50;

    %% Jacobian init
    f11=zeros(1,n); f12=zeros(1,n); f13=zeros(1,n);
    f21=zeros(1,n); f22=zeros(1,n); f23=zeros(1,n);
    f31=zeros(1,n); f32=zeros(1,n); f33=zeros(1,n);

    % Phi(0) = I_3
    f11(1)=1;
    f22(1)=1;
    f33(1)=1;

    CUM = zeros(3,1);
    gsr_count = 0;

    peak_buffer=[];

    num_steps = n-1;

    % Transient duration aligned with a complete GSR interval
    tr_steps = round(0.30*num_steps/N)*N;


    for k=2:n

        %% ----- State update -----
        x(k)=(SIGMA*(y(k-1)-x(k-1)))*h^q1 - memo(x,c1,k);
        y(k)=(RHO*x(k)-y(k-1)-x(k)*z(k-1))*h^q2 - memo(y,c2,k);
        z(k)=(x(k)*y(k)-BETA*z(k-1))*h^q3 - memo(z,c3,k);

        %% ----- Jacobian update -----
        f11(k)=(SIGMA*(f21(k-1)-f11(k-1)))*h^q1 - memo(f11,c1,k);
        f12(k)=(SIGMA*(f22(k-1)-f12(k-1)))*h^q1 - memo(f12,c1,k);
        f13(k)=(SIGMA*(f23(k-1)-f13(k-1)))*h^q1 - memo(f13,c1,k);

        f21(k)=(RHO*f11(k-1)-f21(k-1)-f11(k-1)*z(k-1)-x(k-1)*f31(k-1))*h^q2 - memo(f21,c2,k);
        f22(k)=(RHO*f12(k-1)-f22(k-1)-f12(k-1)*z(k-1)-x(k-1)*f32(k-1))*h^q2 - memo(f22,c2,k);
        f23(k)=(RHO*f13(k-1)-f23(k-1)-f13(k-1)*z(k-1)-x(k-1)*f33(k-1))*h^q2 - memo(f23,c2,k);

        f31(k)=(f11(k-1)*y(k-1)+x(k-1)*f21(k-1)-BETA*f31(k-1))*h^q3 - memo(f31,c3,k);
        f32(k)=(f12(k-1)*y(k-1)+x(k-1)*f22(k-1)-BETA*f32(k-1))*h^q3 - memo(f32,c3,k);
        f33(k)=(f13(k-1)*y(k-1)+x(k-1)*f23(k-1)-BETA*f33(k-1))*h^q3 - memo(f33,c3,k);

        %% ----- GSR -----
        % MATLAB index k corresponds to mathematical time step k-1
        if mod(k-1,N) == 0

        J = [f11(k) f12(k) f13(k);
         f21(k) f22(k) f23(k);
         f31(k) f32(k) f33(k)];

        [J,E] = GSR(J);

        % Store reorthonormalized tangent vectors
        f11(k)=J(1,1); f12(k)=J(1,2); f13(k)=J(1,3);
        f21(k)=J(2,1); f22(k)=J(2,2); f23(k)=J(2,3);
        f31(k)=J(3,1); f32(k)=J(3,2); f33(k)=J(3,3);

        % Accumulate growth factors only after the transient
        if (k-1) > tr_steps
        CUM = CUM + log(E(:));
        gsr_count = gsr_count + 1;
        end
        end

        %% ----- Bifurcation peak detection -----
    if (k-1) > tr_steps && k > 2
    if x(k-1)>x(k-2) && x(k-1)>x(k)
        peak_buffer(end+1)=x(k-1);
    end
    end
    end

    if gsr_count > 0
    LE_last = CUM/(gsr_count*h_norm);

    L1(kpar) = LE_last(1);
    L2(kpar) = LE_last(2);
    L3(kpar) = LE_last(3);
    end

    bif_cell{kpar}=peak_buffer;
fprintf('Iteration = %d | BETA = %.4f | Elapsed time = %.2f s\n', ...
        kpar, BETA, toc(iter_tic));
end

%% ================= BIFURCATION MATRIX =================
D = zeros(0,2);
for k=1:K
    pks = bif_cell{k};
    if isempty(pks), continue; end
    D = [D; BETA_vals(k)*ones(numel(pks),1), pks(:)];
end

%% ================= PLOTS =================
figure;
plot(BETA_vals,L1,'k','LineWidth',1.2); hold on;
plot(BETA_vals,L2,'--','LineWidth',1.2);
plot(BETA_vals,L3,':','LineWidth',1.2);
xlabel('\beta'); ylabel('Lyapunov Exponents');
legend('\lambda_1','\lambda_2','\lambda_3'); grid on;

figure;
plot(D(:,1),D(:,2),'.','MarkerSize',3);
xlabel('\beta'); ylabel('x_{peak}');
grid on;
