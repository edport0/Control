%Script da matéria de controle e aplicações, continuação de
%modelagem de sistemas dinâmicos.

%comecemos escrevendo os valores dos nossos parâmetros

%parâmetros
B_b = 2.5;
D_b = 0.4;
Loa = 8.4;
Rho_w = 997;
H_b = 1.7;
g = 9.81;
M_d = 150;
M_e = 10;
r_d = 0.230;
r_e = 0.090;
e_d = 0.150;
l_e = 0.500;
c_1 = 0.428;
c_2 = 0.006;
c_3 = 0.149;
w = 0.50/(sqrt(2.5/(2*g)));
Omega = 8000*pi/30;
C_g = 500;
B_g = 2500;
a = 2.0;


A = [
    0, 0, 1, 0;
    0, 0, 0, 1;
    -g*Loa*(B_b^2-6*D_b^2)/(D_b*(B_b^2*(Loa+12*c_3)+Loa*H_b^2)), 0, -(12*sqrt(2)*g^2*(B_b/g)^(3/2)*c_2)/(B_b^2*(Loa+12*c_3)+Loa*H_b^2), -(6*Omega*M_d*r_d^2)/(B_b*D_b*Rho_w*(B_b^2*(Loa+12*c_3)+Loa*H_b^2));
    0, -12*C_g/(e_d^2*(M_d-M_e)+3*M_d*r_d^2+M_e*(l_e^2+3*r_e^2)), Omega*(-1+(12*M_d*r_d^2)/(e_d^2*(M_d-M_e)+3*M_d*r_d^2+M_e*(l_e^2+3*r_e^2))), -(12*B_g)/(e_d^2*(M_d-M_e)+3*M_d*r_d^2+M_e*(l_e^2+3*r_e^2))
];


E = [
    0;
    0;
    g * Loa * (B_b^2 - 6 * D_b^2) / (D_b * (B_b^2 * (Loa + 12 * c_3) + Loa * H_b^2));
    0;
];

C = [
    1,0,0,0;
    0,1,0,0;
    0,0,0,1;
];

B = [
    0;
    0;
    0;
    1/(1/12*M_e*(3*r_e^2+l_e^2-e_d^2)+1/12*M_d*(3*r_d^2+e_d^2))
];
% B = [
%     0;
%     0;
%     0;
%     1/(2*(e_d^2*(M_d-M_e)+3*M_d*r_d^2+M_e*(l_e^2+3*r_e^2)));
% ];

D = 0;

%após definir as matrizes para o espaço de estados, podemos buscar as
%informações da matriz de controlabilidade e observabbilidade.
%a partir dessas matrizes conseguimos conluir se o sistema  é ou não
%controlável e observável.

CM = ctrb(A,B);

if (rank(CM) == size(A,1))
    disp("Modelo Controlável")

end

OB = obsv(A,C);

if (rank(OB) == size(A,1))
    disp("Modelo Observável")

end


%feito isso, vamos buscar entender mais sobre o problema do regulador(2
%tipos - alocação de polos e lqr) comecemos pelo de alocação.

sys = ss(A,B,C,D); %passa para espaço de estados
polos = pole(sys); %calcula os polos do sistema em malha aberta
p = poly(A); %me dá o polinomio caracteristico da matriz A


% Escolher polos desejados para malha fechada
polos_desejados = [-2.4,-4.7601+5*0.346198i,-4.7601000-5*0.346198i,-2.5]; %polos desejadas conforme Ogata[-1903,-0.47601+5*0.346198i,-0.47601000-5*0.346198i,-23161]

% referencia:[-1.903,-4.7601+5*0.346198i,-4.7601000-5*0.346198i,-2.3161];
%para voltar ao oscilatorio dividir a parte real dos num complex por 10


delta_t = 0.01;
ti = 0;
tf = 100;
t_array = ti:delta_t:tf;
x_0 = [(10*pi/180) 0 0 0];
x_a = zeros(size(A, 2), length(t_array));
u_a = zeros(size(B, 2), length(t_array));
x_a(:, 1) = x_0; % Certifique-se de que x_0 está definido
% Projetar a matriz de ganho K
K = place(A, B, polos_desejados);

A_ = A - B * K;

for i = 1:(length(t_array) - 1)
    x_a(:, i + 1) = expm(A_ * delta_t) * x_a(:, i);
    u_a(:, i) = -K * x_a(:, i);
end

fig = 1; % Inicializar a variável fig
% for i = 1:1:size(x_a,1)
%     figure(fig) % Certifique-se de que isso está sendo usado como pretendido
%     plot(t_array, x_a(i,:), 'r', 'LineWidth', 2)
%     title('Alocação de Polos')
%     xlim([0 10])
%     xticks(0:5:40)
%     set(gca, 'FontSize', 22)
%     fig = fig + 1; % Incrementar fig para a próxima figura
% 
%     if i == 1
%         xlabel('Tempo [s]')
%         ylabel('Rolagem \phi(t) [rad]')
%     elseif i == 2
%         xlabel('Tempo [s]')
%         ylabel('Nutação \theta(t) [rad]')
%     elseif i == 3
%         xlabel('Tempo [s]')
%         ylabel('Velocidade de roll [rad/s] ')
%     elseif i == 4
%         xlabel('Tempo [s]')
%         ylabel('Velocidade de nutação [rad/s]', 'Interpreter', 'latex')
%     end
% end
   
poly_aloc = poly(A_); %polinomio caracteristico malha fechada
sys_a = ss(A_,B,C,D);
polos_aloc = pole(sys_a); %polos malha fechada


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%agora vamos fazer LQR
%a matriz Q deve ser 4x4 uma vez que nosso vetor de estados tem 4
%componentes
%R é definido com um único valor pois nossa entrada é apenas um "momento de
%controle"
%as diagonais são escolhidas conforme o método de Bryson, q_i = 1/x_i_max^2
%e u_i = 1/u_i_max^2
%precisamos definir roll max, nutação max, e suas respectivas velocidades
%os dados vão ser colocados em rad ou rad/s (posição e velocidades angulare
%na matriz Q)



% Define as matrizes de ponderação Q e R para o LQR
R =  0.5;
Q = [ 1, 0, 0, 0; %17566.8
      0, 1, 0, 0;   %905.2
      0, 0, 100000000, 0;   % 100000000
      0, 0, 0, 100000];  %68.8

% Calcula a matriz de ganho do controlador LQR
K_lqr = lqr(A,B,Q,R);

% Define o novo sistema A para o LQR
F_lqr = A - B*K_lqr;

% Inicializa o vetor de estado para o LQR
x_lqr = zeros(size(A,2), length(t_array));
x_lqr(:,1) = x_0;

% Inicializa o vetor de controle para o LQR
u_lqr = zeros(size(B,2), length(t_array));


for i = 1:(length(t_array) - 1)
    x_lqr(:, i + 1) = expm(F_lqr * delta_t) * x_lqr(:, i);
    u_lqr(:, i) = -K * x_lqr(:, i);
end


% Cria um sistema de espaço de estados para o controlador LQR
sys_lqr = ss(F_lqr, B, C, D);
polos_lqr = pole(sys_lqr);


figure(fig)
pzmap(sys_lqr,'r',sys_a,'b')
legend('LQR','Aloca o')
fig = fig +1;

%analise grafica

cor_lqr = 'r'; % Cor vermelha para LQR
cor_alocacao = 'b'; % Cor azul para Alocação

for i = 1:size(x_lqr, 1)
    figure(fig)
    plot(t_array, x_lqr(i,:), cor_lqr, 'LineWidth', 2)
    hold on
    plot(t_array, x_a(i,:), cor_alocacao, 'LineWidth', 2)
    title('LQR')
    xlim([0 20])
    xticks(0:2:20)
    legend('LQR','Alocação','Location','best')
    set(gca, 'FontSize', 22)
    fig = fig + 1;

    if i == 1
        xlabel('Tempo [s]')
        ylabel('Rolagem \phi(t) [rad]')
    elseif i == 2
        xlabel('Tempo [s]')
        ylabel('Nutação \theta(t) [rad]')
    elseif i == 3
        xlabel('Tempo [s]')
        ylabel('Velocidade de roll [rad/s]')
    elseif i == 4
        xlabel('Tempo [s]')
        ylabel('Velocidade de nutação [rad/s]', 'Interpreter', 'latex')
    end
    hold off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bode(sys_a);
% bode(sys_lqr);

%%%%%%%%% vamos analisar as entradas para sabermos qual método


for i=1:size(u_lqr,1)
    figure(fig)
    plot(t_array, u_lqr(i,:),'m','LineWidth', 3 );
    hold on
    plot(t_array, u_a(i,:),'g','LineWidth', 3 );
    %title('LQR')
    xlim([0 10])
    xticks([0:2:10])
    legend('LQR',' Alocação ')
    set(gca,'FontSize',22)
    fig = fig + 1;
    if i==1
    xlabel('Tempo [s]')
    ylabel('Entrada $u1$ [Nm]','Interpreter','latex')
    end
end 

%----------------------------------------------------------------------------------


% OBSERVADOR

%ALOCAÇÃO OBSERVADOR
% Polos desejados para o observador (ajustados para serem o dobro da parte real dos polos originais)
polos_desejados_obs = [-4.8, -9.5202 + 5*0.346198i, -9.5202 - 5*0.346198i, -5.0];

% Projetar a matriz de ganho L para o observador
L = place(A', C', polos_desejados_obs)'; % Correção aqui, A' e C' no place e transpor o resultado para obter L

% Define a nova matriz A para o observador
A_o = A - L * C;

% Simulação de observador
% Inicializa o vetor de estado estimado
x_est = zeros(size(A, 2), length(t_array));
x_est(:, 1) = x_0; % Considerando que você deseja começar com a mesma condição inicial

% Loop de simulação para o observador
for i = 1:(length(t_array) - 1)
    % Simulação do próximo estado com estimativa corrente
    x_est(:, i + 1) = expm(A_o * delta_t) * x_est(:, i) + L * (C * x_a(:, i) - C * x_est(:, i));
end

% Calcula polinômio característico da matriz A_o e os polos do observador
poly_char_obs = poly(A_o); 
sys_obs = ss(A_o, B, C, D);
polos_obs = pole(sys_obs);

% Plot dos polos do sistema controlador e observador
figure(fig)
pzmap(sys_obs, 'r', sys_a, 'b')
legend('Observador', 'Controlador')
title('Mapa de Polos: Observador e Controlador por Alocação')
fig = fig + 1;

%---------------------------------------------------------------------------------
%LQR OBSERVADOR - ISSO AQUI VAI DAR UM BELO TRABALHO


%REFERENCIAS (Q E R USADAS PARA O CONTROLADOR)
% R =  0.5;
% Q = [ 1, 0, 0, 0;
%       0, 1, 0, 0;   
%       0, 0, 100000000, 0;   
%       0, 0, 0, 100000];  

R_o = 100;

Q_o = [ 1000, 0, 0, 0; 
      0, 1000, 0, 0;   
      0, 0, 1000, 0;   
      0, 0, 0, 1000];  

L_lqr = lqr(A',C',Q_o,R_o);

% Define o novo sistema A para o LQR
A_lqr = A - L_lqr'*C;

sys_lqro=ss(A_lqr,B,C,D);
polos_lqro = pole(sys_lqro);


% Plot dos polos do sistema controlador e observador
figure(fig)
pzmap(sys_lqro, 'r', sys_a, 'b')
legend('Observador', 'Controlador')
title('Mapa de Polos: Observador e Controlador por LQR')
fig = fig + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Erros de observação Inicial
% Valores iniciais dos erros
e1 = 0.1; %0.1
e2 = 0; %0
e3 = 0.1/180; %0.1/180  
e4 = 0.3; %0.3

% Definição inicial do vetor de erro
e_o0 = [e1; e2; e3; e4];  % Vetor coluna

% Alocação do espaço para o erro
e_ao = zeros(size(e_o0, 1), length(t_array));  % Alocação de espaço para erro de alocação de polos
e_lqo = zeros(size(e_o0, 1), length(t_array)); %mesma coisa para LQR
% Condição inicial de erro
e_ao(:, 1) = e_o0;  % Definindo a condição inicial para erro de alocação de polos
e_lqo(:,1) = e_o0;


% Simulação dos erros ao longo do tempo
for i = 1:floor(tf / delta_t)  % Uso de floor para garantir correta terminação
    e_ao(:, i + 1) = expm((A - L * C) * delta_t) * e_ao(:, i);  % Propagação do erro para alocação de polos
    e_lqo(:, i + 1) = expm((A - L_lqr' * C) * delta_t) * e_lqo(:, i);
end

% Análise gráfica do erro de observação controlado por observador por alocação de polos
figure(fig);
for j = 1:size(e_ao, 1)  % Loop para plotar todos os estados
    plot(t_array, e_ao(j, :), 'LineWidth', 2);
    hold on;
end
title('Erro do Observador por Alocação de Polos');
xlabel('Tempo [s]');
ylabel('Erro do observador e(t)');
xlim([0, 10]);
xticks(0:1:10);
legend({'e_1 [rad]', 'e_2 [rad]', 'e_3 [rad/s]', 'e_4 [rad/s]'}, 'Location', 'northeast');  % Atualizado para refletir todos os estados
set(gca, 'FontSize', 22);
hold off;
fig = fig + 1;  % Incrementar o contador de figura

%---------------------------------------------------------------------------------
%ANALOGO PARA LQR ASSIM QUE FAZER A PARTE INCIAL DE LQR

figure(fig);
for j = 1:size(e_lqo, 1)  % Loop para plotar todos os estados
    plot(t_array, e_lqo(j, :), 'LineWidth', 2);
    hold on;
end
title('Erro do Observador por LQR');
xlabel('Tempo [s]');
ylabel('Erro do observador e(t)');
xlim([0, 10]);
xticks(0:1:10);
legend({'e_1 [rad]', 'e_2 [rad]', 'e_3 [rad/s]', 'e_4 [rad/s]'}, 'Location', 'northeast');  % Atualizado para refletir todos os estados
set(gca, 'FontSize', 22);
hold off;
fig = fig + 1;  % Incrementar o contador de figura

%--------------------------------------------------------------------------------
% PRINCIPIO DA SEPARAÇÃO
%--------------------------------------------------------------------------------

%ALOCAÇÃO

K_a = K;
K_o = L;
K_ao = L;

Lambda_a = [A, -B*K_a; K_ao*C, A-B*K_a-K_ao*C];
Lambda_lq = [A, -B*K_lqr; L_lqr'*C, A-B*K_lqr-L_lqr'*C];

Phi_ao = expm(Lambda_a*delta_t); %matriz de ttransição para alocação de polos
Phi_lqo = expm(Lambda_lq*delta_t);

x_chap0= x_0' - e_o0; %estado observador inicial

z_a = zeros(length(x_0) + length(x_chap0),length(t_array));
z_lq = zeros(length(x_0) + length(x_chap0),length(t_array));

z_a(:,1) = [x_0';x_chap0];
z_lq(:,1) = [x_0';x_chap0];

for i=1:tf/delta_t
    z_a(:,i+1) = Phi_ao*z_a(:,i);
    z_lq(:,i+1) = Phi_lqo*z_lq(:,1);
end

% Alocação do espaço para o erro (princípio da separação)
e_ao_ps = zeros(size(z_a, 1) / 2, length(t_array));  % Alocação de espaço para erro

% Cálculo do erro baseado no princípio da separação
e_ao_ps = z_a(1:size(z_a, 1) / 2, :) - z_a(size(z_a, 1) / 2 + 1:size(z_a, 1), :);

e_lqo_ps= zeros(size(z_lq, 1) / 2, length(t_array));

e_lqo_ps = z_lq(1:size(z_lq, 1) / 2, :) - z_lq(size(z_lq, 1) / 2 + 1:size(z_lq, 1), :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Observador por Alocação (Posições)
for i=5:6
figure(fig)
plot(t_array,z_a(i,:),'LineWidth',2)
hold on
title('0bservador por Alocação')
xlabel("Tempo (s)")
ylabel("Estado (rad)")
xlim([0 20])
xticks([0:5:20])
legend('Phi','Theta')
set(gca,'FontSize',22)
end
hold off
fig = fig +1;

for i=1:2
figure(fig)
plot(t_array,e_ao(i,:),'LineWidth',2)
hold on
title('ErroObs - Alocação')
xlabel("Tempo (s)")
ylabel("Estado (rad)")
xlim([0 20])
xticks([0:5:20])
legend('Phi','Theta')
set(gca,'FontSize',22)
end
hold off
fig = fig +1;

%-----------------------------------------------------------------
%Observador por Alocação (Derivada das Posições)
for i=7:8
figure(fig)
plot(t_array,z_a(i,:),'LineWidth',2)
hold on
title('0bservador por Alocação')
xlabel("Tempo (s)")
ylabel("Estado (rad/s)")
xlim([0 20])
xticks([0:5:20])
legend('DPhi/dt','DTheta/dt')
set(gca,'FontSize',22)
end
hold off
fig = fig +1;

for i=3:4
figure(fig)
plot(t_array,e_ao(i,:),'LineWidth',2)
hold on
title('ErroObs - Alocação')
xlabel("Tempo (s)")
ylabel("Estado (rad/s)")
xlim([0 20])
xticks([0:5:20])
legend('DPhi/dt','DTheta/dt')
set(gca,'FontSize',22)
end
hold off
fig = fig +1;

%-------------------------------------------------------------------
%COMPARAÇÃO ESTADO REAL E ESTADO OBSERVADO (ALOCAÇÃO)
for i = 1:size(x_a, 1)
    figure(fig)
    plot(t_array, z_a(i+4,:), cor_lqr, 'LineWidth', 2)
    hold on
    plot(t_array, x_a(i,:), cor_alocacao, 'LineWidth', 2)
    title('Comparativo Alocação')
    xlim([0 20])
    xticks(0:2:20)
    legend('Observado','Real','Location','best')
    set(gca, 'FontSize', 22)
    fig = fig + 1;

    if i == 1
        xlabel('Tempo [s]')
        ylabel('Rolagem \phi(t) [rad]')
    elseif i == 2
        xlabel('Tempo [s]')
        ylabel('Nutação \theta(t) [rad]')
    elseif i == 3
        xlabel('Tempo [s]')
        ylabel('Velocidade de roll [rad/s]')
    elseif i == 4
        xlabel('Tempo [s]')
        ylabel('Velocidade de nutação [rad/s]', 'Interpreter', 'latex')
    end
    hold off
end

%AJEITAR TODA LÓGICA DAQUI PARA BAIXO!!!!!!!!!!!!!!!!!!!!!!!!
%-----------------------------------------------------------------------
% Observador por LQR (Posições)

% for i=5:6
% figure(fig)
% plot(t_array,z_lq(i,:),'LineWidth',2)
% hold on
% title('0bservador por LQR')
% xlabel("Tempo (s)")
% ylabel("Estado (rad)")
% xlim([0 20])
% xticks([0:5:20])
% legend('Phi','Theta')
% set(gca,'FontSize',22)
% end
% hold off
% fig = fig +1;
% 
% for i=1:2
% figure(fig)
% plot(t_array,e_lqo(i,:),'LineWidth',2)
% hold on
% title('ErroObs - LQR')
% xlabel("Tempo (s)")
% ylabel("Estado (rad)")
% xlim([0 20])
% xticks([0:5:20])
% legend('Phi','Theta')
% set(gca,'FontSize',22)
% end
% hold off
% fig = fig +1;
% 
% 
% %Observador por LQR (Derivada das Posições)
% for i=7:8
% figure(fig)
% plot(t_array,z_lq(i,:),'LineWidth',2)
% hold on
% title('0bservador por LQR')
% xlabel("Tempo (s)")
% ylabel("Estado (rad/s)")
% xlim([0 20])
% xticks([0:5:20])
% legend('DPhi/dt','DTheta/dt')
% set(gca,'FontSize',22)
% end
% hold off
% fig = fig +1;
% 
% for i=3:4
% figure(fig)
% plot(t_array,e_lqo(i,:),'LineWidth',2)
% hold on
% title('ErroObs - LQR')
% xlabel("Tempo (s)")
% ylabel("Estado (rad/s)")
% xlim([0 20])
% xticks([0:5:20])
% legend('DPhi/dt','DTheta/dt')
% set(gca,'FontSize',22)
% end
% hold off
% fig = fig +1;
