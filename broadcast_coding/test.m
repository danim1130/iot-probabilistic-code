%clearvars

%x = sym('x',[1 9]);
%syms p;
%A = [0      p   (1-p)^2   p*(1-p)  0    0       0       0         0;
%     0      p   (1-p)^2   p*(1-p)  0    0       0       0         0;
%     0      p   (1-p)^2   p*(1-p)  0    0       0       0         0;        
%     0      0    0        0        p   (1-p)^2  p*(1-p) 0         0;
%     0      0    0        0        p   (1-p)^2  p*(1-p) 0         0;
%     0      0    0        0        p   (1-p)^2  p*(1-p) 0         0;
%     p^2    0    0        0        0    0       0      (1-p)^2    2*p*(1-p);
%     p^2    0    0        0        0    0       0      (1-p)^2    2*p*(1-p);
%     0      0    0        0        p   (1-p)^2  p*(1-p) 0         0];

x = sym('x',[1 17]);
syms p1;
syms p2;
%p2 = p1;
pl = (1-p1)*(1-p2);
A = [%  1           2           3           4           5           6           7           8           9           10         11           12          13            14           15           16           17
        0           p1          pl          (1-p1)*p2   0           0           0           0           0           0          0            0           0             0            0            0            0;          %m1 state (1)
        1           0           0           0           0           0           0           0           0           0          0            0           0             0            0            0            0;          %1 received (2)
        1           0           0           0           0           0           0           0           0           0          0            0           0             0            0            0            0;          %none received (3)
        0           0           0           0           1           0           0           0           0           0          0            0           0             0            0            0            0;          %2 received 1 (4)
        0           0           0           0           0           p2          pl          (1-p2)*p1   0           0          0            0           0             0            0            0            0;          %m2 state (5)
        0           0           0           0           1           0           0           0           0           0          0            0           0             0            0            0            0;          %2 received (6)
        0           0           0           0           1           0           0           0           0           0          0            0           0             0            0            0            0;          %none received (7)
        0           0           0           0           0           0           0           0           1           0          0            0           0             0            0            0            0;          %1 received 2 (8)
        0           0           0           0           0           0           0           0           0           pl         (1-p1)*p2    (1-p2)*p1   p1*p2         0            0            0            0;          %comb state (9)
        0           0           0           0           0           0           0           0           1           0          0            0           0             0            0            0            0;          %none received comb (10)
        0           0           0           0           1           0           0           0           0           0          0            0           0             0            0            0            0;          %2 received comb (11)
        0           0           0           0           0           0           0           0           0           0          0            0           0             1            0            0            0;          %1 received comb (12)
        1           0           0           0           0           0           0           0           0           0          0            0           0             0            0            0            0;          %both received comb (13)
        0           0           0           0           0           0           0           0           0           0          0            0           0             0            p1           pl           (1-p1)*p2;  %m1* state (14)
        0           0           0           0           0           0           0           0           0           0          0            0           0             1            0            0            0;          %1 received (15)
        0           0           0           0           0           0           0           0           0           0          0            0           0             1            0            0            0;          %none received (16)
        0           0           0           0           0           0           0           0           1           0          0            0           0             0            0            0            0;          %2 received 1 (17)
];

res = cell2sym(struct2cell(solve([x * A == x, sum(x) == 1], x)));
res_norm = res ./ sum(res);
received_msg_weight = [0 1 0 0 0 1 0 0 0 0 1 1 2 0 1 0 0];
total_msg_weight = [0 1 1 1 0 1 1 1 0 1 1 1 1 0 1 1 1];
goodput = simplify((received_msg_weight * res_norm) / (total_msg_weight * res_norm));
%state_values = [2 1 0 0 1 0 0 0 1];
%average_value = state_values * res_norm;
%figure;
%hold on;
%fplot(goodput, [0,1], 'k--');
%fplot(p1, [0, 1], 'k');
%xlabel('Transfer probability (p)');
%ylabel('Efficiency');
legend('1 to 2 routing', 'Simple routing', 'Location', 'southeast')
%fplot(1 - (1-p)^2, [0, 1]);
%figure;
%fplot(average_value / p, [0,1]);
%figure;
%fplot(average_value - p, [0,1]);