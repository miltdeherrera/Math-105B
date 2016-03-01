% Pt 1

% the function being estimated
syms x
f(x) = (1 + x^2)^-1;
fp(x) = diff(f);

% initializing the arrays

% step size
h1 = 10/11;
h2 = 10/21;
h3 = 10/51;
h4 = 10/101;

X1 = -5:h1:5;
X2 = -5:h2:5;
X3 = -5:h3:5;
X4 = -5:h4:5;

Y1 = double(f(X1));
Y2 = double(f(X2));
Y3 = double(f(X3));
Y4 = double(f(X4));






% setting up central differences arrays

for i = 1:numel(X1)
    buffer = double(f(X1(i)+h1*.5)-f(X1(i)-h1*.5));
    Y1cdiff(i) = buffer/h1;
end

for i = 1:numel(X2)
    buffer = double(f(X2(i)+h2*.5)-f(X2(i)-h2*.5));
    Y2cdiff(i) = buffer/h2;
end

for i = 1:numel(X3)
    buffer = double(f(X3(i)+h3*.5)-f(X3(i)-h3*.5));
    Y3cdiff(i) = buffer/h3;
end

for i = 1:numel(X4)
    buffer = double(f(X4(i)+h4*.5)-f(X4(i)-h4*.5));
    Y4cdiff(i) = buffer/h4;
end

% Calculating fwd differences
for i=1:numel(X1)-1
    Y1fdiff(i) = FwdDiff(Y1,i,h1);
end

for i=1:numel(X2)-1
    Y2fdiff(i) = FwdDiff(Y2,i,h2);
end

for i=1:numel(X3)-1
    Y3fdiff(i) = FwdDiff(Y3,i,h3);
end

for i=1:numel(X4)-1
    Y4fdiff(i) = FwdDiff(Y4,i,h4);
end

% Calculating backward differences
for i=2:numel(X1)
    Y1bdiff(i) = BckDiff(Y1,i,h1);
end

for i=2:numel(X2)
    Y2bdiff(i) = BckDiff(Y2,i,h2);
end

for i=2:numel(X3)
    Y3bdiff(i) = BckDiff(Y3,i,h3);
end

for i=2:numel(X4)
    Y4bdiff(i) = BckDiff(Y4,i,h4);
end




% Graphing the estimated derivatives


figure
title('Pt 1: Approximations of fprime(x)')
subplot(2,2,1)
plot(X1(1:numel(X1)-1),Y1fdiff)
hold on
plot(X1,Y1bdiff)
plot(X1,Y1cdiff)
title('Subplot 1: n = 11')
legend('Fwd diff','Back diff', 'Cdiff')

subplot(2,2,2)
plot(X2(1:numel(X2)-1),Y2fdiff)
hold on
plot(X2,Y2bdiff)
plot(X2,Y2cdiff)
title('Subplot 2: n = 21')
legend('Fwd diff','Back diff', 'Cdiff')

subplot(2,2,3)
plot(X3(1:numel(X3)-1),Y3fdiff)
hold on
plot(X3,Y3bdiff)
plot(X3,Y3cdiff)
title('Subplot 3: n = 51')
legend('Fwd diff','Back diff', 'Cdiff')

subplot(2,2,4)
plot(X4(1:numel(X4)-1),Y4fdiff)
hold on
plot(X4,Y4bdiff)
plot(X4,Y4cdiff)
title('Subplot 4: n = 101')
legend('Fwd diff','Back diff', 'Cdiff')


% Pt 2
% Discussion: As expected, the absolute error goes down as n increases. The
% step size becomes lower, increasing the denominator of the error term
% h*f''(xi)/2. The actual error starts out about half as larger as the
% theoretical error, and the ratio of actual error:theoretical error 
% becomes smaller as n becomes larger. The discrepancy seems to be caused
% by the (relatively) large value of M.

% setting up arrays with the values of the actual derivatives

Y1_actual = fp(X1);
Y2_actual = fp(X2);
Y3_actual = fp(X3);
Y4_actual = fp(X4);

% calculating the max absolute value of error

Y1_abs_error = max(abs(Y1_actual-Y1cdiff));
Y2_abs_error = max(abs(Y2_actual-Y2cdiff));
Y3_abs_error = max(abs(Y3_actual-Y3cdiff));
Y4_abs_error = max(abs(Y4_actual-Y4cdiff));

step_values = [11,21,51,101];
abs_error_array = [Y1_abs_error;Y2_abs_error;Y3_abs_error;Y4_abs_error];

% finding the max value of fpp (using same algorithm I used in HW4): I find
% the inflection points, calculate the value of f^2(x) for each one, and
% the max value equals f^2(xi)

fpp(x) = diff(fp);
fppp(x) = diff(fpp);

% finding the inflection points for f^3(x)
crit_points = double(solve(fppp));

size = numel(crit_points);

% finding the value of f^3(x) at each crit point to find the max value
% thereof
t = double(subs(fpp,crit_points(1)));

% since crit points all fall within [-5,5], in this case we don't need to implement
% additional criteria to make sure they aren't skewing the true value of M

for i = 2 : size-1
    y3rd = subs(fpp,crit_points(i));
    num3 = double(y3rd);
    
    if num3 > t
        t = num3;
    end
end

% t = M. Now we can output the final value of the error calculation

% calculating theoretical errors for each step size

th_error_array = [h1*t/2; h2*t/2; h3*t/2; h4*t/2];

C = horzcat(abs_error_array, th_error_array);


figure
title('Pt 2: error as n increases)')
bar(step_values,C)
xlabel('n')
ylabel('error estimate')
legend('Absolute Error','Theoretical Error')

% PT 3
% fpp(x) est := (1/h^2)f(x_0-h)-2f(x_0)+f(x_0+h)
% called "second derivative midpoint formula" in textbook

% fpp(x) was already constructed for pt 2.
% first will construct new central difference arrays. the difference
% between these and the ones we used for calculating central difference of
% pt. 1 is that there is a full step size between x_0 and the points we use
% to calculate it, vs half a step.

h1 = 10/11;
h2 = 10/21;
h3 = 10/51;
h4 = 10/101;

X1 = -5-h1:h1:5+h1;
X2 = -5-h2:h2:5+h2;
X3 = -5-h3:h3:5+h3;
X4 = -5-h4:h4:5+h4;

% Y values

Y1 = double(f(X1));
Y2 = double(f(X2));
Y3 = double(f(X3));
Y4 = double(f(X4));

% initialize second diff arrays

Y1seconddiff = zeros(numel(X1));
Y2seconddiff = zeros(numel(X2));
Y3seconddiff = zeros(numel(X3));
Y4seconddiff = zeros(numel(X4));

Y1seconddiff = Y1seconddiff(:,1);
Y2seconddiff = Y2seconddiff(:,1);
Y3seconddiff = Y3seconddiff(:,1);
Y4seconddiff = Y4seconddiff(:,1);


% we are only interested in the middle (n-2) values from each X array, so
% that will be implemented accordingly in the FOR loops

for i = 2:numel(X1)-1
    Y1seconddiff(i) = (1/h1)^2*((Y1(i-1) - 2*Y1(i) + Y1(i+1)));
end

for i = 2:numel(X2)-1
    Y2seconddiff(i) = (1/h2)^2*((Y2(i-1) - 2*Y2(i) + Y2(i+1)));
end

for i = 2:numel(X3)-1
    Y3seconddiff(i) = (1/h3)^2*((Y3(i-1) - 2*Y3(i) + Y3(i+1)));
end

for i = 2:numel(X4)-1
    Y4seconddiff(i) = (1/h4)^2*((Y4(i-1) - 2*Y4(i) + Y4(i+1)));
end

% discard first and last entries in each X and Y array

X1 = X1(2:end-1);
X2 = X2(2:end-1);
X3 = X3(2:end-1);
X4 = X4(2:end-1);

Y1seconddiff = Y1seconddiff(2:end-1);
Y2seconddiff = Y2seconddiff(2:end-1);
Y3seconddiff = Y3seconddiff(2:end-1);
Y4seconddiff = Y4seconddiff(2:end-1);

% Graphing -- basically the same as in pt 1



figure
title('Pt 3: Approximations of f^(2)(x)')
subplot(2,2,1)
plot(X1,Y1seconddiff)
title('Subplot 1: n = 11')

subplot(2,2,2)
plot(X2,Y2seconddiff)
hold on
title('Subplot 2: n = 21')

subplot(2,2,3)
plot(X3,Y3seconddiff)
title('Subplot 3: n = 51')

subplot(2,2,4)
plot(X4,Y4seconddiff)
title('Subplot 4: n = 101')


% PT 4
% Absolute error goes down somewhat as n increases, but gradually draws
% level with the theoretical error. This has to do with round-off error
% caused by small values of h, as mentioned in pg. 181 of the textbook.


% First, we find absolute errors.

Y1_actual = double(fpp(X1));
Y2_actual = double(fpp(X2));
Y3_actual = double(fpp(X3));
Y4_actual = double(fpp(X4));

Y1_abs_error = max(abs(Y1_actual-Y1seconddiff.'));
Y2_abs_error = max(abs(Y2_actual-Y2seconddiff.'));
Y3_abs_error = max(abs(Y3_actual-Y3seconddiff.'));
Y4_abs_error = max(abs(Y4_actual-Y4seconddiff.'));

step_values = [11,21,51,101];
abs_error_array = [Y1_abs_error;Y2_abs_error;Y3_abs_error;Y4_abs_error];

% Now we find theoretical error limits

% finding the max value of fpppp (using same algorithm I used in HW4)
fpp(x) = diff(fp);
fppp(x) = diff(fpp);
fpppp(x) = diff(fppp);
fppppp(x) = diff(fpppp);

% finding the inflection points for f^5(x)
crit_points_tentative = double(solve(fppppp));

size = numel(crit_points_tentative);

% getting rid of crit points outside of [-5,5]

crit_points = [];


for i = 1:size
    if crit_points_tentative(i) >= -5 & crit_points_tentative(i) <= 5
        crit_points = [crit_points crit_points_tentative(i)];
    end
end

% finding the value of f^4(x) at each crit point to find the max value

t = double(subs(fpppp,crit_points(1)));

for i = 2 : size-1
    y4th = subs(fpppp,crit_points(i));
    num4 = double(y4th);
    
    if num4 > t
        x_true = crit_points(i);
    end
end

% t = M. Now we can output the final value of the error calculation

% calculating theoretical errors for each step size

th_error_array = [h1^2*t/12; h2^2*t/12; h3^2*t/12; h4^2*t/12];

C = horzcat(abs_error_array, th_error_array);


figure
title('Pt 4: error as n increases)')
bar(step_values,C)
xlabel('n')
ylabel('error estimate')
legend('Absolute Error','Theoretical Error')


% Pt 5
% As N goes up, the error becomes small to the point of only being visible
% on logarithmic axes. 

% Let g(t) = [f(10^-t) - f(0)]/10^-t. This equates to g(t) =
% [(1+10^-2t)^(-1)-1]*10^t. We end up with:
%
%   g(t) = -10^t/(10^2t+1)
%
% lim g(t) as t -> +inf = 0. Thus we would expect g(t) to get very small as
% t becomes large, and in fact this is what we observe.

N = [1 2 5 10 20 40];

for i = 1:6
    NDiff(i) = (f(10^-N(i)) - f(0))/10^-N(i);
end

figure
title('Pt 5: scatter plot of N vs estimate of fp(x)')
scatter(N,NDiff)
set(gca,'yscale','log')
xlabel('N')
ylabel('estimate of d/dx*[f(0)]')
