%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATING THE INTERPOLATION WINDOW %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Calculates the interpolation window G according to J. Selva,
% "Interpolation of Bounded Bandlimited Signals and Applications", see eq. (12).
% The interpolation window is the convolution, in the frequency domain, of
% G0(f) and G1(f). G0(f) is a rectangular window of bandwidth B0, while G1
% is an approximation of the prolate spheroidal wave function.
% The interpolation window G is the same appearing in eq. (25) of J. Selva,
% "Convolution-Based Trigonometric Interpolation of Band-Limited Signals".
function [G, kg] = KnabWindow(T, B, P, N)

% --- B:        bilateral bandwidth of the signal to be interpolated

% --- Differential frequency
df = 1/(N*T);

% --- kg factor
kg = floor((1/T-B/2)/df);   

% --- Frequencies where calculating G = G0 * G1;
f = -df * kg + df * (0: kg);

B0 = 1/T;                           % --- Bilateral bandwidth of the sinc function
B1 = B0 - B;                        % --- Bilateral bandwidth of the approximate prolate spheroidal wave function
I = (B1 / 2 <= (f + B0 / 2));

% --- Definition of G
G = zeros(2*size(f)-1);

% --- Forces G=1 in the signal bandwidth
G(I)=1;

Int = @(x)ConvWindow(x, B1, 2 * P * T);

for k = 1:find(not(I), 1, 'last' );

   IInt = IntervalIntersection([-B1/2, B1/2], f(k) + [-B0 / 2, B0 / 2]);                
   if ~isempty(IInt)
       G(k) = quad(Int, IInt(1), IInt(2), 1e-12);
   end
   
end

G = G + flipdim(G, 2);
G(kg + 1) = G(kg + 1) - 1;
G = G / N;

