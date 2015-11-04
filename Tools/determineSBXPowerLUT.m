

%% Collect data
x = [0,1,2,8:8:224,108,110,111,196,204,206,203,205];
y = [.242,.2417,.2419,.544,2.35,4.25,5.76,16.91,24.8,29.56,34.5,39.41,43.94,48.25,52.15,55.56,112.3,123.8,134.4,145.0,154.5,163.1,170.9,177.2,183.1,186.9,189.7,191.7,191.3,190.7,189.1,57.2,57.99,111.2,190.6,191.4,191.6,191.9,192.2];


%% Interpolate curve

% interpolate curve
[~,ind] = max(y);
DACregistry = 0:x(ind);
values = interp1(x,y,DACregistry);

% check interpolation
[~,order] = sort(x);
figure;
plot(x(order),y(order), 'bx');
hold on;
plot(DACregistry, values, 'r--');
legend('Data', 'Interpolation','Location','SouthEast');
xlabel('DAC registry value');
ylabel('Power (mW)');


%% Create LUT
N = 256;

ideal = min(y):range(y)/(N-1):max(y);

LUT = nan(1,256);
for index = 1:N
    [~,ind] = min(abs(ideal(index)-values));
    LUT(index) = DACregistry(ind);
end

