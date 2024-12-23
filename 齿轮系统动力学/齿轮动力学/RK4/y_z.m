function zhuanjiao = y_z()
%% 压力角导入
alpha = readmatrix('K_betai_K_M.csv', 'Range', 'A2:A1002');
%% alpha = 20 * ones(301, 1) + noise_level * randn(1, 301)';
zhuanjiao=tand(alpha)-tand(11.4950);
end

