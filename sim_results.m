% Benchmark Script: SNG vs RSC
numRuns = 10;
res_sng = []; % [time, dist, nAll, nRoute]
res_rsc = [];

fprintf('Benchmark baslatildi. 10 basarili run hedefleniyor...\n');

while size(res_sng, 1) < numRuns || size(res_rsc, 1) < numRuns
    currentSeed = randi(100000);
    
    % --- SNG TESTİ ---
    if size(res_sng, 1) < numRuns
        rng(currentSeed);
        tic;
        try
            run('SNG.m'); % Dosya adinin SNG.m oldugundan emin ol
            t = toc;
            if exist('pathIds','var') && ~isempty(pathIds) && exist('route_points','var')
                d = 0;
                for k = 1:size(route_points, 1)-1
                    d = d + norm(route_points(k+1, :) - route_points(k, :));
                end
                res_sng(end+1, :) = [t, d, length(nodes), length(pathIds)];
                fprintf('SNG Basarili (%d/%d)\n', size(res_sng,1), numRuns);
            end
        catch
            % Hata durumunda sessizce devam et
        end
        clearvars -except res_sng res_rsc numRuns currentSeed % Verileri koru
    end

    % --- RSC TESTİ ---
    if size(res_rsc, 1) < numRuns
        rng(currentSeed);
        tic;
        try
            run('RSC.m');
            t = toc;
            if exist('pathIds','var') && ~isempty(pathIds) && exist('route_points','var')
                d = 0;
                for k = 1:size(route_points, 1)-1
                    d = d + norm(route_points(k+1, :) - route_points(k, :));
                end
                res_rsc(end+1, :) = [t, d, length(nodes), length(pathIds)];
                fprintf('RSC Basarili (%d/%d)\n', size(res_rsc,1), numRuns);
            end
        catch
            % Hata durumunda sessizce devam et
        end
        clearvars -except res_sng res_rsc numRuns currentSeed
    end
end

% Ortalamalari Hesapla
avg_sng = mean(res_sng, 1);
avg_rsc = mean(res_rsc, 1);

% Nihai Tablo Çıktısı
fprintf('\n============================================================\n');
fprintf('                 10 TEST ORTALAMA SONUÇLARI\n');
fprintf('============================================================\n');
fprintf('%-15s | %-10s | %-10s | %-10s | %-10s\n', 'Algoritma', 'Süre(s)', 'Mesafe(m)', 'Düğüm(Top)', 'Düğüm(Yol)');
fprintf('------------------------------------------------------------\n');
fprintf('%-15s | %-10.4f | %-10.2f | %-10.1f | %-10.1f\n', 'Dairesel SNG', avg_sng(1), avg_sng(2), avg_sng(3), avg_sng(4));
fprintf('%-15s | %-10.4f | %-10.2f | %-10.1f | %-10.1f\n', 'Dairesel RSC', avg_rsc(1), avg_rsc(2), avg_rsc(3), avg_rsc(4));
fprintf('============================================================\n');