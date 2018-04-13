function [mu, sigma2] = cgp_test()

    X = [0,0.02,0.075,0.08,0.14,0.15,0.155,0.156,0.18,0.22,0.29,0.32,0.36,...
        0.37,0.42,0.5,0.57,0.63,0.72,0.785,0.8,0.84,0.925,1;...
        0.29,0.02,0.12,0.58,0.38,0.87,0.01,0.12,0.22,0.08,0.34,0.185,0.64,...
        0.02,0.93,0.15,0.42,0.71,1,0,0.21,0.5,0.785,0.21]';
    y = [0.33699176 -0.65759451  0.89324920 -0.74273446 -0.96290377  0.41256824 ...
        0.99282109  0.08077027 -0.89687895 -0.09583807 -0.52584450 -0.96401324 ...
        0.65914050 -0.55305251  0.98017623 -0.61080573  0.66921334  0.99251078 ...
        0.94700532 -0.70490141  0.51436857  0.98695415  0.94655705  0.78608913]';
    x_hats = [0,0.2; ...
              0,0.289;
              0.9251,0.786];
 
    X = rand([30,4]);
    y = rand([30,1]);
    x_hats = rand([100000,4]);

    Sig = 0.1*diag(ones([size(X,1),1]));
    v = 0.1*ones(size(x_hats,1),1);
    beta = 0;
    lambda_lower = 0.2;
    [mu, sigma2] = composite_gp_estimation(X,y,x_hats,Sig,v,beta,lambda_lower);
    if size(x_hats,1) < 10000
        csvwrite('tmp_X.csv',X);
        csvwrite('tmp_y.csv',y);
        csvwrite('tmp_x_hats.csv',x_hats);
        csvwrite('tmp_Sig.csv',Sig);
        csvwrite('tmp_v.csv',v);
        csvwrite('tmp_beta.csv',beta);
        csvwrite('tmp_lambda_lower.csv',lambda_lower);

        system('rm -f tmp_Yp.csv tmp_Varp.csv');
        system('/usr/local/bin/R CMD BATCH ../run_cgp.R');
        mu_R = csvread('tmp_Yp.csv',1,0)
        sigma2_R = csvread('tmp_Varp.csv',1,0)
        system('rm -f tmp_*.csv');

        fprintf('max differences:\n');
        max(mu-mu_R)
        max(sigma2-sigma2_R)
    end
end
