function [fit_result] = PointByPointExponentialFit(data,T)
%assume data is of form: x,y,z, and we're fitting along z.
dimX = size(data,1);
dimY = size(data,2);
fit_result = 0*data;

for k=1:dimX;
    for j=1:dimY;
              [~,temp,~] = single_exp_fit(T', squeeze(data(k,j,:))');
              fit_result(k,j,:) = temp;
    end
end
