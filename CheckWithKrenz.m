% Passive wall
% elastin

col2=jet(Newgen);
figure;
for k=1:Newgen
    
    kc = mechanical_properties_PA(Pmid(k),Radius(k));
    
    
    %     figure;hold on;
    
    PP = [0:Pmid(k)/10:45*133.32];
    D0 = [];
    for j=1:length(PP)
        D0(j) = 2*ZeroP(Me(k),Mt(k),Radius(k),Pmid(k),PP(j));
    end
    scatter(D0./D0(1),PP/133.32,'MarkerFaceColor',col2(k,:));hold on;
end
plot(1 + 0.012*PP/133.32,PP/133.32);

load('data_loaded_2.mat');

data = data_fit(:,2);
scatter(data,data_fit(:,1));

% axis([1 D0(end)./D0(1) 0 PP(end)/133.32])
% figure;scatter(2*Radius./D0,Pmid/133.32);
% hold on;
% scatter(2*Radius./D0,1+0.02*Pmid/133.32);
