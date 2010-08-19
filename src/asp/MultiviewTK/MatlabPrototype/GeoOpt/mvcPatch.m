% read multiple-view data 
Mviews=mvread();
gDEM = imread('../../ground-DEM.tif');
sDEM = imread('../../initial-DEM.tif');

Patch.center=[500; 500];
Patch.width=[30 6]; % correlation and smoothing window size
Patch.height=double(gDEM(512, 512))+10;
Patch.range=Patch.height+[-60 60];     % crop the region of interest
Patch.fill = 1e6;   % scale of FillValues for imtransform
Patch.plane=[pix2dir(Mviews.georef,Patch.center); ...
    Mviews.radius+Patch.height];

Patch=mpcrop(Mviews,Patch);
Patch=mporthoproj(Patch);

p = Patch.plane;
mpgeoGaus(p,Patch);

options = optimset('LargeScale','on','Display','iter');
options = optimset(options,'GradObj','off','GradConstr','on');
options = optimset(options,'TolFun',1e-9,'TolX',1e-9);
options = optimset(options,'MaxFunEvals',120000,'MaxIter',150);
options = optimset(options,'MaxTime',600);

before = now;2
datestr(before)

[p,fval,exitflag,output] = fmincon(@(p)mpgeoGaus(p,Patch),p,[],[],[],[], ...
    p-[0.5*[1;1;1]; 100],p+[0.5*[1;1;1]; 100],@(p)UnitNorm(p),options);
%    p-[0.005*[1;1;1]; 100],p+[0.005*[1;1;1]; 100],@(p)UnitNorm(p),options);

after = now;
datestr(before)
datestr(after)
(after-before)*24*3600

oPatch=Patch;
oPatch.plane=p;
oPatch.height=p(4)/(p(1:3)'*pix2dir(Mviews.georef,Patch.center))-Mviews.radius
oPatch=mporthoproj(oPatch);
for i=1:numel(oPatch.ortho)
     figure(1), subplot(2,2,i), imshow(oPatch.image{i});
     figure(2), subplot(2,2,i), imshow(oPatch.ortho{i});
end
e = pix2dir(Mviews.georef,Patch.center);
p(4)/(p(1:3)'*e)-Mviews.radius
double(gDEM(512, 512))