function Patch = mporthoproj(Patch)
% Patch.plane: plane parmeters Patch.plane(1:3) normal vector Patch.plane(4) distance from the origin
% I: Orbital images
% P: Projection matrices
% H: Homgraphy transformation
% Patch.center: center point [x; y]

% normalized plane
p = Patch.plane/norm(Patch.plane(1:3));
v = p(1:3); % normal vector
d = p(4);   % distance from the origin

% direction vector and its Jacobian
[e, de] = pix2dir(Patch.georef,Patch.center);

n=numel(Patch.camera);
for i=1:n
    Q = d*Patch.camera{i}(:,1:3)+Patch.camera{i}(:,4)*v';
    u = Q*e;
    S = Q*de*Patch.georef(1:2,1:2);

    % homography
    Hp = [S u-S*Patch.center];

    tform = maketform('projective',inv(Hp)');
    Patch.ortho{i} = imtransform(Patch.image{i},tform,'bicubic',...
        'xdata',Patch.center(1)+Patch.width(1)*[-1 1]/2,...
        'ydata',Patch.center(2)+Patch.width(1)*[-1 1]/2,...
        'FillValues',Patch.fill*i);
end