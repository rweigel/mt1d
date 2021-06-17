% Pixel coordinates derived using addDigitizer.jy in
% Autoplot on the image locations.png from USGS.
% http://geomag.usgs.gov/conductivity/

clear;more off;
pixelsreference = load('location_grid_pixels.txt');

regions = {'CP1_Coordinates_Pixel.txt',...
           'PT1_Coordinates_Pixel.txt'};

pixelsregion = load(regions{2});

F = TriScatteredInterp(pixelsreference(:,1),pixelsreference(:,2),pixelsreference(:,3));
lats = F(pixelsregion(:,1),pixelsregion(:,2));
F = TriScatteredInterp(pixelsreference(:,1),pixelsreference(:,2),pixelsreference(:,4));
lons = F(pixelsregion(:,1),pixelsregion(:,2));

lats(end+1) = lats(1);
lons(end+1) = lons(1);

clf;
plot(lons,lats,'r.','MarkerSize',20);hold on;
plot(lons,lats)
[X,Y] = meshgrid([-75:-0.25:-95],[30:0.25:45]);

I = inpolygon(X,Y,lons,lats);

X = X(:);
Y = Y(:);
I = I(:);

plot(X,Y,'.');
plot(X(I),Y(I),'g.');