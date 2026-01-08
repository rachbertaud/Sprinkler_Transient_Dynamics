fig = openfig('forward_800.fig'); %open figure

%%
dataObjs = findobj(fig, 'Type', 'line'); %find the line

xData = get(dataObjs, 'XData'); %get xData from the line
yData = get(dataObjs, 'YData'); %get yData from the line

close all %close the figure

%pull out peaks data
x_peaks = xData{1}; 
y_peaks = yData{1};

%pull out all data
full_X = xData{2};
full_y = yData{2};