%% App with Auto-Reflow that Updates Plot Based on User Selections
% This app shows how to define controls and tabs within the panels of an 
% app with auto-reflow. The controls are in an anchored panel on the left. 
% The right panel that reflows contains two tabs. One tab displays a chart 
% and user interface components for adjusting the chart. The other tab 
% contains a table with the data used to make the chart. User selections 
% update both the plot and the table. The app responds to resizing by
% automatically growing, shrinking, and reflowing the app content.
% 
% The app includes these components:
% 
% * Check boxes &mdash; used to update the plot and table when the user 
% selects or clears a check box.
% * Switch &mdash; used to toggle the data that is visualized in the chart
% * Button group containing radio buttons &mdash; used to manage exclusive 
% selection of radio buttons. When the user selects a radio button, the button 
% group executes a callback function to update the plot with the appropriate data.
% * Slider &mdash; used to adjust histogram bin width. This slider only appears when the 
%  *Histogram* plotting option is selected in the button group.
% * Table &mdash; used to view the data associated with the chart.
% 
% <<../patientsdisplay_screenshot_19a.png>>  

% Copyright 2015 The MathWorks, Inc.