% Upslope Area Toolbox
% Version 1.4.2 17-Sep-2009
%
% Requires Image Processing Toolbox(TM).
%
% border_nans          - Find NaNs connected to DEM border
% dem_flow             - Downslope flow direction for a DEM
% dependence_map       - Dependence map for pixel flow in a DEM
% facet_flow           - Facet flow direction
% fill_sinks           - Fill interior sinks in a DEM
% flow_matrix          - Linear equations representing pixel flow
% influence_map        - Influence map for pixel flow in a DEM
% pixel_flow           - Downslope flow direction for DEM pixels
% postprocess_plateaus - Replace upslope areas for plateaus with mean value 
% upslope_area         - Upslope area measurements for a DEM
% vis_dem_flow         - Visualize flow directions in a DEM
% vis_map              - Visualize influence or dependence map for a DEM
%
% milford_ma_dem.mat   - Sample DEM data provided by USGS and distributed
%                        via Geo Community (geoworld.com), a USGS data
%                        distribution partner.  The data set is a 1:24,000-scale
%                        raster profile digital elevation model.  Download the
%                        "Milford" file from the "Digital Elevation Models (DEM)
%                        - 24K Middlesex County, Massachusetts, United States"
%                        page: 
%
%                        http://data.geocomm.com/catalog/US/61059/526/group4-3.html
%                        
%
% Steven L. Eddins
% Copyright 2007-2009 The MathWorks, Inc.
