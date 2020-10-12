# fier_ee

Repository to replicate Forecasting of Inundation Extent using REOF with Earth Engine

Workflow is an adaptation from https://doi.org/10.1016/j.rse.2020.111732

## Example use

```julia
# script for running the end-to-end fier process
include("src/FIEREE.jl")
using PyCall

ee = pyimport("ee");
ee.Initialize();

cloud_project = "<cloud-project-name>";
session = get_session("<path/to/ee/secret_key>");

# box = [104.45,15.1,105.1,15.35];
basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7");
test_basin = basins.filter(ee.Filter.eq("HYBAS_ID",4071080170));

box = geom_to_bbox(test_basin)
resolution = 0.005; 

start_time = "2017-01-01";
end_time = "2020-01-01";

domain = get_domain(box,resolution);

data_arr = get_s1_data(cloud_project,session,domain,start_time,end_time);
dates = get_s1_dates(domain,start_time,end_time);

spatial_modes,temporal_modes = reof(data_arr);

df = get_spt_table(geom=test_basin);
df_matched = match_dates(df,dates);

results = find_fits(df_matched,spatial_modes,temporal_modes,data_arr)


```
