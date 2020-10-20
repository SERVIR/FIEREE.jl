using FIEREE, PyCall, JSON, Dates

cloud_project = "ee-sandbox"
secret_key = "/home/Socrates/kmarkert/fier_testing/ee-sandbox.json"

# initialize Earth Engine
ee = pyimport("ee")
ee.Initialize()

session = FIEREE.get_session(secret_key)

# define parameters for space and time info
resolution = 0.005
start_time = "2015-01-01";
end_time = "2020-01-01";


kh = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017").filter(ee.Filter.eq("country_na","Cambodia"));

basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7").filterBounds(kh);
basin_codes = basins.aggregate_array("HYBAS_ID").getInfo()

cnt = 1
total = length(basin_codes)

for code in basin_codes
    now = Dates.now()
    println("$now -- ($cnt/$total): Processing basin code $code")
    basin = basins.filter(ee.Filter.eq("HYBAS_ID",code));

    box = FIEREE.geom_to_bbox(basin)

    domain = FIEREE.get_domain(box,resolution);
    println(domain)

    dates = FIEREE.get_s1_dates(domain,start_time,end_time);
    data_arr = FIEREE.get_s1_data(cloud_project,session,domain,start_time,end_time);

    spatial_modes,temporal_modes = FIEREE.reof(data_arr;n_modes=4);

    df = FIEREE.get_spt_table(basin);
    df_matched = FIEREE.match_dates(df,dates);

    results = FIEREE.find_fits(df_matched,spatial_modes,temporal_modes,data_arr);

    open("fier_results_$code.json","w") do io
        JSON.print(io,results,4)
    end

end