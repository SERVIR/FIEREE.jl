using PyCall, JSON

ee = pyimport("ee")
ee.Initialize()

kh = ee.FeatureCollection("USDOS/LSIB_SIMPLE/2017").filter(ee.Filter.eq("country_na","Cambodia"));
basins = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7").filterBounds(kh)

output_dir = "../testing"

files = filter(x->endswith(x, ".json"), readdir(output_dir))

features = Array{PyObject}(undef,length(files))

for i in 1:length(files)
    file = files[i]

    code = parse(Int,split(split(file,"_")[end],".")[1])

    f = open(file,"r")
    data_dict = JSON.parse(f)
    close(f)

    features[i] = ee.Feature(basins.filter(ee.Filter.eq("HYBAS_ID",code)).first()).set(data_dict)

end 

fc = ee.FeatureCollection(features)

task = ee.batch.Export.table.toAsset(collection=fc, description="fier_testing_results", assetId="users/kelmarkert/fier_testing_resutls")
task.start()
