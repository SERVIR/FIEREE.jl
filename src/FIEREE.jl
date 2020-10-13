module FIEREE

using PyCall, JSON, Statistics, LinearAlgebra, DataFrames, Dates, HTTP, ProgressMeter, NCDatasets, Polynomials

ee = pyimport("ee")

struct Domain
    bbox::Vector
    resolution::Float64
    shape::Vector
    crs::String
end
    
function round_out(bbox::Vector,resolution::Float64)
    min_buff = bbox .- (bbox .% resolution)
    min_buff = map(x->round.(x,digits=10), min_buff)

    max_buff = bbox .+ (resolution .- (bbox .% resolution))
    max_buff = map(x->round.(x,digits=10), max_buff)

    new_bbox = vcat(min_buff[1:2],max_buff[3:end])
end

function get_domain(bbox::Vector,resolution::Float64; crs::String="EPSG:4326")
    if any((bbox .% resolution) .> 0 )
        bbox = round_out(bbox,resolution)
    end

    xdim = convert(Int,round((bbox[3]-bbox[1]) / resolution))
    ydim = convert(Int,round((bbox[4]-bbox[2]) / resolution))
    shape = [xdim,ydim]

    Domain(bbox,resolution,shape,crs)
end

function domain_to_coords(domain::Domain)
    res = domain.resolution
    west,east = domain.bbox[1], (domain.bbox[3] - res)
    south,north = (domain.bbox[2] - res), domain.bbox[4]
    xx = collect(Array{Float32}(west:res:east))
    yy = collect(Array{Float32}(south:res:north))
    return xx,yy
end

function get_ee_region(bbox::Vector)
    minx, miny, maxx, maxy = bbox
    coordinates = [[minx,miny],[minx,maxy],[maxx,maxy],[maxx,miny],[minx,miny]]

    ee.Geometry.Polygon(coordinates)
end

function geom_to_bbox(geom)
    coords = geom.geometry().bounds().coordinates().getInfo()[1,:,:]
    minx, maxx = min(coords[:,1]...), max(coords[:,1]...)
    miny, maxy = min(coords[:,2]...), max(coords[:,2]...)

    return [minx,miny,maxx,maxy]

end

function response_to_array(response,width,height;channels=1,dtype=Float64)
    buffer = IOBuffer()
    write(buffer,response.content)
    foo = reinterpret(dtype,take!(buffer))
    close(buffer) # close the io buffer for clean memory management

    # find the number of pixels we expect based on img 
    data_size::Int64 = width*height*channels
    # find the offset of what was returned vs what is expected
    offset::Int64 = size(foo,1)-data_size + 1 

    arr::Array{dtype} = dtype.(foo)[offset:end]
    
    if channels > 1
        # reshape array from 1-d to correct 3-d dimensions
        arr3d::Array{dtype} = reshape(arr,(channels,height,width))
        # reorder the dimensions for row-major to column-major
        img = permutedims(arr3d,(1,3,2))
    else
        arr2d::Array{dtype} = reshape(arr,(width,height))
        img = rotl90(arr2d)
    end

    # cast the result as a float32 array 
    img = convert(Array{Float32}, img)

end 

function exponential_backoff(c::Int)
    N = 2^c -1
    backoffTime = (1/(N+1)) * sum(Array(1:1:N))
    bt = max(backoffTime,0.1) # force backoffs to be at least 0.1s
end

function get_s1_dates(domain::Domain,start_time::String,end_time::String)
    geom = get_ee_region(domain.bbox)
    s1 = (ee.ImageCollection("COPERNICUS/S1_GRD").
        filterDate(start_time,end_time).
        filterBounds(geom).
        filter(ee.Filter.eq("orbitProperties_pass","ASCENDING")).
        select("VV").
        sort("system:time_start")
    )
    
    date_str = ee.List(s1.map(img_to_date).aggregate_array("isodate")).distinct().getInfo()
    # convert epoch milliseconds to date
    # Julia considers start of epoch Jan 1, 0000 so we add 1970 years to convert to correct times
    dates = map(x->DateTime(x),date_str)

end

# helper function to set date info, used for getting list of dates
function img_to_date(img)
    return img.set("isodate",img.date().format("YYYY-MM-dd"))
end

# helper function to mask high view angle areas
function s1_qa(img)
    angle = img.select("angle")
    return img.updateMask(angle.gt(32).And(angle.lt(45)))
end

function get_s1_data(project::String,session,domain::Domain,start_time::String,end_time::String)
    compute_image_url = "https://earthengine.googleapis.com/v1alpha/projects/$project/image:computePixels"
    res = domain.resolution
    xorigin = domain.bbox[1]
    yorigin = domain.bbox[4]
    w,h = domain.shape
    crs = domain.crs
    revisit = 6 # days

    max_retries = 5
    retry_offset = 2 # used to control where to start exponential_backoff
    n_retries = max_retries+retry_offset

    geom = get_ee_region(domain.bbox)

    s1 = (ee.ImageCollection("COPERNICUS/S1_GRD").
        filterDate(start_time,end_time).
        filterBounds(geom).
        filter(ee.Filter.eq("orbitProperties_pass","ASCENDING")).
        sort("system:time_start").
        map(s1_qa).
        select("VV")
    )
    
    dates = ee.List(s1.map(img_to_date).aggregate_array("isodate")).distinct().getInfo()
    n = length(dates)

    jrc = ee.Image("JRC/GSW1_2/GlobalSurfaceWater").select("occurrence").unmask(0).lt(80)

    # n::Int64 = s1.size().getInfo()
    # img_list = s1.toList(n)

    result_arr::Array{Float32} = zeros(Float32,n,h,w,1)
    # result_arr =fill!(Array{Union{Missing, Float64}}(missing,n,h,w,c),-9999)

    @showprogress for i in 1:n
        eeDate = ee.Date(dates[i])
        img = ee.Image(s1.filterDate(eeDate,eeDate.advance(revisit,"day")).mean()).
            updateMask(jrc)

        serialized = ee.serializer.encode(img, for_cloud_api=true)
        
        payload = Dict(
            "expression" => serialized, 
            "fileFormat" => "NPY", 
            "grid" => Dict(
                "affineTransform" => Dict(
                    "scaleX" => res,
                    "scaleY" => -res,
                    "translateX" => xorigin,
                    "translateY" => yorigin
                ),
                "dimensions" => Dict(
                    "width" => w,
                    "height" => h
                ),
                "crsCode"=> crs
            )
        )

        for j in retry_offset:n_retries
            response = session.post(
                url = compute_image_url,
                data = JSON.json(payload)
            )

            if response.status_code == 200
                arr = response_to_array(response,w,h)
                result_arr[i,:,:,:] = reshape(arr,(h,w,1))
                break
            else
                sleep(exponential_backoff(j))
            end

            if j == n_retries
                println(response.text)
                throw(Base.ProcessFailedException)
                return
            end
        end

    end

    # find no data values and set to missing
    #result_arr[result_arr.==-9999.] .= missing

    return result_arr  

end


function get_session(secret_key)
    # import gcloud authentication packages for interacting with the EE REST API 
    service_account = pyimport("google.oauth2.service_account")
    gauth = pyimport("google.auth.transport.requests")
    # define the credentials to authorize
    credentials = service_account.Credentials.from_service_account_file(secret_key)
    scoped_credentials = credentials.with_scopes(
        ["https://www.googleapis.com/auth/cloud-platform"])

    # authorize a session to start sending requests
    session = gauth.AuthorizedSession(scoped_credentials)
end

function reof(array::Array{Float32}; variance_threshold::Float64=0.727,n_modes::Int=-1,max_rotations=50)
    spatial_shape = size(array)[2:3]
    flattened = reshape(array,(size(array)[1],prod(spatial_shape)))
    time_shape = size(flattened)[1]

    spatial_mean = mean(flattened,dims=1)
    # time_mean = mean(flattened,dims=2)
    centered = flattened .- spatial_mean

    cov_mat = cov(centered)

    eigen_values,eigen_vectors = eigen(cov_mat)

    # reverse the order of the eigen values and vectors so variance explained is in decreasing order
    eigen_values = eigen_values[end:-1:1]
    eigen_vectors = eigen_vectors[:,end:-1:1]
    
    if n_modes <= 0
        δλ = eigen_values .* sqrt(2 / time_shape)
        explained_var = cumsum(δλ) / sum(δλ)
        n_modes = sum(explained_var .< variance_threshold)
    end

    eofs = eigen_vectors[:,1:n_modes]

    pcs = Array{Float64}(undef,(time_shape,n_modes))
    for i in 1:n_modes
        pcs[:,i] = centered * eigen_vectors[:,i]
    end

    spatial_modes, temporal_modes = varimax_rotation(eofs,pcs)

    spatial_modes = reshape(spatial_modes,(vcat(spatial_shape...,n_modes)...))

    return spatial_modes, temporal_modes

end

function varimax_rotation(U,A; max_iter=50)
    Z = A * U' # reconstructed data

    #begin varimax rotations
    B = Z * U

    lim = 1
    tot = 1e-6
    cnt = 1

    # need to just initialize V and Ginv
    modes = size(U)[2]
    V = zeros(size(U))
    Ginv = zeros((modes,modes))

    for i in 1:max_iter
    # while lim > tot
        D = diagm(mean(B.^2,dims=1)[1,:])
        C = B.^3
        V = Z' * (C - B * D)
        k,w = eigen((V'*V))
        k = k.^(-1/2)
        for j in 1:size(k)[1]
            if !isfinite(k[j])
                k[j] = 0
            end
        end

        kin = diagm(k)

        Ginv = w * kin * w'

        Uold = U
        U = V * Ginv
        
        Bold = B
        B = Z*U

        if max((U - Uold)...) < tot
            break
            # return U,B
        end
        cnt += 1
    end

    if cnt >= max_iter
        println("rotation never converged...")
    else
        println("rotation coverged after $cnt iterations")
    end

    U = V * Ginv
    B = Z*U

    return U, B

end

function dualvarimax(Z,Phi)
    U = Phi
    B = Z' * U
    # Unew = zeros(size(Phi))

    lim = 1
    tot = 1e-6
    max_iter = 50
    cnt = 1

    # need to just initialize V and Ginv
    modes = size(U)[2]
    V = zeros(size(U))
    Ginv = zeros((modes,modes))

    for i in 1:max_iter
    # while lim > tot
        D = diagm(mean(B.^2,dims=1)[1,:])
        C = B.^3
        V = Z * (C - B * D)
        k,w = eigen((V'*V))

        kin = diagm(k.^(-1/2))

        Ginv = w * kin * w'

        Uold = U
        U = V * Ginv
        
        Bold = B
        B = Z'*U

        if max((U - Uold)...) < tot
            break
            # return U,B
        end
        cnt += 1
    end

    U = V* Ginv
    B = Z' * U

    return B, U

end

function get_spt_table(geom=nothing,domain=nothing, use_reach_id=false)
    # define the api service url and endpoint
    spt_api_url = "https://tethys2.byu.edu/localsptapi/api"
    endpoint = "HistoricSimulation"

    if geom === nothing
        geom = get_ee_region(domain.bbox)
    end

    # geom = get_ee_region(domain.bbox)
    reaches = ee.FeatureCollection("users/kelmarkert/public/geoglows/global-geoglows-drainageline")
    roi_reaches = reaches.filterBounds(geom).sort("to_node",false)
    # need to update filtering process
    # currently gets first reach but needs to be outlet reach
    if use_reach_id
        reach_id = ee.Feature(roi_reaches.first()).get("arcid").getInfo()
        request_url = "$spt_api_url/$endpoint/?reach_id=$reach_id&return_format=json"
    else
        x,y = ee.List(ee.Feature(roi_reaches.first()).geometry().coordinates()).get(0).getInfo()
        request_url = "$spt_api_url/$endpoint/?lon=$x&lat=$y&return_format=json"
    end

    response = HTTP.get(request_url)
    if response.status == 200
        data_dict = JSON.parse(String(response.body))

        df = DataFrame(data_dict["time_series"])
        df[:,"datetime"] = DateTime.(SubString.(df[:,"datetime"],1,19),DateFormat("Y-m-dTH:M:Si"))
        df[:,"flow"] = convert(Array{Float64},df.flow)
    else
        throw(Base.ProcessFailedException)
        return
    end
    return df
end

function match_dates(df::DataFrame,dates::Vector{DateTime})
    # TODO: check to make sure that dates range is within df[:,"datetime"] range
    # produces inconsistent results if not ... not sure why...
    mask = sum(map(x -> x.==df[:,"datetime"],dates),dims=1)[1]
    indices = findall(!iszero,mask)

    df[indices,:]

end

function find_fits(df::DataFrame,sm::Array{Float64},tcp::Array{Float64},stack::Array{Float32})
    # TODO: best fitting???
    # simplisting fitting and will return correlation of tcp to streamflow with rmse of predictions
    x = convert(Array{Float64},df.flow)

    output = Dict{String,Float64}()

    spatial_mean = mean(stack[:,:,:],dims=1)[1,:,:]
    flat_shape = prod(size(spatial_mean))

    errors = Array{Any}(undef,size(df)[1])
    cors = Array{Any}(undef,size(df)[1])

    for i in 1:size(tcp)[2]
        mode_name = "tcp$i"
        y = tcp[:,i]

        for j in 1:3
            model = fit(x,y,j)

            test_preds = syntesize(model,df.flow,sm[:,:,i],spatial_mean)
            
            for k in 1:size(test_preds)[1]
                diff = test_preds[k,:,:] - stack[k,:,:,1]
                errors = mean(abs.(diff))
                cors = cor(reshape(test_preds[k,:,:],flat_shape),reshape(stack[k,:,:,1],flat_shape))
            end

            push!(output, "$(mode_name)_fit_correlation_order$j" => cor(y,model.(x)))
            push!(output, "$(mode_name)_pred_correlation_order$j" => mean(cors))
            push!(output, "$(mode_name)_pred_rmse_order$j"=>mean(errors))
        end

    end

    return output
    
end

function syntesize(model,streamflow,spatial_mode,spatial_mean)
    z = model.(streamflow)

    preds = Array{Float64}(undef,(vcat(length(z),size(spatial_mode)...)...))

    for i in 1:length(z)
        preds[i,:,:] = z[i]*spatial_mode+spatial_mean
    end

    return preds
end

function persist_data(fname::String, dataset, domain, dates)
    # TODO: write data persistance function to save either sar or synthesized data

end

# end module
end
