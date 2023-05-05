from fastapi import FastAPI, APIRouter, Form, Query, HTTPException, Request, File, UploadFile
from fastapi.templating import Jinja2Templates
from fastapi.staticfiles import StaticFiles
from fastapi.responses import HTMLResponse, FileResponse, StreamingResponse
from pydantic import BaseModel
from typing import Optional, Any, List
from pathlib import Path

import os
from os import listdir
from os.path import isfile, join
import netCDF4
import shutil
import zipfile
from io import StringIO, BytesIO

# custom
from schemas import *
#from user_settings import configuration, uploaded_files
#from satellite_getter import *

from pyroscope.time_conv import *
from pyroscope.grid_ncf import *
from pyroscope.gridding import *

#grididng

BASE_PATH = Path(__file__).resolve().parent
TEMPLATES = Jinja2Templates(directory=str(BASE_PATH))


app = FastAPI() #(title="Recipe API", openapi_url="/openapi.json")
api_router = APIRouter()

#pathing
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
UPLOAD_DIR = os.path.join(BASE_DIR, "uploads")
STATIC_DIR = os.path.join(BASE_DIR, "static")

#saved data
configuration = {
    "gridsize":0.25,
    "limit": [-89.875, 89.875, -179.875, 179.875] ,
    "fill_value": -9999., #default is NAN
    "time_interval": 30, #minute intervals that will be split
    "time_start": "2020/01/01/00/00", # year/month/day/hour/min
    "time_end": "2020/01/01/00/59",   # year/month/day/hour/min
    "geo_var": ["latitude", "longitude"] #default (what geophysical variables are mapped to)
    ,"phy_var": ["Sensor_Zenith", "Scattering_Angle", "Image_Optical_Depth_Land_And_Ocean", "Optical_Depth_Land_And_Ocean"] # output names and master list
    ,"phy_var_nc": ["sensor_zenith_angle", "Scattering_Angle", "Image_Optical_Depth_Land_And_Ocean", "Optical_Depth_Land_And_Ocean"] # geophysical variables netCDF
    ,"phy_var_hdf": ["Sensor_Zenith", "Scattering_Angle", "Image_Optical_Depth_Land_And_Ocean","Optical_Depth_Land_And_Ocean"] # geophysical variables hdf
    ,"pixel_range": [0, 500] # Range for pixel count at single pixel
}
filenames = []
#filenames = ['AERDT_L2_ABI_G16.A2020001.0000.001.2022003073056.nc', 'AERDT_L2_ABI_G16.A2020001.0040.001.2022003074938.nc']
uploaded_files = []

@app.get("/")
async def index(request: Request):
    return TEMPLATES.TemplateResponse("index.html", {"request": request})

# get user request information
@app.post("/config")
async def submit_config(request: Request,
                 gridsize: float = Form(...),
                 limit: str = Form(...),
                 fill_value: float = Form(...),
                 time_interval: float = Form(...),
                 start_date: str = Form(...),
                 start_time: str = Form(...),
                 end_date: str = Form(...),
                 end_time: str = Form(...),
                 
                 geo_var: str = Form(...),
                 phy_var: str = Form(...),
                 phy_var_nc: str = Form(...),
                 phy_var_hdf: str = Form(...),
                 pixel_range: str = Form(...)
                 ):
    
    #time string edit
    time_start = start_date + "/" + start_time.replace(":", "/")
    time_end = end_date + "/" + end_time.replace(":", "/")
    
    #save into configuration
    configuration["gridsize"] = float(gridsize)
    configuration["limit"] = [float(num) for num in limit.strip('][').split(', ')]
    configuration["fill_value"] = float(fill_value)
    configuration["time_start"] = time_start
    configuration["time_end"] = time_end
    configuration["time_interval"] = int(time_interval)
    configuration["geo_var"] = geo_var.strip('][').split(', ')
    configuration["phy_var"] = phy_var.strip('][').split(', ')
    configuration["phy_var_nc"] = phy_var_nc.strip('][').split(', ')
    configuration["phy_var_hdf"] = phy_var_hdf.strip('][').split(', ')
    configuration["pixel_range"] = [int(num) for num in pixel_range.strip('][').split(', ')]
    
    user_config_html = """
    <html>
        <head>
            <link href="https://unpkg.com/tailwindcss@^2/dist/tailwind.min.css" rel="stylesheet">
        </head>
        <body>
            <h2 class="mb-8 text-2xl font-medium text-white">
                User Configurations:
            </h2>
            <br>
            <p id = "user_config" class="text-base text-white">"""
    user_config_html = user_config_html + str(configuration)
    user_config_html = user_config_html + """
            </p>
            <br><br>
        </body>
    </html>
    """
    return HTMLResponse(content=user_config_html, status_code=200)

#file upload
@app.post("/file")
async def submit_file(request: Request,
                      file: List[UploadFile]):
    
    if file[0].filename != "":
        filenames.extend([f.filename for f in file])
        uploaded_files.extend(file)
    else:
        print("No files")
    
    file_upload_html =  """
    <html>
        <head>
            <link href="https://unpkg.com/tailwindcss@^2/dist/tailwind.min.css" rel="stylesheet">
        </head>
        <body>
            <h2 class="mb-8 text-2xl font-medium text-white">
                Files Uploaded:
            </h2>
            <br>
            <p id="file_uploads" class="text-base text-white">"""
    file_upload_html = file_upload_html + str(filenames)
    file_upload_html = file_upload_html + """
            </p>
            <br><br>
        </body>
    </html>
    """
    
    # save file to uploads folder
    for i in range(len(filenames)):
        file_location = os.path.join(UPLOAD_DIR, filenames[i])
        with open(file_location, "wb+") as file_object:
            shutil.copyfileobj(uploaded_files[i].file, file_object)   
            
    print("files uploaded")
    
    return HTMLResponse(content=file_upload_html, status_code=200)

# generate response
@app.post("/generate")
async def generate_file(request: Request):
    print("Generating for ... ", filenames)
    
    # configuration read in 
    filelist = [str(os.path.join(UPLOAD_DIR, filename)) for filename in filenames]
    gsize = configuration["gridsize"]
    geo_list = configuration["geo_var"]
    phy_list = configuration["phy_var"]
    phy_nc = configuration["phy_var_nc"]
    phy_hdf = configuration["phy_var_hdf"]
    limit = configuration["limit"]
    pixel_range = configuration["pixel_range"]
    static_file = str(os.path.join(STATIC_DIR, "LSM_ELV_QDEG_FIXED.nc"))
    start = to_datetime(configuration["time_start"])
    end = to_datetime(configuration["time_end"])
    time_interval = configuration["time_interval"]
    
    #time processing    
    #print(filelist)
    split_files = split_filetimes(filelist, start, end, int(time_interval))
    split_files = split_filenames(split_files)
    
    curr = start
    
    #begin processing time period by time period
    #as bucketed in split_files
    for f in split_files:
        #print("FILELIST:",f)
        
        #XAERDT_L3_MEASURES_QD_HH.YYYYMMDD.HHMM.V0.ProcessingDate.nc
        time_string = str(curr).replace(":", "").replace("-", "").replace(" ", "")
        time_string = time_string[:8] + '.' + time_string[8:-2]
        processing_date = sys_time()
        filename = "XAERDT_L3_MEASURES_QD_HH." + time_string+ ".V0." +processing_date+ ".nc"
        filelocation = str(os.path.join(STATIC_DIR, filename))
        print(filelocation)
        
        ds = grid_nc_sensor_statistics_metadata(limit, gsize, geo_list, phy_list, f, filelocation, curr, time_interval, phy_nc, phy_hdf, static_file, pixel_range)
        ds.close()
        
        curr = curr + datetime.timedelta(minutes = int(time_interval))
    
    """
    for filename in filenames:
        
        file_location = os.path.join(UPLOAD_DIR, filename)
        
        L2FID = netCDF4.Dataset(file_location,'r',format='NETCDF4')
        print(L2FID)
        L2FID.close()
        
        print("opened and close")
    """
    generate_file_html =  """
    <html>
        <head>
            <link href="https://unpkg.com/tailwindcss@^2/dist/tailwind.min.css" rel="stylesheet">
        </head>
        <body>
            <h2 class="mb-8 text-2xl font-medium text-white">
                Generating files
            </h2>
        </body>
    </html>
    """
    
    return HTMLResponse(content=generate_file_html, status_code=200)

# download file response
@app.post("/download")
async def download_file(request: Request):
    print("Downloading zip file... ")
    
    #f = open(str(BASE_PATH)+"/testing.nc", "x")
    
    #SAVE_FILE_PATH = os.path.join(BASE_PATH, "testing.nc")
    #print("created")
    
    # Return as a download
    #headers = {'Content-Disposition': 'attachment; filename="testing.nc"'}
    
    #return FileResponse(
    #    path=SAVE_FILE_PATH,
    #    media_type="application/netcdf", #	application/netcdf #text/plain
    #    headers=headers
    #)
    
    zip_subdir = "archive"
    zip_io = BytesIO()
    # The zip compressor
    #zf = zipfile.ZipFile(zip_io, "w")
    with zipfile.ZipFile(zip_io, mode='w', compression=zipfile.ZIP_DEFLATED) as temp_zip:
        for file_path in [f for f in listdir(STATIC_DIR) if isfile(join(STATIC_DIR, f))]:
            # Calculate path for file in zip
            fdir, fname = os.path.split(file_path)
            file_path = os.path.join(STATIC_DIR, fname)
            zip_path = os.path.join(zip_subdir, fname)
            
            #print(file_path)
            #print(zip_path)
            
            # Add file, at correct path
            temp_zip.write(file_path, zip_path)
    
    # Must close zip for all contents to be written
    #zf.close()
    
    return StreamingResponse(
        iter([zip_io.getvalue()]), 
        media_type="application/x-zip-compressed", 
        headers = { 'Content-Disposition': f'attachment; filename="datafusion.zip"'} 
    ) 

app.include_router(api_router)

if __name__ == "__main__":
    # Use this for debugging purposes only
    import uvicorn

    uvicorn.run(app, host="0.0.0.0", port=8001, log_level="debug")
