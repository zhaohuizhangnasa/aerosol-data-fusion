from pydantic import BaseModel, HttpUrl

from typing import Sequence

#user config request
class Config(BaseModel):
    #gridsettings
    gridsize: float
    limit: list
    fill_value: float
    time_interval: int
    start_date: float
    start_time: float
    end_date: float
    end_time: float
    #variables
    geo_var: list
    phy_var: list
    phy_var_nc: list
    phy_var_hdf: list
    pixel_range: list
