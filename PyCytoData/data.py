import fcsparser
from pandas import DataFrame
import numpy as np
import cytomulate 

from zipfile import ZipFile
from urllib.request import urlopen
from io import BytesIO

import os
import pkg_resources
import glob


from typing import Optional, List, Dict, Literal, Any


class DataLoader():

    # Package data directory and Path
    _data_dir = pkg_resources.resource_filename("PyCytoData", "data/")
    _data_path: Dict[str, str] = {"levine13": _data_dir + "levine13/"}
    # Check data status
    _data_status: Dict[str, bool] = {}
    for d in _data_path.keys():
        _data_status[d] = os.path.exists(_data_path[d])

    @classmethod    
    def load_dataset(cls, dataset: Literal["levine13"], force_download: bool = False) -> List[np.ndarray]:
        
        if not cls._data_status[dataset]:
            cls._download_data(dataset = dataset, force_download = force_download)
            
        return cytomulate.utilities.FileIO.load_data(cls._data_path[dataset]+dataset+".txt", col_names = True)
    
    
    @classmethod
    def load_cell_types(cls, dataset: Literal["levine13"], force_download: bool = False) -> np.ndarray:
        
        if not cls._data_status[dataset]:
            cls._download_data(dataset = dataset, force_download = force_download)
            
        cell_type_path: str = cls._data_path[dataset] + dataset + "_cell_types.txt"
        return cytomulate.utilities.FileIO.load_data(cell_type_path, col_names = True)[1]
            

    @classmethod
    def _download_data(cls,
                      dataset: Literal["levine13"],
                      force_download: Optional[bool]=False) -> int:
        
        """Method to download datasets.

        """
        urls: Dict[str, str] = {"levine13": "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Levine_13dim/Levine_13dim_fcs_files.zip"}

        if not force_download:
            value = input(f"Would you like to download {dataset}? [y/n]")
            
            if value.lower() != "y":
                message_1 = f"\nYou have declined to download {dataset}.\n"
                print(message_1)

                return 1

        # Download message
        message_2 = "\nDownload in progress...\n"
        message_2 += "This may take quite a while, "
        message_2 += "go grab a coffee or cytomulate it!\n"
        print(message_2)
        
        contents = urlopen(urls[dataset])
        contents = contents.read()
        zip_file = ZipFile(BytesIO(contents))
        zip_file.extractall(cls._data_path[dataset])
        
        exprs: List[np.ndarray] = cls._preprocess(dataset)
        path: str = cls._data_path[dataset] + dataset + ".txt"
        labels_path: str = cls._data_path[dataset] + dataset + "_cell_types.txt"
        cytomulate.utilities.FileIO.save_np_array(exprs[2], path, col_names = exprs[0])
        cytomulate.utilities.FileIO.save_np_array(exprs[1], labels_path)
        return 0
    
    
    @classmethod
    def _preprocess(cls, dataset: Literal["levine13"]) -> List[np.ndarray]:
        fcss: List[str] = glob.glob(cls._data_path[dataset] + "*.fcs")
        fcs_cell_types: List[str] = [types.split("_")[1] for types in fcss]

        cell_types: np.ndarray = np.array([])
        exprs: DataFrame = DataFrame()
        meta: Dict[str, Any] = {}
        
        i: int
        p: str
        for i, p in enumerate(fcss):
            temp: DataFrame
            meta, temp = fcsparser.parse(p, reformat_meta=True)
            exprs = exprs.append(temp)
            cell_types = np.concatenate((cell_types, np.repeat(fcs_cell_types[i], temp.shape[0])))
            
        exprs_array: np.ndarray = exprs.to_numpy()
        colnames: np.ndarray = meta['_channels_']["$PnN"].to_numpy()
        
        return [colnames, cell_types, exprs_array]