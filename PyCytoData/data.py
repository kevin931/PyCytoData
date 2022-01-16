import pandas as pd
import numpy as np
from PyCytoData import exceptions

from zipfile import ZipFile
from urllib.request import urlopen
from io import BytesIO

import fcsparser
import _csv
import csv
import os
import pkg_resources
import glob

from numpy.typing import ArrayLike
from typing import Optional, List, Dict, Literal, Any, Union


class CytoData():
    
    def __init__(self,
                 expression_matrix: ArrayLike,
                 features: Optional[ArrayLike]=None,
                 cell_types: Optional[ArrayLike]=None,
                 sample_index: Optional[ArrayLike]=None):
        
        
        self._expression_matrix: np.ndarray = np.array(expression_matrix)
        if len(self._expression_matrix.shape) != 2:
            raise exceptions.ExpressionMatrixDimensionError(shape=self._expression_matrix.shape)
        
        self._n_cells = self._expression_matrix.shape[0]
        self._n_features = self._expression_matrix.shape[1]
        
        if features is not None:
            self._features = np.array(features)
        else:
            self._features = np.full(self.n_features, None)
        
        if cell_types is None:
            self._cell_types = np.full(self.n_cells, None)
        else:
            self._cell_types = np.array(cell_types)
            
        if sample_index is None:
            self._sample_index = np.repeat(0, self.n_cells)
        else:
            self._sample_index = np.array(sample_index)
        
        self._n_samples: int = len(set(self._sample_index))
        self._n_cell_types: int = len(set(self._cell_types))
        
        if self._features.shape[0] != self.n_features:
            raise exceptions.DimensionMismatchError(n=self.n_features, var = "features")
        if self._cell_types.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "cell_types")
        if self._sample_index.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "sample_index")
    
    
    def add_sample(self, expression_matrix: ArrayLike, sample_index: ArrayLike, cell_types: Optional[ArrayLike]=None):
        expression_matrix = np.array(expression_matrix)
        sample_index = np.array(sample_index)
        
        if len(expression_matrix.shape) != 2:
            raise exceptions.ExpressionMatrixDimensionError(expression_matrix.shape)
        if sample_index.shape[0] != expression_matrix.shape[0]:
            raise exceptions.DimensionMismatchError(n=expression_matrix.shape[0], var = "sample_index")
        if cell_types is not None and np.array(cell_types).shape[0] != expression_matrix.shape[0]:
            raise exceptions.DimensionMismatchError(n=expression_matrix.shape[0], var = "cell_types")
        
        self.expression_matrix = np.concatenate((self.expression_matrix, expression_matrix))
        self.sample_index = np.concatenate((self.sample_index, sample_index))
        
        if cell_types is None:
            self.cell_types = np.concatenate((self.cell_types, np.full(expression_matrix.shape[0], None)))
        else:
            self.cell_types = np.concatenate((self.cell_types, np.array(cell_types)))
         
         
    @property
    def expression_matrix(self):
        return self._expression_matrix   
    
    
    @expression_matrix.setter
    def expression_matrix(self, expression_matrix: ArrayLike):
        expression_matrix = np.array(expression_matrix)
        if len(expression_matrix.shape) != 2:
            raise exceptions.ExpressionMatrixDimensionError(expression_matrix.shape)
        self.n_cells = expression_matrix.shape[0]
        self.n_features = expression_matrix.shape[1]
        self._expression_matrix = expression_matrix
        
  
    @property
    def sample_index(self):
        return self._sample_index
    
    
    @sample_index.setter
    def sample_index(self, sample_index: ArrayLike):
        sample_index = np.array(sample_index)
        if sample_index.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "sample_index")
        self._sample_index = sample_index
        self.n_samples = len(set(self._sample_index))
        
        
    @property
    def cell_types(self):
        return self._cell_types
    
    
    @cell_types.setter
    def cell_types(self, cell_types: ArrayLike):
        cell_types = np.array(cell_types)
        if cell_types.shape[0] != self.n_cells:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "cell_types")
        self._cell_types = cell_types
        self.n_cell_types = len(set(self.cell_types))
        
        
    @property
    def features(self):
        return self._features
    
    
    @features.setter
    def features(self, features: ArrayLike):
        features = np.array(features)
        if features.shape[0] != self.n_features:
            raise exceptions.DimensionMismatchError(n=self.n_cells, var = "features")
        self._features = features
        
        
    @property
    def n_cells(self):
        return self._n_cells
    
    
    @n_cells.setter
    def n_cells(self, n_cells: int):
        if not isinstance(n_cells, int):
            raise TypeError(f"'n_cells' has to be 'int' instead of {type(n_cells)}")
        self._n_cells = n_cells
        

    @property
    def n_features(self):
        return self._n_features
    
    
    @n_features.setter
    def n_features(self, n_features: int):
        if not isinstance(n_features, int):
            raise TypeError(f"'n_features' has to be 'int' instead of {type(n_features)}")
        self._n_features = n_features
        
        
    @property
    def n_samples(self):
        return self._n_samples
    
    
    @n_samples.setter
    def n_samples(self, n_samples: int):
        if not isinstance(n_samples, int):
            raise TypeError(f"'n_samples' has to be 'int' instead of {type(n_samples)}")
        self._n_samples = n_samples
        
        
    @property
    def n_cell_types(self):
        return self._n_cell_types
    
    
    @n_cell_types.setter
    def n_cell_types(self, n_cell_types: int):
        if not isinstance(n_cell_types, int):
            raise TypeError(f"'n_samples' has to be 'int' instead of {type(n_cell_types)}")
        self._n_cell_types = n_cell_types
        


class DataLoader():

    # Package data directory and Path
    _data_dir = pkg_resources.resource_filename("PyCytoData", "data/")
    _data_path: Dict[str, str] = {"levine13": _data_dir + "levine13/",
                                  "levine32": _data_dir + "levine32/",
                                  "samusik": _data_dir + "samusik/"}
    # Check data status
    _data_status: Dict[str, bool] = {}
    for d in _data_path.keys():
        _data_status[d] = os.path.exists(_data_path[d])

    @classmethod    
    def load_dataset(cls, dataset: Literal["levine13", "levine32", "samusik"], force_download: bool = False) -> CytoData:
        
        if not cls._data_status[dataset]:
            cls._download_data(dataset = dataset, force_download = force_download)
            
        data: CytoData = FileIO.load_delim(cls._data_path[dataset]+dataset+".txt", col_names = True)
        
        cell_type_path: str = cls._data_path[dataset] + dataset + "_cell_types.txt"
        if os.path.exists(cell_type_path):
            cell_types: np.ndarray = np.loadtxt(fname=cell_type_path, dtype ="str", delimiter="\t")
            data.cell_types = cell_types
            
        return data
            

    @classmethod
    def _download_data(cls,
                      dataset: Literal["levine13", "levine32", "samusik"],
                      force_download: bool=False) -> int:
        
        """Method to download datasets.

        """
        urls: Dict[str, List[str]] = {"levine13": ["http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Levine_13dim/Levine_13dim_fcs_files.zip"],
                                       "levine32": ["http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Levine_32dim/Levine_32dim_fcs_files.zip"],
                                       "samusik": ["http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Samusik/Samusik_fcs_files.zip",
                                                   "http://imlspenticton.uzh.ch/robinson_lab/HDCytoData/Samusik/Samusik_population_IDs.zip"]}

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
        
        for url in urls[dataset]:
            contents = urlopen(url)
            contents = contents.read()
            zip_file = ZipFile(BytesIO(contents))
            print("Here")
            zip_file.extractall(cls._data_path[dataset])
        
        data: CytoData = cls._preprocess(dataset)
        path: str = cls._data_path[dataset] + dataset + ".txt"
        cell_types: np.ndarray = data.cell_types.reshape(data.n_cells, 1)
        sample_index: np.ndarray = data.sample_index.reshape(data.n_cells, 1)
        meta_data: np.ndarray = np.concatenate((cell_types, sample_index), axis=1)
        metadata_path: str = cls._data_path[dataset] + dataset + "_metadata.txt"
        FileIO.save_np_array(data.expression_matrix, path, col_names = data.features)
        FileIO.save_np_array(meta_data, metadata_path, dtype="%s")
        return 0
    
    
    @classmethod
    def _preprocess(cls, dataset: Literal["levine13", "levine32", "samusik"]) -> CytoData:
        
        fcss: List[str] = glob.glob(cls._data_path[dataset] + "*.fcs")
        exprs: pd.DataFrame = pd.DataFrame()
        meta: Dict[str, Any] = {}
        sample_length = []
        temp: pd.DataFrame
        fcs: str
        for fcs in fcss:
            meta, temp = fcsparser.parse(fcs, reformat_meta=True)
            sample_length.append(len(temp))
            exprs = exprs.append(temp)
        
        if dataset == "levine13":
            data: CytoData = cls._preprocess_levine13(fcss, exprs, meta, sample_length)
        elif dataset == "levine32":
            data: CytoData = cls._preprocess_levine32(fcss, exprs, meta, sample_length)
        elif dataset == "samusik":
            data: CytoData = cls._preprocess_samusik(fcss, exprs, meta, sample_length)
        
        return data
    
    
    @classmethod
    def _preprocess_levine13(cls, fcss: List[str], exprs: pd.DataFrame, meta: Dict[str, Any], sample_length: List[int]) -> CytoData:
        expression_matrix: np.ndarray = exprs.to_numpy() 
        # Cell Types
        fcs_cell_types: List[str] = [types.split("_")[-2] for types in fcss]
        cell_types: np.ndarray = np.array([])
        for s in range(len(sample_length)):
            cell_types = np.concatenate((cell_types, np.repeat(fcs_cell_types[s], sample_length[s])))
        # Features
        colnames: np.ndarray = meta['_channels_']["$PnN"].to_numpy()
        # Construct data
        data: CytoData = CytoData(expression_matrix, colnames, cell_types)
        
        return data
    
    
    @classmethod
    def _preprocess_levine32(cls, fcss: List[str], exprs: pd.DataFrame, meta: Dict[str, Any], sample_length: List[int]) -> CytoData:
        colnames: np.ndarray = np.array(exprs.columns)
        expression_matrix: np.ndarray = exprs.to_numpy()
        # Cell Types and Sample Index
        fcs_cell_types: List[str] = [types.split("_")[-2] for types in fcss]
        cell_types: np.ndarray = np.array([])
        sample_index: np.ndarray = np.array([])
        sample_names: List[str] = [types.split("-")[3] for types in fcss]
        for s in range(len(sample_length)):
            cell_types = np.concatenate((cell_types, np.repeat(fcs_cell_types[s], sample_length[s])))
            sample_index = np.concatenate((sample_index, np.repeat(sample_names[s], sample_length[s])))
        # Construct data
        data: CytoData = CytoData(expression_matrix, colnames, cell_types, sample_index)
        
        return data
    
    
    @classmethod
    def _preprocess_samusik(cls, fcss: List[str], exprs: pd.DataFrame, meta: Dict[str, Any], sample_length: List[int]) -> CytoData:
        colnames: np.ndarray = np.array(exprs.columns)
        expression_matrix: np.ndarray = exprs.to_numpy()
        # Sample Index
        sample_index: np.ndarray = np.array([])
        sample_names: List[str] = [types.split("_")[3] for types in fcss]
        for s in range(len(sample_length)):
            sample_index = np.concatenate((sample_index, np.repeat(sample_names[s], sample_length[s])))
            
        # Cell Types
        meta_data: pd.DataFrame = pd.read_csv(cls._data_path["samusik"] + "population_assignments.txt", delimiter = "\t", header = None) #type: ignore
        temp: pd.DataFrame = meta_data[0].str.split("_| ", expand = True) #type: ignore
        meta_data = pd.concat([meta_data, temp], axis = 1)
        meta_array: np.ndarray = meta_data.iloc[:,[1,5,8]].to_numpy()
        
        cell_types: List[str] = ["unassigned"]*expression_matrix.shape[0]
        unique_sample: np.ndarray = np.unique(meta_array[:,1])
        counter: int = 0
        for i, s in enumerate(sample_length):
            subset: np.ndarray = meta_array[meta_array[:,1] == unique_sample[i]]
            for j in range(len(subset)):
                try:
                    cell_types[counter + int(subset[j,2])] = subset[j, 0]
                except TypeError:
                    continue
            counter += s
        
        # Construct data
        data: CytoData = CytoData(expression_matrix, colnames, cell_types, sample_index)
        
        return data


class FileIO():
    
    @staticmethod
    def load_delim(files: Union[List[str], str],
                   col_names: bool=True,
                   drop_columns: Optional[Union[int, List[int]]]=None,
                   delim: str="\t",
                   dtype = float
                   ) -> CytoData:
        
        if not isinstance(files, str) and not isinstance(files, list):
            raise TypeError("'files' has to be str or a list of str as paths.")
        elif not isinstance(files, list):
            files = [files]
        
        skiprows: int = 1 if col_names else 0
        colnames: Optional[np.ndarray] = None
        data: Optional[CytoData] = None
        
        for i, file in enumerate(files):                            
            # Load Data
            f: "np.ndarray" = np.loadtxt(fname=file, dtype=dtype, skiprows=skiprows, delimiter=delim)
            sample_index: np.ndarray = np.full(f.shape[0], i)
            if drop_columns is not None:
                f = np.delete(f, drop_columns, axis=1)
            if i==0:
                data = CytoData(expression_matrix=f, sample_index=sample_index)
            else:
                assert data is not None
                data.add_sample(expression_matrix=f, sample_index=sample_index)
                
            # Load column names
            if i==0:
                if col_names:
                    colnames = np.loadtxt(fname=file, dtype ="str", max_rows=1, delimiter=delim)
                    if drop_columns is not None:
                        colnames = np.delete(colnames, drop_columns)
                    skiprows = 1
                else:
                    colnames = np.full(data.n_features, None)
        
        assert data is not None
        assert colnames is not None
        data.features = colnames
        
        return data
    
    
    @staticmethod
    def save_2d_list_to_csv(data: List[List[Any]], path: str, overwrite: bool = False):
        """Save a nested list to a CSV file.

        :param data: The nested list to be written to disk
        :type data: List[List[Any]]
        :param path: Path to save the CSV file
        :type path: str
        
        ..note:: 
        
            By default, this method does not overwrite existing files. In case a file exists,
            a ``FileExistsError`` is thrown.
        """
        if os.path.exists(path) and not overwrite:
            raise FileExistsError()
        
        i: int
        j: int
        
        with open(path, "w") as f:      
            w: "_csv._writer" = csv.writer(f)
            for i in range(len(data[0])):
                row: List[Any] = []
                for j in range(len(data)):
                    row.append(data[j][i])
                w.writerow(row)
            
            
    @staticmethod
    def save_np_array(array: "np.ndarray",
                      path: str,
                      col_names: Optional["np.ndarray"]=None,
                      dtype: str="%.18e",
                      overwrite: bool = False) -> None:
        """Save a NumPy array to a plain text file

        :param array: The NumPy array to be saved
        :type array: np.ndarray
        :param file: Path to save the plain text file
        :type file: str
        :param col_names: Column names to be save as the first row, defaults to None
        :type col_names: np.ndarray, optional
        :param dtype: NumPy data type, defaults to "%.18e"
        :type dtype: str, optional
        
        ..note:: By default, this method does not overwrite existing files. In case a file exists,
            a ``FileExistsError`` is thrown.
        """
        if os.path.exists(path) and not overwrite:
            raise FileExistsError()
            
        with open(path, "w") as f:
            if col_names is not None:
                f.write("\t".join(list(map(str, col_names))))
                f.write("\n")
            np.savetxt(f, array, delimiter="\t", fmt=dtype)
            
    
    @staticmethod
    def make_dir(dir_path: str, add_number_if_dir_exists: bool = False, _counter: int=0) -> str:
        """Create a new directory

        :param dir_path: Path to the new directory to be created
        :type dir_path: str
        :param add_number_if_dir_exists: If the directory already exists, append a number to the
            name until creation of directory is successful, defaults to True
        :type add_number_if_dir_exists: bool, optional
        :return: The path to the new directory
        :rtype: str
        
        .. Warning:: 
        
            The ``add_number_if_dir_exists`` can be dangerous when this method is run inside
            of a loop. This behavior may be deprecated and removed in the future.
        """
        dir_path = dir_path.rstrip("/")
        if _counter==0:
            new_dir_path = dir_path
        else:
            new_dir_path = dir_path + str(_counter)
        
        try:
            os.makedirs(new_dir_path)
        except FileExistsError:
            if add_number_if_dir_exists: 
                new_dir_path = FileIO.make_dir(dir_path, _counter = _counter+1)
            else:
                raise
            
        return new_dir_path