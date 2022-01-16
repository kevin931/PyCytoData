from PyCytoData import FileIO, CytoData, DataLoader, exceptions
import numpy as np
import pandas as pd

import pytest
import os
import sys
from io import StringIO
import shutil
import csv
import _csv

from typing import List, Any, Dict, Tuple, Literal


class TestCytoData():
    
    @classmethod
    def setup_class(cls):
        expression_matrix: np.ndarray = np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])
        features: np.ndarray = np.array(["feature1", "feature2", "feature3"])
        cell_types: np.ndarray = np.array(["TypeA", "TypeB"])
        sample_index: np.ndarray = np.array(["SampleA", "SampleA"])
        cls.dataset = CytoData(expression_matrix=expression_matrix, cell_types=cell_types, features=features, sample_index=sample_index)
        
        
    @pytest.mark.parametrize("attr,expected",
            [("expression_matrix", np.ndarray),
             ("features", np.ndarray),
             ("cell_types", np.ndarray),
             ("sample_index", np.ndarray),
             ("n_features", int),
             ("n_cells", int),
             ("n_samples", int),
             ("n_cell_types", int)])
    def test_attributes_type(self, attr: str, expected: Any):
        assert isinstance(getattr(self.dataset, attr), expected)
        
        
    @pytest.mark.parametrize("attr,expected",
            [("n_features", 3),
             ("n_cells", 2),
             ("n_samples", 1),
             ("n_cell_types", 2)])
    def test_metadata_attributes_values(self, attr: str, expected: Any):
        assert getattr(self.dataset, attr) == expected
        
    
    @pytest.mark.parametrize("attr,expected",
        [("expression_matrix", np.array([[1.1, 2.2, 3.3], [4.4, 5.5, 6.6]])),
         ("features", np.array(["feature1", "feature2", "feature3"])),
         ("cell_types", np.array(["TypeA", "TypeB"])),
         ("sample_index", np.array(["SampleA", "SampleA"]))])
    def test_array_attributes_values(self, attr: str, expected: Any):
        assert np.all(getattr(self.dataset, attr) == expected)
        
        
    def test_add_sample(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        cell_types: np.ndarray = np.array(["TypeC", "TypeD"])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB"])
        self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index, cell_types=cell_types)
        assert self.dataset.expression_matrix.shape == (4,3)
        assert self.dataset.cell_types.shape[0] == 4
        assert self.dataset.sample_index.shape[0] == 4
        assert "SampleB" in self.dataset.sample_index
        assert "TypeC" in self.dataset.cell_types
        assert "TypeD" in self.dataset.cell_types
        
        
    def test_add_sample_sample_index_error(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        cell_types: np.ndarray = np.array(["TypeC", "TypeD"])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index, cell_types=cell_types)
        except exceptions.DimensionMismatchError as e:
            message = str(e)
            assert "The `sample_index` attribute has to be of length 2." in message
    
    
    def test_add_sample_cell_types_error(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        cell_types: np.ndarray = np.array(["TypeC", "TypeD", "TypeD"])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index, cell_types=cell_types)
        except exceptions.DimensionMismatchError as e:
            message = str(e)
            assert "The `cell_types` attribute has to be of length 2." in message
            
            
    def test_add_sample_expression_matrix_dimension_error(self):
        new_sample: np.ndarray = np.array([1.2, 2.3, 3.4])
        sample_index: np.ndarray = np.array(["SampleB", "SampleB"])
        try:
            self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index)
        except exceptions.ExpressionMatrixDimensionError as e:
            message = str(e)
            assert "The shape (3,) is unsupported. Please reshape it to or use a two-dimensional array." in message
            
            
    def test_add_sample_no_cell_types(self):
        new_sample: np.ndarray = np.array([[1.2, 2.3, 3.4], [4.5, 5.6, 6.7]])
        sample_index: np.ndarray = np.array(["SampleC", "SampleC"])
        self.dataset.add_sample(expression_matrix=new_sample, sample_index=sample_index)
        assert None in self.dataset.cell_types
        
        
    def test_setters(self):
        expression_matrix = np.array([[1.2, 2.3, 3.4, 4.5], [4.5, 5.6, 6.7, 7.8]])
        self.dataset.expression_matrix = expression_matrix
        self.dataset.sample_index = np.array(["SampleA", "SampleB"])
        self.dataset.cell_types = np.array(["TypeA", "TypeB"])
        assert self.dataset.n_cells == 2
        assert self.dataset.n_features == 4
        assert self.dataset.n_samples == 2
        assert self.dataset.n_cell_types == 2
        
        
    def test_setters_n_values(self):
        self.dataset.n_cells = 2
        self.dataset.n_features = 4
        self.dataset.n_samples = 2
        self.dataset.n_cell_types = 2
        assert self.dataset.n_cells == 2
        assert self.dataset.n_features == 4
        assert self.dataset.n_samples == 2
        assert self.dataset.n_cell_types == 2
    
    
    @pytest.mark.parametrize("attr",
        ["n_features","n_cells","n_samples"])
    def test_setter_type_error(self, attr: str):
        try:
            setattr(self.dataset, attr, "test")
        except TypeError as e:
            f"'{attr}' has to be 'int' instead of str" in str(e)
        
        
    @pytest.mark.parametrize("attr,input_array,length",
        [("sample_index", np.array(["a", "a", "a"]), 2),
         ("cell_types", np.array(["a", "a", "a"]), 2),
         ("features", np.array(["a", "a", "a"]), 4)])
    def test_setter_dimension_mismatch_error(self, attr: str, input_array: np.ndarray, length: int):
        try:
            setattr(self.dataset, attr, input_array)
        except exceptions.DimensionMismatchError as e:
            f"The `{attr}` attribute has to be of length {length}."
    
    
    @classmethod
    def teardown_class(cls):
        pass
        

class TestDataLoader():
    
    @classmethod
    def setup_class(cls):
        os.makedirs("./tmp_pytest/data/levine13")
        os.makedirs("./tmp_pytest/data/levine32")
        os.makedirs("./tmp_pytest/data/samusik")
            
        for dataset in ["levine13", "levine32", "samusik"]:
            exprs_path: str = "./tmp_pytest/data/"+ dataset + "/" + dataset +".txt"
            types_path: str = "./tmp_pytest/data/"+ dataset + "/" + dataset +"_cell_types.txt"
            
            with open(exprs_path, "w") as f:
                tsv_writer: "_csv._writer" = csv.writer(f, delimiter="\t")
                tsv_writer.writerow(["CD1", "CD2", "CD3"])
                tsv_writer.writerow([4.4, 5.5, 6.6])
                tsv_writer.writerow([1.2, 1.3, 1.4])
                
            with open(types_path, "w") as f:
                tsv_writer: "_csv._writer" = csv.writer(f, delimiter="\t")
                tsv_writer.writerow(["B"])
                tsv_writer.writerow(["A"])
        
        with open("./tmp_pytest/data/samusik/population_assignments.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter="\t")
            tsv_writer.writerow(["a_a_a_01_a.fcs Event 000", "TypeA"])
            tsv_writer.writerow(["a_a_a_01_a.fcs Event 001", "TypeB"])
            tsv_writer.writerow(["a_a_a_02_a.fcs Event 000", "TypeA"])


    def test_preprocess_levine32(self):
        exprs: pd.DataFrame = pd.DataFrame({"FeatureA": [1.1, 2.2, 3.3],
                                            "FeatureB": [4.4, 5.5, 6.6]})
        fcss: List[str] = ["a-a-a-AML08-a_a_a_TypeA_a.fcs", "a-a-a-AML09-a_a_a_TypeB_a.fcs"]
        sample_length: List[int] = [2, 1]
        meta: Dict[str, Any] = {}
        
        data: CytoData = DataLoader._preprocess_levine32(fcss, exprs,meta, sample_length)
        assert isinstance(data, CytoData)
        assert data.n_cells == 3
        assert data.n_features == 2
        assert data.n_samples == 2


    def test_preprocess_levine13(self):
        exprs: pd.DataFrame = pd.DataFrame({"FeatureA": [1.1, 2.2, 3.3],
                                            "FeatureB": [4.4, 5.5, 6.6]})
        fcss: List[str] = ["a_a_a_TypeA_a.fcs", "a_a_a_TypeB_a.fcs"]
        sample_length: List[int] = [2, 1]
        meta: Dict[str, Any] = {"_channels_": {"$PnN": pd.Series(["FeatureA", "FeatureB"])}}
        
        data: CytoData = DataLoader._preprocess_levine13(fcss, exprs,meta, sample_length)
        assert isinstance(data, CytoData)
        assert data.n_cells == 3
        assert data.n_features == 2
        assert data.n_samples == 1
        
        
    def test_preprocess_samusik(self, mocker):
        exprs: pd.DataFrame = pd.DataFrame({"FeatureA": [1.1, 2.2, 3.3],
                                            "FeatureB": [4.4, 5.5, 6.6]})
        fcss: List[str] = ["a_a_a_01_a.fcs", "a_a_a_02_a.fcs"]
        sample_length: List[int] = [2, 1]
        meta: Dict[str, Any] = {}
        
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/" + "levine13/",
                                                          "levine32": "./tmp_pytest/data/" + "levine32/",
                                                          "samusik": "./tmp_pytest/data/" + "samusik/"})
        
        data: CytoData = DataLoader._preprocess_samusik(fcss, exprs,meta, sample_length)
        assert isinstance(data, CytoData)
        assert data.n_cells == 3
        assert data.n_features == 2
        assert data.n_samples == 2
        
    
    @pytest.mark.parametrize("dataset, expected",
                             [("levine13", "levine13"),
                              ("levine32", "levine32"),
                              ("samusik", "samusik")]
                             )   
    def test_preprocess(self, mocker, dataset: Literal["levine13", "levine32", "samusik"], expected: str):
        fcss: List[str] = ["a_a_a_01_a.fcs", "a_a_a_02_a.fcs"]
        exprs: pd.DataFrame = pd.DataFrame({"FeatureA": [1.1, 2.2, 3.3],
                                            "FeatureB": [4.4, 5.5, 6.6]})
        fcs: Tuple[Dict, pd.DataFrame] = ({}, exprs)
        mocker.patch("PyCytoData.data.glob.glob", return_value=fcss)
        mocker.patch("PyCytoData.data.fcsparser.parse",return_value=fcs)
        mocker.patch("PyCytoData.DataLoader._preprocess_levine13", return_value="levine13")
        mocker.patch("PyCytoData.DataLoader._preprocess_levine32", return_value="levine32")
        mocker.patch("PyCytoData.DataLoader._preprocess_samusik", return_value="samusik")
        
        out: Any = DataLoader._preprocess(dataset)
        assert out == expected
        
    
    @pytest.mark.parametrize("dataset",
                            ["levine13", "levine32", "samusik"])     
    def test_load_dataset(self, mocker, dataset: Literal["levine13", "levine32", "samusik"]):
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/" + "levine13/",
                                                          "levine32": "./tmp_pytest/data/" + "levine32/",
                                                          "samusik": "./tmp_pytest/data/" + "samusik/"})
        mocker.patch("PyCytoData.DataLoader._data_status", {dataset: True})
        
        data: CytoData = DataLoader.load_dataset(dataset=dataset)
        assert isinstance(data, CytoData)
        assert None not in data.cell_types
        assert data.n_cells == 2
        assert data.n_features == 3
       
        
    def test_load_dataset_no_cell_types(self, mocker): 
        os.remove("./tmp_pytest/data/levine13/levine13_cell_types.txt")
        mocker.patch("PyCytoData.DataLoader._data_status", {"levine13": True})
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/" + "levine13/",
                                                          "levine32": "./tmp_pytest/data/" + "levine32/",
                                                          "samusik": "./tmp_pytest/data/" + "samusik/"})
        
        data: CytoData = DataLoader.load_dataset(dataset="levine13")
        assert isinstance(data, CytoData)
        assert None in data.cell_types
        assert data.n_cells == 2
        assert data.n_features == 3
        
    
    @pytest.mark.parametrize("force_download",
                            [True, False]
                            )   
    def test_load_dataset_download_dataset(self, mocker, force_download: bool):         
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        mocker.patch("PyCytoData.DataLoader._data_path", {"levine13": "./tmp_pytest/data/" + "levine13/",
                                                          "levine32": "./tmp_pytest/data/" + "levine32/",
                                                          "samusik": "./tmp_pytest/data/" + "samusik/"})
        mocker.patch("PyCytoData.DataLoader._data_status", {"levine13": False})
        mock_download: mocker.MagicMock = mocker.MagicMock()
        mocker.patch("PyCytoData.DataLoader._download_data", mock_download)
        
        data: CytoData = DataLoader.load_dataset(dataset="levine13", force_download=force_download)
        assert isinstance(data, CytoData)
        assert None in data.cell_types
        assert data.n_cells == 2
        assert data.n_features == 3
        mock_download.assert_called_with(dataset="levine13", force_download=force_download)

            
    @pytest.mark.parametrize("input_value",
                             ["n", "N"]
                             )     
    def test_download_assets_force_and_input_no(self, mocker, input_value):
        
        screen_stdout = sys.stdout
        string_stdout = StringIO()
        sys.stdout = string_stdout
        
        mocker.patch("builtins.input", return_value=input_value)
        out:int = DataLoader._download_data(dataset="levine13", force_download=False)
        
        output = string_stdout.getvalue()
        sys.stdout = screen_stdout
        
        assert "You have declined to download" in output    
        assert "Download in progress" not in output
        assert out == 1   
            
    
    @pytest.mark.parametrize("input_value,force_download",
                            [("y", False),
                             ("Y", False),
                             ("y", True)])   
    def test_download_assets_url_download_extract_mock(self, mocker, input_value: str, force_download: bool):
        contents_mock = mocker.MagicMock()
        contents_mock.read.return_value = b"This is a test."
        zip_mock = mocker.MagicMock()
        zip_mock.extractall.return_value = None
        
        mocker.patch("builtins.input", return_value=input_value)
        mocker.patch("PyCytoData.data.urlopen", return_value = contents_mock)
        mocker.patch("PyCytoData.data.ZipFile", return_value = zip_mock)
        mocker.patch("PyCytoData.DataLoader._data_dir", "./tmp_pytest/data/")
        expression_matrix: np.ndarray = np.array([[1, 2],[2, 3]])
        data_mock: CytoData = CytoData(expression_matrix)
        mocker.patch("PyCytoData.DataLoader._preprocess", return_value=data_mock)
        save_array_mock: mocker.MagicMock = mocker.MagicMock()
        mocker.patch("PyCytoData.FileIO.save_np_array", save_array_mock)
        
        out: int = DataLoader._download_data(dataset="levine13", force_download=force_download)
        
        save_array_mock.assert_called()
        zip_mock.extractall.assert_called_once()
        contents_mock.read.assert_called_once()
        assert out == 0
    
    
    @classmethod
    def teardown_class(cls):
        shutil.rmtree("./tmp_pytest/")
        

class TestFileIO():
    
    @classmethod
    def setup_class(cls):
        os.mkdir("./tmp_pytest/")
        
        with open("./tmp_pytest/file_read_test_csv.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter=",")
            tsv_writer.writerow(["col1", "col2", "col3"])
            tsv_writer.writerow([1, 2, 3])
            tsv_writer.writerow([4, 5, 6])
            
        with open("./tmp_pytest/file_read_test_tsv.txt", "w") as f:
            tsv_writer: "_csv._writer" = csv.writer(f, delimiter="\t")
            tsv_writer.writerow([1.1, 2.2, 3.3])
            tsv_writer.writerow([4.4, 5.5, 6.6])
        
    
    @pytest.mark.parametrize("path,delim,col_names,dtype",
                [("./tmp_pytest/file_read_test_csv.txt", ",", True, int),
                 ("./tmp_pytest/file_read_test_tsv.txt", "\t", False, float)]
                )
    def test_load_delim_filetype(self, path: str, delim: str, col_names: bool, dtype):
        out_file: CytoData = FileIO.load_delim(path, col_names, delim=delim, dtype = dtype)
        assert isinstance(out_file, CytoData)
        assert out_file.n_cells == 2
        assert out_file.n_features == 3
        
        if path == "./tmp_pytest/file_read_test_csv.txt":
            assert "col1" in out_file.features
            assert out_file.expression_matrix.dtype == np.dtype("int64")
        else:
            assert None in out_file.features
            assert out_file.expression_matrix.dtype == np.dtype("float64")
            
    
    @pytest.mark.parametrize("drop_cols,expected_shape",
            [([0, 1], (2, 1)), (1, (2,2))]
            )      
    def test_load_delim_drop_col(self, drop_cols, expected_shape):
        path: str = "./tmp_pytest/file_read_test_csv.txt"
        out_file: CytoData = FileIO.load_delim(path, True, drop_columns=drop_cols, delim=",", dtype = int)
        assert out_file.expression_matrix.shape == expected_shape
        
        
    @pytest.mark.parametrize("save_list,path",
                    [([[1, 2, 3], [1, 2, 3]], "./tmp_pytest/2d_list_1.csv"),
                    ([[1.4, 2.2], ["a", "b"]], "./tmp_pytest/2d_list_2.csv"),
                    ([[1.4, 1], [1, "a"]], "./tmp_pytest/2d_list_3.csv")]
                    )   
    def test_save_2d_list_to_csv(self, save_list: List[List[Any]], path: str):
        FileIO.save_2d_list_to_csv(save_list, path)
        assert os.path.exists(path)


    def test_save_2d_list_to_csv_exception(self):
        path: str = "./tmp_pytest/2d_list_3.csv"
        save_list: List[List[Any]] = [[1.4, 1], [1, "a"]]
        try:
            FileIO.save_2d_list_to_csv(save_list, path)
        except FileExistsError:
            assert True
            

    def test_save_np_array(self):
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        path: str = "./tmp_pytest/nparr.txt"
        FileIO.save_np_array(arr, path)
        assert os.path.exists("./tmp_pytest/nparr.txt")
        
        
    def test_save_np_array_colnames(self): 
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        header: "np.ndarray" = np.array(["A", "B"])
        path: str = "./tmp_pytest/nparr_colnames.txt"
        FileIO.save_np_array(arr, path, col_names=header)
        assert os.path.exists("./tmp_pytest/nparr_colnames.txt")
        
    
    def test_save_np_array_exception(self):
        arr: "np.ndarray" = np.array([[2.1, 2.2, 2.3], [1.1, 2.2, 3.3]])
        path: str = "./tmp_pytest/nparr.txt"
        try:
            FileIO.save_np_array(arr, path)
        except FileExistsError:
            assert True
        
        
    def test_dir_create(self):
        FileIO.make_dir("./tmp_pytest/test_create")
        assert os.path.exists("./tmp_pytest/test_create")
        
     
    @pytest.mark.parametrize("path_in,counter,path_out",
                        [("./tmp_pytest/test_recur", 0, "./tmp_pytest/test_recur1"),
                        ("./tmp_pytest/test_recur_counter", 2, "./tmp_pytest/test_recur_counter2")]
                        )   
    def test_dir_create_recursive(self, path_in: str, counter: int, path_out: str):
        os.mkdir(path_in)
        FileIO.make_dir(path_in, add_number_if_dir_exists=True, _counter=counter)
        assert os.path.exists(path_out)
    
    
    @pytest.mark.parametrize("path_in,recursive,expected",
                    [("./tmp_pytest/test_recur_out", False, "./tmp_pytest/test_recur_out"),
                    ("./tmp_pytest/test_recur_out", True, "./tmp_pytest/test_recur_out1")]
                    )  
    def test_dir_create_return(self, path_in: str, recursive: bool, expected: str):
        path_out: str = FileIO.make_dir(path_in, add_number_if_dir_exists=recursive)
        assert path_out == expected
    
    
    def test_dir_create_exception(self):
        try:
            FileIO.make_dir("./tmp_pytest/test_create")
        except FileExistsError:
            assert True
    
    
    @classmethod
    def teardown_class(cls):
        shutil.rmtree("./tmp_pytest/")