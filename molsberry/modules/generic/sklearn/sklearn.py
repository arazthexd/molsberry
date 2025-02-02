from ....core import (
    SimpleBlock, NumericData, FloatRep, NpData, NpArrayRep,
    StringData, StringRep
)

from ....core import Representation, generate_random_str

import numpy as np
from typing import Dict, List, Tuple  
import os

import joblib 
from sklearn.linear_model import LinearRegression

class NumericArrayGenerator(SimpleBlock):  
    name = "numeric_array_generator"  
    display_name = "Numeric Array Generator"  
    
    outputs = [  
        ("array_out", NpData, NpArrayRep, False)  
    ]  
    batch_groups = []  
    
    def __init__(self, num_inputs: int, debug=False, save_output=False, num_workers=None):  
        super().__init__(debug, save_output, num_workers)  
        self.num_inputs = num_inputs  
        self._inputs = self._generate_inputs()  

    @property  
    def inputs(self) -> List[tuple]:  
        return self._inputs  

    def _generate_inputs(self) -> List[tuple]:  
        return [(f"num{i+1}", NumericData, FloatRep, False) for i in range(self.num_inputs)]  

    def create_numpy_array(self, input_dict: Dict[str, FloatRep]) -> np.ndarray:  
        """Creates a NumPy array from the input dictionary."""  
        array_data = []  
        for i in range(self.num_inputs):  
            key = self.input_keys[i]  
            array_data.append(input_dict[key].content) 

        return np.array(array_data)  

    def operate(self, input_dict: Dict[str, Representation]) -> Dict[str, Representation]:  
        """Process inputs and generate a NumPy array."""  
        numpy_array = self.create_numpy_array(input_dict) 

        main_out_key = self.output_keys[0]  

        return {main_out_key: self._get_out_rep(main_out_key)(numpy_array)}  
    
class CorrCoef(SimpleBlock):
    name = "CorrCoefFinder"
    display_name = "CorrCoeff calculator"
    inputs = [
        ("npa1", NpData, NpArrayRep, False),
        ("npa2", NpData, NpArrayRep, False)
    ]
    outputs = [
        ("num_out", NumericData, FloatRep, False)
    ]
    batch_groups = []

    def __init__(self, debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        
    def coorelationfind(self, npa1, npa2):
        correlation_matrix = np.corrcoef(npa1, npa2)  
        return correlation_matrix[0, 1]  

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        npa1 = input_dict[self.input_keys[0]].content
        npa2 = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        calc = self.coorelationfind(npa1, npa2)
        solv = generate_random_str(6)
        text_solv = os.path.join(self.base_dir, f"{solv}.txt")
        with open(text_solv, 'w') as f:
            f.write(str(calc))
        return {main_out_key: self._get_out_rep(main_out_key)(calc)}

class SKLinearRegression(SimpleBlock):
    name = "SKLinearRegressionModel"
    display_name = "sk-learn linear regression model"
    inputs = [
        ("npa1", NpData, NpArrayRep, False),
        ("npa2", NpData, NpArrayRep, False)
    ]
    outputs = [
        ("model", StringData, StringRep, False)
    ]
    batch_groups = []

    def __init__(self, debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        
    def linearregression(self, npa1, npa2, model):
        npa1 = npa1.reshape(-1, 1) if npa1.ndim == 1 else npa1  
        self.model = LinearRegression()  
        self.model.fit(npa1, npa2)  
        joblib.dump(self.model, model)  
        return model

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        npa1 = input_dict[self.input_keys[0]].content
        npa2 = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        model = generate_random_str(6)
        model = os.path.join(self.base_dir, f"{model}.pkl")
        calc = self.linearregression(npa1, npa2, model)
        return {main_out_key: self._get_out_rep(main_out_key)(calc)}
    
class SKPredictor(SimpleBlock):
    name = 'Skpredictor'
    display_name = 'Sk-learn predictor with the model'
    inputs = [
        ('model', StringData, StringRep, False),
        ('num', NumericData, FloatRep, False)
    ]
    outputs = [
        ("num_out", NumericData, FloatRep, False)
    ]
    batch_groups = []

    def __init__(self, debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        
    def predictor(self, model, num):
        model = joblib.load(model)
        num = np.array(num).reshape(-1, 1)
        prediction = model.predict(num)
        return prediction

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        model = input_dict[self.input_keys[0]].content
        num = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        calc = self.predictor(model, num)
        return {main_out_key: self._get_out_rep(main_out_key)(calc)}