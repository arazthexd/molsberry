from ....core import (
    SimpleBlock, NumericData, FloatRep, NpData, NpArrayRep,
    StringData, StringRep
)

from ....core import Representation, generate_random_str

import numpy as np
from typing import Dict, List, Tuple  
import os

class Multiplier(SimpleBlock):
    name = "multiplier"
    display_name = "Numeric Multiplier"
    inputs = [
        ("num1", NumericData, FloatRep, False),
        ("num2", NumericData, FloatRep, False)
    ]
    outputs = [
        ("num_out", NumericData, FloatRep, False)
    ]
    batch_groups = []

    def __init__(self, debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        
    def multiplier(self, num1, num2):
        return num1 * num2

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        num1 = input_dict[self.input_keys[0]].content
        num2 = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        calc = self.multiplier(num1, num2)
        solv = generate_random_str(6)
        text_solv = os.path.join(self.base_dir, f"{solv}.txt")
        with open(text_solv, 'w') as f:
            f.write(str(calc))
        return {main_out_key: self._get_out_rep(main_out_key)(calc)}
    
class Subtractor(SimpleBlock):
    name = "subtractor"
    display_name = "Numeric Subtractor"
    inputs = [
        ("num1", NumericData, FloatRep, False),
        ("num2", NumericData, FloatRep, False)
    ]
    outputs = [
        ("num_out", NumericData, FloatRep, False)
    ]

    @property
    def batch_groups(self) -> List[Tuple[str]]:
        return [("num1", "num2")]

    def __init__(self, debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        
    def subtractor(self, num1, num2):
        return num1 - num2

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        num1 = input_dict[self.input_keys[0]].content
        num2 = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        calc = self.subtractor(num1, num2)
        solv = generate_random_str(6)
        text_solv = os.path.join(self.base_dir, f"{solv}.txt")
        with open(text_solv, 'w') as f:
            f.write(str(calc))
        return {main_out_key: self._get_out_rep(main_out_key)(calc)}

class Adder(SimpleBlock):  
    name = "adder"  
    display_name = "Numeric Adder"  
    
    outputs = [  
        ("num_out", NumericData, FloatRep, False)  
    ]  
    
    def __init__(self, num_inputs: int, debug=False, save_output=False, num_workers=None):  
        super().__init__(debug, save_output, num_workers)  
        self.num_inputs = num_inputs  
        self._inputs = self._generate_inputs()  

    @property  
    def inputs(self) -> List[tuple]:  
        return self._inputs  
    
    @property
    def batch_groups(self) -> List[Tuple[str]]:
        return [tuple(f"num{i+1}" for i in range(self.num_inputs))]

    def _generate_inputs(self) -> List[tuple]:  
        return [(f"num{i+1}", NumericData, FloatRep, False) 
                for i in range(self.num_inputs)] 


    def adder(self, *args): 
        return sum(args)

    def operate(self, input_dict: Dict[str, Representation]) \
            -> Dict[str, Representation]:  
        
        numbers = [input_dict[self.input_keys[i]].content 
                   for i in range(self.num_inputs)]  
        main_out_key = self.output_keys[0]  

        calc = self.adder(*numbers)  
 
        solv = generate_random_str(6)  
        text_solv = os.path.join(self.base_dir, f"{solv}.txt")  
 
        with open(text_solv, 'w') as f:  # TODO save for rep
            f.write(str(calc))  

        return {main_out_key: self._get_out_rep(main_out_key)(calc)}  

    
class Divider(SimpleBlock):
    name = "divider"
    display_name = "Numeric divider"
    inputs = [
        ("num1", NumericData, FloatRep, False),
        ("num2", NumericData, FloatRep, False)
    ]
    outputs = [
        ("num_out", NumericData, FloatRep, False)
    ]
    batch_groups = []

    def __init__(self, debug = False, save_output = False, num_workers = None):
        super().__init__(debug, save_output, num_workers)
        
    def divider(self, num1, num2):
        return num1 / num2

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        num1 = input_dict[self.input_keys[0]].content
        num2 = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        calc = self.divider(num1, num2)
        solv = generate_random_str(6)
        text_solv = os.path.join(self.base_dir, f"{solv}.txt")
        with open(text_solv, 'w') as f:
            f.write(str(calc))
        return {main_out_key: self._get_out_rep(main_out_key)(calc)}
    

