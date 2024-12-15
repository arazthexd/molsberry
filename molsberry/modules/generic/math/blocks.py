from ....core import (
    SimpleBlock, NumericData, FloatRep
)

from .mathrep import NumRep
from ....core import Representation, generate_random_str

from typing import Dict
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
    batch_groups = []

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
        
    def adder(self, num1, num2):
        return num1 + num2

    def operate(self, input_dict: Dict[str, Representation]) \
        -> Dict[str, Representation]:
        num1 = input_dict[self.input_keys[0]].content
        num2 = input_dict[self.input_keys[1]].content
        main_out_key = self.output_keys[0]
        calc = self.adder(num1, num2)
        solv = generate_random_str(6)
        text_solv = os.path.join(self.base_dir, f"{solv}.txt")
        with open(text_solv, 'w') as f:
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