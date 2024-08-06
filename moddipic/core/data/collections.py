from typing import List, Type

from .abstract import Data

class Batched(Data):
    def __init__(self, data_list: List[Data]):
        self._data_list = data_list 
        self.obtain_batch_type()
        # TODO: Assert if with every branch of the batched tree, type is same.

    def obtain_batch_type(self):
        selected_data = self._data_list[0]
        if isinstance(selected_data, Batched):
            btype_child = selected_data.batch_type
        elif isinstance(selected_data, Data):
            btype_child = [type(selected_data)]
        else:
            raise ValueError("Batch should only consist of Data class.")

        self._batch_type: List[Type] = [Batched] + btype_child
    
    @property
    def batch_type(self):
        return self._batch_type # TODO: return in a format easier to understand
    
    @property
    def shallow_dtype(self):
        return self._batch_type[1]
    
    @property
    def basic_dtype(self):
        return self._batch_type[-1]
    
    @property
    def depth(self):
        return len(self._batch_type) - 1
    
    def __len__(self):
        return len(self._data_list)
    
    def __iter__(self):
        return iter(self._data_list)
    
    def __getitem__(self, index):
        return self._data_list[index]

#     def get_basic_data_type(self):
#         for i, element in enumerate(self.data):
#             if isinstance(element, Batched):
#                 el_type = element.get_basic_data_type()
#             else:
#                 el_type = type(element)
#             if i == 0:
#                 ref_el_type = el_type
#                 continue
#             # print(ref_el_type, el_type)
#             assert ref_el_type == el_type
#         return ref_el_type
    
#     def get_depth(self):
#         max_depth = 0
#         for i, element in enumerate(self.data):
#             if isinstance(element, Batched):
#                 depth = 1 + element.get_depth()
#             else:
#                 depth = 0
#             max_depth = max(max_depth, depth)
#         return max_depth
      
    
    
    

    

# class Grouped:
#     def __init__(self, members: list):
#         assert isinstance(members, list)
#         self.members = members
#     # TODO: Complete Grouped if needed...
    
