class Batched:
    def __init__(self, data: list):
        assert isinstance(data, list)
        self.data = data

    def get_basic_data_type(self):
        for i, element in enumerate(self.data):
            if isinstance(element, Batched):
                el_type = element.get_basic_data_type()
            else:
                el_type = type(element)
            if i == 0:
                ref_el_type = el_type
                continue
            print(ref_el_type, el_type)
            assert ref_el_type == el_type
        return ref_el_type
    
    def get_depth(self):
        max_depth = 0
        for i, element in enumerate(self.data):
            if isinstance(element, Batched):
                depth = 1 + element.get_depth()
            else:
                depth = 0
            max_depth = max(max_depth, depth)
        return max_depth
      
    def __getitem__(self, index):
        return self.data[index]
    
    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

class Grouped:
    def __init__(self, members: list):
        assert isinstance(members, list)
        self.members = members
    # TODO: Complete Grouped if needed...
    
