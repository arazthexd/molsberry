from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List, Tuple, Any, Dict

import os
import pathlib

class _Connection():
    def __init__(self, parent: Pipeline, source_block: PipelineBlock, source_key: str, 
                 target_block: PipelineBlock, target_key: str) -> None:
        source_block._add_connection(self)
        target_block._add_connection(self)
        source_block._cyclic_check_pass([])
        self.source_block = source_block
        self.source_key = source_key
        self.target_block = target_block
        self.target_key = target_key
        self.parent = parent
    
    def get_data_from_source(self, input_dict):
        self.source_block._auto_execute(input_dict)
        return self.source_block._output[self.source_key]

class PipelineBlock(ABC):
    def __init__(self, debug: bool = False) -> None:
        self.debug = debug
        self._executed = False
        self._latest_input = None
        self._output = dict()
        self._in_connections: Dict[str, _Connection] = dict()
        self._out_connections: Dict[str, List[_Connection]] = dict()
        self._parent = None
        self._name_in_parent = None
    
    def reset(self, removeOutput: bool = True, removeConnections: bool = True) -> None:
        if removeOutput:
            self._executed = False
            self._latest_input = None
            self._output = dict()
        if removeConnections:
            self._in_connections: Dict[str, _Connection] = dict()
            self._out_connections: Dict[str, List[_Connection]] = dict()
    
    def _add_connection(self, connection: _Connection) -> None:
        if connection.source_block == self:
            assert connection.source_key in self.output_keys
            if connection.source_key not in self._out_connections.keys():
                self._out_connections[connection.source_key] = []
            self._out_connections[connection.source_key].append(connection)
        if connection.target_block == self:
            assert connection.target_key in self.required_input_keys + self.optional_input_keys
            assert connection.target_key not in self._in_connections.keys()
            self._in_connections[connection.target_key] = connection
    
    def _next_blocks(self) -> List[PipelineBlock]:
        next_blocks = []
        for connections in self._out_connections.values():
            for connection in connections:
                if connection.target_block not in next_blocks:
                    next_blocks.append(connection.target_block)
        return next_blocks
        
    def _cyclic_check_pass(self, prev_blocks: List[PipelineBlock]) -> None:
        assert self not in prev_blocks
        prev_blocks.append(self)
        for block in self._next_blocks():
            block._cyclic_check_pass(prev_blocks)

    @abstractmethod
    def execute(self, input_dict: Dict[str, Any]) -> Dict[str, Any]:
        pass
    
    def _auto_execute(self, pipeline_dict: Dict[str, Any]) -> None:
        connection_incoming_input = {name: connection.get_data_from_source(pipeline_dict) 
                                     for name, connection in self._in_connections.items()}
        block_input = pipeline_dict.copy()
        block_input.update(connection_incoming_input)
        
        if self._executed and self._latest_input == block_input:
            return
        
        print()
        print(f"[Running Pipe Block: ({self._name_in_parent}) {self.name}]")
        self.output = self.execute(block_input)

    @property
    def output(self) -> dict:
        return self._output

    @output.setter
    def output(self, d: dict):
        assert set(d.keys()) == set(self.output_keys)
        self._output = d

    @property
    def executed(self) -> bool:
        return self._executed

    @property
    @abstractmethod
    def output_keys(self) -> List[str]:
        pass
    
    @property
    @abstractmethod
    def required_input_keys(self) -> List[str]:
        pass
    
    @property
    @abstractmethod
    def optional_input_keys(self) -> List[str]:
        pass
        
    @property
    @abstractmethod
    def name(self) -> str:
        pass

class IOBlock(PipelineBlock):
    def __init__(self, keys: List[str], debug: bool = False) -> None:
        super().__init__(debug)
        self._keys = keys

    def execute(self, input_dict: Dict[str, Any]) -> Dict[str, Any]:
        return input_dict
    
    @property
    def output_keys(self) -> List[str]:
        return self._keys
    
    @property
    def required_input_keys(self) -> List[str]:
        return self._keys
    
    @property
    def optional_input_keys(self) -> List[str]:
        return self._keys
        
    @property
    def name(self) -> str:
        return "IO Block"
    
class Pipeline(ABC):
    def __init__(self, base_dir: str = ".", initiate: bool = True):
        if not os.path.exists(base_dir):
            os.mkdir(base_dir)
        
        self._initiated = False
        self.base_dir = pathlib.Path(base_dir)
        self._blocks: Dict[str, PipelineBlock] = dict()
        self._connections: Dict[str, _Connection] = dict()

        # helper
        self._connect_naming_i = 0
        self._block_naming_i = 0

        # for initiation
        self._required_input_keys = set()
        self._optional_input_keys = set()

        self.build()
        if initiate:
            self.initiate()
    
    def reset(self):
        raise NotImplementedError()
    
    def add_block(self, block: PipelineBlock, block_name: str | None = None) -> None:
        assert not self._initiated
        if block_name is None:
            block_name = "block{i}"
            while True:
                self._block_naming_i += 1
                name_suggest = block_name.format(i=self._block_naming_i)
                if name_suggest not in self._blocks.keys():
                    block_name = name_suggest
                    break
        assert block_name not in self._blocks.keys()
        
        block._parent = self
        block._name_in_parent = block_name # TODO
        self._blocks[block_name] = block
    
    def add_blocks(self, blocks: List[PipelineBlock] | Dict[str, PipelineBlock]) -> None:
        if isinstance(blocks, list):
            [self.add_block(block) for block in blocks]
        if isinstance(blocks, dict):
            [self.add_block(block, name) for name, block in blocks.items()]
    
    def add_connection(self, source_name: str, source_key: str, target_name: str, target_key: str,
                       connection_name: str | None = None) -> None:
        assert not self._initiated
        if connection_name is None:
            connection_name = "connection{i}"
            while True:
                self._connect_naming_i += 1
                name_suggest = connection_name.format(i=self._connect_naming_i)
                if name_suggest not in self._connections.keys():
                    connection_name = name_suggest
                    break
        assert connection_name not in self._connections.keys()
        assert source_name in self._blocks.keys()
        assert target_name in self._blocks.keys()
        
        connection = _Connection(self, self._blocks[source_name], source_key, 
                                 self._blocks[target_name], target_key)
        self._connections[connection_name] = connection

    @abstractmethod
    def build(self):
        pass

    def initiate(self):
        for name, block in self._blocks.items():
            block_incoming = block._in_connections.values()
            [self._required_input_keys.add(inp) for inp in block.required_input_keys if inp not in block_incoming]
            [self._optional_input_keys.add(inp) for inp in block.optional_input_keys if inp not in block_incoming]
        
        assert "output" not in self._blocks.keys()
        self.add_block(IOBlock(self.output_keys), "output")
        [self.add_connection(*okey.split("."), "output", okey) for okey in self.output_keys]

        self._initiated = True
    
    def execute(self, input_dict: Dict[str, Any]) -> Dict[str, Any]:
        assert self._initiated

        l = 55
        print()
        print(">"*l)
        print(">>" + str.center("STARTED: " + self.name, l-5), ">>")
        print(">"*l)

        self._blocks["output"]._auto_execute(input_dict)

        print()
        print("<"*l)
        print("<<" + str.center("ENDED: " + self.name, l-5), "<<")
        print("<"*l)

        return self._blocks["output"]._output

    def _auto_execute(self, pipeline_dict: Dict[str, Any]):
        raise NotImplementedError() # TODO: Implement a method to let pipelines act as a pipe block as well

    @property
    def initiated(self):
        return self._initiated
    
    @property
    def required_input_keys(self) -> List[str]:
        assert self._initiated
        return self._required_input_keys
    
    @property
    def optional_input_keys(self) -> List[str]:
        assert self._initiated
        return self._optional_input_keys
    
    @property
    @abstractmethod
    def output_keys(self) -> List[str]:
        pass

    @property
    @abstractmethod
    def name(self) -> str:
        return "pipeline"
    