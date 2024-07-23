from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List, Tuple, Any, Dict

import os
import pathlib
import pickle

class _Connection():
    def __init__(self, parent: Pipeline, 
                 source_block: PipelineBlock, source_key: str,
                 target_block: PipelineBlock, target_key: str) -> None:
        self.source_block = source_block
        self.source_key = source_key
        self.target_block = target_block
        self.target_key = target_key
        self.source_block._add_connection(self)
        self.target_block._add_connection(self)
        self.source_block._cyclic_check_pass([])
        
        self.parent = parent

    def get_data_from_source(self, 
                             global_dict: Dict[str, Any]) -> Dict[str, Any]:
        
        # Execute the source block and get return its output.
        # NOTE: If it has been executed it won't be executed again.
        self.source_block._auto_execute(global_dict)
        return self.source_block._output[self.source_key]

class PipelineBlock(ABC):
    # TODO: Do we need to add remove_connection to blocks?
    # TODO: Save directory property... vs Pipeline.
    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        self.debug = debug
        self.save_output = save_output

        # Whenever the block runs, _executed will be set to True. Also,
        # _latest_input and _output will be set accordingly. In the next run
        # if the input is the same as before, the output will be returned
        # without actually running the block.
        self._executed = False
        self._latest_input = None
        self._output = dict()

        # Every block will have only one input connection for each input key
        # but can have multiple output connections for each output key. The 
        # connection will have names. These could also be named automatically.
        self._in_connections: Dict[str, _Connection] = dict()
        self._out_connections: Dict[str, List[_Connection]] = dict()
        
        # Every block must be part of a `Pipeline`. The parent pipeline as
        # well as the name of the block in the parent blocks dict must be set
        # before attempting to execute the block.
        self._parent: Pipeline = None
        self._name_in_parent = None
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Name of the pipeline block."""
        pass

    @property
    @abstractmethod
    def required_input_keys(self) -> List[str]:
        """Required keys in block's input (keywords/dictionary)"""
        pass
    
    @property
    @abstractmethod
    def optional_input_keys(self) -> List[str]:
        """Optional keys in block's input (keywords/dictionary)"""
        pass

    @property
    @abstractmethod
    def output_keys(self) -> List[str]:
        """Keys in block's output dictionary."""
        pass

    def check_input(self, input_dict: Dict[str, Any]):
        """Checks validity of input dictionary (input to `execute`).

        Args:
            input_dict (Dict[str, Any]): input dict with needed input keys.
        """

        assert all(k in input_dict.keys() for k in 
                    self.required_input_keys + self.optional_input_keys)
        return
    
    def check_output(self, output_dict: Dict[str, Any]):
        """Checks validity of output dictionary (output from `execute`).

        Args:
            output_dict (Dict[str, Any]): output dict with needed output keys.
        """
        assert set(output_dict.keys()) == set(self.output_keys)
        return

    @abstractmethod
    def execute(self, **kwargs) -> Dict[str, Any]:
        """Abstract method for executing a block.

        This is the only method that needs defining for a new pipeline block.
        In the implementation, every input should be taken in the form of a 
        keyword and at the end a dictionary containing its outputs and their
        keys should be returned.

        Returns:
            Dict[str, Any]: Dictionary containing output keys
        """
        pass

    def save(self) -> str:

        pkl_path = os.path.join(self.base_dir, 
                                (self.name+".pkl").replace(" ", "_").lower())
        with open(pkl_path, "wb") as f:
            pickle.dump(self._output, f)
        
        return pkl_path
    
    def _auto_execute(self, global_dict: Dict[str, Any]) -> None:
        """Acts as a wrapper around `execute` for automatic runs.

        Args:
            global_dict (Dict[str, Any]): Keys and data not provided using the
                block's connections or overwritten by input connections.
        """

        # Initial Checks...
        assert isinstance(self._parent, Pipeline)
        assert isinstance(self._name_in_parent, str)

        # Given the global_dict, it updates it using any input coming from its 
        # input connections and then executes with that info.
        connection_incoming_input = {
            name: connection.get_data_from_source(global_dict) 
            for name, connection in self._in_connections.items()}
        block_input = global_dict.copy()
        block_input.update(connection_incoming_input)
        
        # Do not execute if input is the same as last time and it has been
        # executed previously.
        if self._executed and self._latest_input == block_input:
            self.save()
            return
        
        # Logging...
        # TODO: Use the logging module of python instead of print?
        print()
        print(f"[Running Pipe Block: ({self._name_in_parent}) {self.name}]")
        self.check_input(block_input)
        self._output = self.execute(**block_input)
        if self.save_output:
            self.save()
    
    def reset(self, removeExecutionHistory: bool = True, 
              removeConnections: bool = True) -> None:
        
        if removeExecutionHistory:
            self._executed = False
            self._latest_input = None
            self._output = dict()
        
        if removeConnections:
            self._in_connections: Dict[str, _Connection] = dict()
            self._out_connections: Dict[str, List[_Connection]] = dict()
    
    def _add_connection(self, connection: _Connection) -> None:

        # NOTE: This method is only meant to be called in the parent.
        
        # If the block is the source of the connection, it should have the 
        # source_key of the connection in its output keys. It should also be
        # added to the list of out connections for the specific key.
        if connection.source_block == self:
            assert connection.source_key in self.output_keys
            if connection.source_key not in self._out_connections.keys():
                self._out_connections[connection.source_key] = []
            self._out_connections[connection.source_key].append(connection)
        
        # If the block is the target of the connection, it should have the 
        # target_key of the connection in either its required or optional input
        # keys and there shouldn't be another connection for its specific key.
        if connection.target_block == self:
            assert connection.target_key in \
                self.required_input_keys + self.optional_input_keys
            assert connection.target_key not in self._in_connections.keys()
            self._in_connections[connection.target_key] = connection
    
    def _next_blocks(self) -> List[PipelineBlock]:

        # Just a helper method for finding what nodes/blocks come after this
        # specific block

        next_blocks = []
        for connections in self._out_connections.values():
            for connection in connections:
                if connection.target_block not in next_blocks:
                    next_blocks.append(connection.target_block)
        return next_blocks
        
    def _cyclic_check_pass(self, prev_blocks: List[PipelineBlock]) -> None:

        # Helper method to check if the block forms cycles with its previous
        # and next blocks. This should not happen.

        assert self not in prev_blocks
        prev_blocks.append(self)
        for block in self._next_blocks():
            block._cyclic_check_pass(prev_blocks)

    @property
    def output(self) -> dict:
        """Output of the latest execution."""
        return self._output

    @output.setter
    def output(self, d: dict):
        """Setter for `output`.

        Args:
            d (dict): output dictionary containing output keys and values.
        """
        self.check_output(d)
        self._output = d

    @property
    def executed(self) -> bool:
        """Whether the block has been previously executed or not."""
        return self._executed
    
    @property
    def base_dir(self) -> str:
        return self._parent.base_dir
    
class OutputBlock(PipelineBlock, ABC):
    name = "Output Block"
    def __init__(self, keys: List[str], debug: bool = False) -> None:
        super().__init__(debug)

        self.keys = keys
    
    @property
    def required_input_keys(self):
        return self.keys
    
    @property
    def optional_input_keys(self):
        return []
    
    @property
    def output_keys(self):
        return self.keys
    
    def execute(self, **kwargs) -> Dict[str, Any]:
        return kwargs
    
    # def save(self):

    #     pkl_path = os.path.join(self.base_dir, 
    #                             (self._parent.name+".pkl")
    #                             .replace(" ", "_").lower())
    #     with open(pkl_path, "wb") as f:
    #         pickle.dump(self._output, f)
        
    #     return pkl_path
    
class Pipeline(PipelineBlock, ABC):

    # NOTE: The reason for inheriting from `PipelineBlock` is to be able to use
    # pipelines as blocks for other pipelines as well.

    def __init__(self, base_dir: str = ".", initiate: bool = True, 
                 debug: bool = False, save: bool = False):
        super().__init__(debug=debug, save_output=False)
        
        self._initiated = False
        self._base_dir = pathlib.Path(base_dir)
        self._blocks: Dict[str, PipelineBlock] = dict()
        self._connections: Dict[str, _Connection] = dict()

        # helper
        self._connect_naming_i = 0
        self._block_naming_i = 0

        # for initiation
        self._required_input_keys = list()
        self._optional_input_keys = list()

        self.build()
        assert "output" in self._blocks.keys()
        if save:
            self._blocks["output"].save_output = True
        if initiate:
            self.initiate()
    
    @property
    @abstractmethod
    def name(self) -> str:
        return "pipeline"
    
    @abstractmethod
    def build(self):
        pass

    def save(self):
        raise NotImplementedError()
    
    def execute(self, **input_dict) -> Dict[str, Any]:
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

        self._output = self._blocks["output"]._output
        return self._output

    def initiate(self):

        for name, block in self._blocks.items():
            block_incoming = block._in_connections.values()
            [self._required_input_keys.append(inp) for inp in 
             block.required_input_keys if inp not in block_incoming]
            [self._optional_input_keys.append(inp) for inp in 
             block.optional_input_keys if inp not in block_incoming]
        
        self._required_input_keys = list(set(self._required_input_keys))
        self._optional_input_keys = list(set(self._optional_input_keys))

        self._initiated = True

    def reset(self):
        raise NotImplementedError()
    
    def add_block(
            self, block: PipelineBlock, block_name: str | None = None
        ) -> None:
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
    
    def add_blocks(
            self, blocks: List[PipelineBlock] | Dict[str, PipelineBlock]
        ) -> None:
        assert not self._initiated

        if isinstance(blocks, list):
            [self.add_block(block) for block in blocks]
        if isinstance(blocks, dict):
            [self.add_block(block, name) for name, block in blocks.items()]
    
    def add_connection(self, source_name: str, source_key: str, 
                       target_name: str, target_key: str,
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
    def output_keys(self) -> List[str]:
        return self._blocks["output"].output_keys
    
    @property
    def base_dir(self) -> str:
        if self._parent:
            return_dir = os.path.join(self._parent._base_dir,
                                      self.name.replace(" ", "_").lower())
        else:
            return_dir = self._base_dir
        
        if not os.path.exists(return_dir):
            os.mkdir(return_dir)
            # TODO: This gets called every time. It shouldn't be so...

        return return_dir