from __future__ import annotations
from abc import ABC, abstractmethod
from typing import List, Tuple, Any, Dict

import os
import pathlib
import pickle

from .data.abstract import Data

# TODO: Optional Input Keys? When to use?

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

    def get_data_from_source(self) -> Data:
        
        # Execute the source block and get return its output.
        # NOTE: If it has been executed it won't be executed again.
        self.source_block._auto_execute()
        data = self.source_block._output[self.source_key]
        assert isinstance(data, Data)
        return data

class PipelineBlock(ABC):
    force_run: bool = False
    # TODO: Do we need to add remove_connection to blocks?
    # TODO: Save directory property... vs Pipeline.
    def __init__(self, debug: bool = False, save_output: bool = False) -> None:
        self.debug: bool = debug
        self.save_output: bool = save_output

        # Whenever the block runs, _executed will be set to True. Also,
        # _latest_input and _output will be set accordingly. In the next run
        # if the input is the same as before, the output will be returned
        # without actually running the block if _executed is set to True.
        self._executed: bool = False
        self._latest_input: Dict[str, Data] = None
        self._output: Dict[str, Data] = dict()

        # Every block will have only one input connection for each input key
        # but can have multiple output connections for each output key. The 
        # connection will have names. Could also be named automatically.
        self._in_connections: Dict[str, _Connection] = dict()
        self._out_connections: Dict[str, List[_Connection]] = dict()
        
        # Every block must be part of a `Pipeline` in order to be able
        # to use its `_auto_execute` method. The `execute` method can 
        # be run independentlyThe parent pipeline as well as the name 
        # of the block in the parent blocks dict must be set before 
        # attempting to _auto_exectute the block.
        self._parent: Pipeline = None
        self._name_in_parent: str = None
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Name of the pipeline block.
        
        In general, this name should be short and abstract. If for displaying
        while execution is in progress a clearer and more complete name is
        needed, the `display_name` property should be set. (defaults to `name`)
        """
        pass

    @property
    def display_name(self) -> str:
        """Name that is printing during execution"""
        return self.name

    @property
    @abstractmethod
    def input_keys(self) -> Tuple[str]:
        """Keys in block's input."""
        pass

    @property
    @abstractmethod
    def output_keys(self) -> Tuple[str]:
        """Keys in block's output dictionary."""
        pass

    @property
    def output(self) -> dict:
        """Output of the latest execution."""
        return self._output

    @output.setter
    def output(self, d: Dict[str, Data]):
        """Setter for `output`.

        Args:
            d (Dict[str, Data]): output dictionary containing output keys 
                and values.
        """
        assert all(isinstance(v, Data) for v in d.values())
        self._output = d

    @property
    def executed(self) -> bool:
        """Whether the block has been previously executed or not."""
        return self._executed
    
    @property
    def base_dir(self) -> str:
        return self._parent.base_dir
    
    @abstractmethod
    def execute(self, **kwargs) -> Dict[str, Data]:
        """Abstract method for executing a block.

        This is the only method that needs defining for a new pipeline block.
        In the implementation, every input should be taken in the form of a 
        argument (optional input keys of the block must have a default in the
        definition of the `execute` method. and at the end a dictionary 
        containing its outputs and their keys should be returned. Each value in 
        the dictionary must be an instance of a subclass of `Data`.

        Returns:
            Dict[str, Data]: Dictionary containing output data.
        """
        pass

    def _auto_execute(self) -> None:
        """Act as a wrapper around `execute` for automatic runs in pipeline."""

        # Pre-execution steps: 1) Is part of pipeline? 2) Get incoming data.
        assert isinstance(self._parent, Pipeline)
        assert isinstance(self._name_in_parent, str)
        block_input: Dict[str, Data] = self._obtain_input_data()
        
        # Do not execute if input is the same as last time and it has been
        # executed previously.
        # TODO: Figure out if repeating exes is because of the second part.
        if self._executed and self._latest_input == block_input:
            # if self.save_output:
            #     self.save() # TODO: Do we need another save here?
            return
        
        # Logging...
        # TODO: Use the logging module of python instead of print
        print()
        print(f"[Running Pipe Block: ({self._name_in_parent}) {self.display_name}]")
        self._carry_exe(block_input)
        self._latest_input = block_input
        self._executed = True
        if self.save_output:
            self.save()

    def _obtain_input_data(self) -> Dict[str, Data]:
        """Pre-execution steps in an automatic run of block in a pipeline."""

        # Obtain incoming data from each of block's connections.
        # NOTE: values' type being `Data` is asserted in `get_data_from_source`
        incoming_input: Dict[str, Data] = {
            name: connection.get_data_from_source() 
            for name, connection in self._in_connections.items()}

        return incoming_input
    
    def _carry_exe(self, block_input: Dict[str, Data]) -> None:
        """Carries the execution of `_auto_execute`"""

        # NOTE: This is separated to be able to modify it in inherited classes.
        self.output: Dict[str, Data] = self.execute(**block_input)
    
    def save(self) -> None:
        """Saves latest output of block if `save_output` is True."""

        block_prefix = os.path.join(self.base_dir, f"{self._name_in_parent}_")
        for output_key, output_data in self.output.items():
            data_prefix = block_prefix + output_key + "_"
            output_data.save(data_prefix=data_prefix)

    def reset(self, removeExecutionHistory: bool = True, 
              removeConnections: bool = True) -> None:
        # TODO: This might need changes depending on intended use.
        
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
        # target_key of the connection in its input keys and there shouldn't 
        # be another connection for its specific key.
        if connection.target_block == self:
            assert connection.target_key in self.input_keys
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

    # def check_input(self, input_dict: Dict[str, Data]):
    #     """Checks validity of input dictionary (input to `execute`).

    #     Args:
    #         input_dict (Dict[str, Any]): input dict with needed input keys.
    #     """

    #     assert all(k in input_dict.keys() for k in 
    #                 self.required_input_keys + self.optional_input_keys)
    #     assert all(isinstance(v, Data) for v in input_dict.values())
    #     return
    
    # def check_output(self, output_dict: Dict[str, Data]):
    #     """Checks validity of output dictionary (output from `execute`).

    #     Args:
    #         output_dict (Dict[str, Any]): output dict with needed output keys.
    #     """
    #     assert set(output_dict.keys()) == set(self.output_keys)
    #     assert all(isinstance(v, Data) for v in output_dict.values())
    #     return

class InputBlock(PipelineBlock):
    name = "Input Block"
    force_run: bool = False

    def __init__(self, keys: List[str], debug: bool = False) -> None:
        super().__init__(debug=debug, save_output=False)
        self._keys = keys.copy()
        self._executed = True
        self._latest_input = dict()
    
    def set_input(self, input_dict: Dict[str, Data]):
        if self.debug:
            assert all(k in self.output_keys for k in input_dict.keys())
        self.output = input_dict.copy()
        # NOTE: A bit of terminology problem here. The output of input block 
        # would be the input of the parent pipeline. input block acts as some 
        # sort of a starting point with no inputs and only output.

    @property
    def input_keys(self):
        return []
    
    @property
    def output_keys(self):
        return self._keys
    
    def execute(self, **kwargs) -> Dict[str, Data]:
        return self.output

class OutputBlock(PipelineBlock):
    name = "Output Block"
    force_run: bool = True

    def __init__(self, keys: List[str], 
                 debug: bool = False, save_output: bool = False) -> None:
        super().__init__(debug=debug, save_output=save_output)
        self._keys = keys
    
    @property
    def input_keys(self):
        return self._keys
    
    @property
    def output_keys(self):
        return self._keys
    
    def execute(self, **kwargs) -> Dict[str, Data]:
        return {k: v for k, v in kwargs.items() if k in self.output_keys}
        # TODO: Why did I have to make this? Test without it.
    
class Pipeline(PipelineBlock, ABC):

    # NOTE: The reason for inheriting from `PipelineBlock` is to be able to use
    # pipelines as blocks for other pipelines as well.

    def __init__(self, base_dir: str = ".", debug: bool = False, 
                 save_output: bool = False, save_all: bool = False):
        super().__init__(debug=debug, save_output=save_output)

        if not os.path.exists(base_dir):
            os.mkdir(base_dir)
        
        self._base_dir = pathlib.Path(base_dir) # TODO: Not use pathlib?
        self._blocks: Dict[str, PipelineBlock | InputBlock | OutputBlock] = {}
        self._connections: Dict[str, _Connection] = dict()

        # Helper for naming unnamed connections and blocks
        self._connect_naming_i = 0
        self._block_naming_i = 0

        # Create the pipeline
        self.build()
        assert "output" in self._blocks.keys()
        assert "input" in self._blocks.keys()
        if save_all:
            for block in self._blocks.values():
                if not isinstance(block, InputBlock) and \
                    not isinstance(block, OutputBlock):
                    block.save_output = True
    
    @property
    def input_keys(self) -> Tuple[str]:
        return self._blocks["input"].output_keys
    
    @property
    def output_keys(self) -> List[str]:
        return self._blocks["output"].output_keys
    
    @property
    def base_dir(self) -> str:
        if self._parent:
            return_dir = os.path.join(self._parent._base_dir, self.name)
        else:
            return_dir = self._base_dir
        
        if not os.path.exists(return_dir):
            os.mkdir(return_dir)
            # TODO: This gets called every time. It shouldn't be so...

        return return_dir

    @abstractmethod
    def build(self):
        pass

    def add_block(
            self, block: PipelineBlock, block_name: str | None = None
        ) -> None:
        """Method for adding single block inside the `build` method."""
        
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
        block._name_in_parent = block_name
        if self.debug:
            block.debug = self.debug
        self._blocks[block_name] = block

    def add_blocks(
            self, blocks: List[PipelineBlock] | Dict[str, PipelineBlock]
        ) -> None:
        """Method for adding multiple blocks inside the `build` method."""

        if isinstance(blocks, list):
            [self.add_block(block) for block in blocks]
        if isinstance(blocks, dict):
            [self.add_block(block, name) for name, block in blocks.items()]
    
    def add_connection(self, source_name: str, source_key: str, 
                       target_name: str, target_key: str,
                       connection_name: str | None = None) -> None:
        """Method for adding single connection inside the `build` method. """
        
        # BAD NAMING WARNING: This method is completely different from the
        # private `_add_connection` method inherited from `PipelineBlock`.
        # This adds a new connection between blocks inside the pipeline while 
        # the inherited one helps establish connections from outside if the
        # pipeline is also acting as a block of another pipeline.

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

    def execute(self, **pipe_input) -> Dict[str, Data]:
        """Execute method for a pipeline which auto executes its blocks."""

        l = 55
        print()
        print(">"*l)
        print(">>" + str.center("STARTED: " + self.display_name, l-5), ">>")
        print(">"*l)

        self._blocks["input"].set_input(pipe_input)

        for block in self._blocks.values():
            if block.force_run:
                block._auto_execute()

        print()
        print("<"*l)
        print("<<" + str.center("ENDED: " + self.display_name, l-5), "<<")
        print("<"*l)

        output = self._blocks["output"]._output
        return output

    def reset(self):
        raise NotImplementedError() 
        # TODO: Implement after making sure about block reset.
    
    
    
    
    
    
    
    