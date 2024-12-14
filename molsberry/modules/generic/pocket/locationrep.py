import numpy as np

from rdkit import Chem

from ....core import Representation

class PocketLocation:
    """Handles pocket location specifications for different systems.

    This class can convert between different types of pocket location formats:
    - point:    Point and Radius, 
    - box2p:    3D Box With Two Opposite Corner Points,
    - boxcxyz:  3D Box Given Center Point and Sizes in X, Y, Z Dimensions,
    - ligand:   Ligand and Radius,
    - residues: Protein Residues
    """

    def __init__(self, method, **kwargs): 
        
        self.params = {}
        if method == "point":
            self._init_point(**kwargs)
        elif method == "box2p":
            self._init_box2p(**kwargs)
        elif method == "boxcxyz":
            self._init_boxcxyz(**kwargs)
        elif method == "ligand":
            self._init_ligand(**kwargs)
        elif method == "residues":
            self._init_residues(**kwargs)
        else:
            raise ValueError("Method: ('point', 'box', 'ligand', 'residues')")
        self.loc_type = method
    
    def point_is_included(self, point):
        if self.loc_type == "point":
            return self._point_is_included_point(point)
        elif self.loc_type == "ligand":
            return self._point_is_included_ligand(point)
        else:
            raise NotImplementedError() # TODO: Write other types...
    
    def get_params(self, param_type: str):
        if param_type not in self.params.keys():
            for format in self.params.keys():
                try:
                    self.loctype_conversion(format, param_type)
                    return self.params[param_type]
                except NotImplementedError:
                    continue
            raise ValueError(f"Could not convert existing formats to {param_type}")
        else:
            return self.params[param_type]
        
    def loctype_conversion(self, source_format, target_format):
        """
        Converts the pocket location to the target format.

        Args:
            target_format (str): The target format ('point', 'box', or 'ligand').
        """
        if source_format not in self.params.keys():
            raise ValueError("Source format is not defined!")
        if source_format == "ligand" and target_format == "boxcxyz":
            coords = self.params["ligand"]["c"]
            radius = self.params["ligand"]["r"]
            lig_center = (np.max(coords, axis=0) + np.min(coords, axis=0)) / 2
            size = np.max(coords, axis=0) - np.min(coords, axis=0) + radius
            self.params["boxcxyz"] = {"c": lig_center, "size": size}
        elif source_format == "point" and target_format == "boxcxyz":
            point = self.params["point"]["p"]
            radius = self.params["point"]["r"]
            self.params["boxcxyz"] = {"c": point, "size": np.array([radius]*3)}
        else:
            raise NotImplementedError()
        # TODO: Implement other conversions between pocket location formats

    def _init_point(self, point=None, radius=10):
        point = np.array(point)
        self.params["point"] = {"p": point, "r": radius}

    def _init_box2p(self, point1, point2):
        pass

    def _init_boxcxyz(self, center, xyz_size):
        self.params["boxcxyz"] = {"c": np.array(center), "size": np.array(xyz_size)}

    def _init_ligand(self, ligand, radius=10):
        if isinstance(ligand, str):
            if ligand.endswith(".pdb"):
                ligand = Chem.MolFromPDBFile(ligand)
            elif ligand.endswith(".mol2"):
                ligand = Chem.MolFromMol2File(ligand)
            elif ligand.endswith(".sdf"):
                ligand = next(Chem.SDMolSupplier(ligand))
            else:
                raise ValueError("ligand has to be .pdb, .mol2 or .sdf...")
        elif isinstance(ligand, Chem.Mol):
            pass
        else:
            raise ValueError("ligand has to be .pdb, .mol2 or .sdf or an rdmol instance...")
        
        coords = ligand.GetConformer().GetPositions()
        self.params["ligand"] = {
            "c": coords, "r": radius, "l": ligand
        }

    def _init_residues(self, protein_pdb, res_ids):
        pass

    def _point_is_included_point(self, point):
        point = np.array(point)
        if np.linalg.norm(point - self.params["point"]["p"]) <= self.params["point"]["r"]:
            return True
        return False
    
    def _point_is_included_ligand(self, point):
        point = np.array(point)
        diff = self.params["ligand"]["c"] - point
        if np.min(np.linalg.norm(diff, axis=1)) <= self.params["ligand"]["r"]:
            return True
        return False
    
class PocketLocationRep(Representation):
    rep_name = "poclocrep"
    def __init__(self, pocket_loc: PocketLocation):
        super().__init__(content=pocket_loc)