import os, requests, time
from typing import Tuple

from rdkit import Chem
from rdkit.Chem import Mol

from ..abstract import ComplexConverter
from ..utils.iotools import save_pl_complex

class ProtossProtoTautoDefiner(ComplexConverter):
    base_url = "https://proteins.plus/api"
    pdb_endpoint = "/pdb_files_rest"
    protoss_endpoint = "/protoss_rest"

    def __init__(self, work_dir: str = "."):
        self.work_dir = work_dir

    def convert(self, ligand: Mol, target_path: str) -> Tuple[Mol, str]:
        input_path = os.path.join(self.work_dir, "tmp_protoss_input.pdb")
        save_pl_complex(ligand, target_path, input_path)
        response = requests.post("https://proteins.plus/api/pdb_files_rest", 
                                 headers={"Accept": "application/json"},
                                 files={"pdb_file[pathvar]": open(input_path, "rb")})

        for i in range(5):
            res2 = requests.get(response.json()["location"])
            if res2.status_code == 200:
                break
            time.sleep(1)

        pdb_id = res2.json()["id"]
        response = requests.post(self.base_url+self.protoss_endpoint,
                                 headers={"Accept": "application/json"},
                                 json={"protoss": {"pdbCode": res2.json()["id"]}})
        
        for i in range(5):
            res2 = requests.get(response.json()["location"])
            if res2.status_code == 200:
                break
            time.sleep(1)
        
        target_pdb_block = requests.get(res2.json()["protein"]).text
        ligand_pdb_block = requests.get(res2.json()["ligands"])



        


