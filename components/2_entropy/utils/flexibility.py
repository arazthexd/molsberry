from rdkit import Chem
from rdkit.Chem import RingInfo, rdmolops

def CalcRingFlexibility(mol):

    rdmolops.FindRingFamilies(mol)
    ring_info: Chem.RingInfo = mol.GetRingInfo()
    ring_atoms_list = ring_info.AtomRingFamilies()
    num_rings = len(ring_atoms_list)
    ring_sizes = [len(ring_atoms) for ring_atoms in ring_atoms_list]
    penalties_per_ring = {
        "num_fused_single": [0] * num_rings,
        "num_nonsingle": [0] * num_rings,
        "num_spiro": [0] * num_rings,
        "num_bridgehead": [0] * num_rings, # TODO: How? What?
        "num_polycyclic": [0] * num_rings, # TODO: How? What?
        "num_specifics": [0] * num_rings,
        # "num_cyclicamide": [0] * num_rings,
        # "num_cyclicthioamide": [0] * num_rings,
    }

    specifics_query = Chem.MolFromSmarts("[C;R](=[O,S])-[O,N;R]")
    specifics_match = mol.GetSubstructMatches(specifics_query)
    for i, ring_atoms1 in enumerate(ring_atoms_list):

        ring_mol = Chem.MolFromSmarts(Chem.MolFragmentToSmarts(mol, ring_atoms1))
        penalties_per_ring["num_nonsingle"][i] = sum([Chem.BondType.SINGLE != bond.GetBondType() for bond in ring_mol.GetBonds()])

        for j, ring_atoms2 in enumerate(ring_atoms_list[i+1:]):
            common = list(set(ring_atoms1).intersection(ring_atoms2))

            if len(common) == 1:
                penalties_per_ring["num_spiro"][i] += 1
                penalties_per_ring["num_spiro"][i+j+1] += 1

            elif len(common) == 2:
                commons_bond = mol.GetBondBetweenAtoms(common[0], common[1])
                if commons_bond.GetBondType() == Chem.BondType.SINGLE:
                    penalties_per_ring["num_fused_single"][i] += 1
                    penalties_per_ring["num_fused_single"][i+j+1] += 1

            elif len(common) > 2:
                penalties_per_ring["num_polycyclic"][i] += 1
                penalties_per_ring["num_polycyclic"][i+j+1] += 1
        
        for match in specifics_match:
            if match[0] in ring_atoms1 and match[2] in ring_atoms1:
                penalties_per_ring["num_specifics"][i] += 1
                break # TODO: Two esters >> p=2 or p=1?

    total_penalties = [sum(ps[i] for ps in penalties_per_ring.values()) for i in range(num_rings)]
    flexibility = sum(max(ring_size-3-total_penalties[i], 0) for i, ring_size in enumerate(ring_sizes))

    return flexibility