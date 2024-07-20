def dissociate_complex(complex_path, lig_name, lig_write, poc_write):
    compl = open(complex_path, "r")
    lig = open(lig_write, "w")
    poc = open(poc_write, "w")
    for line in compl.readlines():
        if lig_name in line:
            lig.write(line)
        else:
            poc.write(line)
    lig.write("END\n")
    [x.close() for x in [compl, lig, poc]]