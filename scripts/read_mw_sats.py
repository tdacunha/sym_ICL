import numpy as np

def read():
    with open("tables/mw_sats.txt") as f: full_text = f.read()
    lines = full_text.split("\n")
    tokens = [[tok.strip() for tok in line.split("&")] for line in lines
              if len(line) > 0]

    max_name_size = max([len(line[0]) for line in tokens])

    dtype = [("name", "S%d" % max_name_size), ("r", np.float64),
             ("r_half", np.float64), ("MV", np.float64),
             ("class", np.int64)]

    out = np.zeros(len(tokens), dtype=dtype)

    for i, line in enumerate(tokens):
        out["name"][i] = tokens[i][0]
        out["class"][i] = int(tokens[i][2])
        out["r"][i] = float(tokens[i][8])
        out["r_half"][i] = float(tokens[i][9])
        out["MV"][i] = float(tokens[i][10])

    return out

if __name__ == "__main__": read()
