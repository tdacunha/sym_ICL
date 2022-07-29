import numpy as np
import astropy.table

def read_ELVES_hosts():
    float_cols = np.loadtxt("tables/ELVES_hosts.txt",
                            delimiter="\t", skiprows=6,
                            usecols=np.arange(1,9, dtype=int)).T
    string_cols = np.loadtxt("tables/ELVES_hosts.txt", delimiter="\t",
                             skiprows=6, usecols=(0, 9, 10), dtype=str).T

    name_len = max(map(len, string_cols[0]))
    source_len = max(map(len, string_cols[1]))
    ref_len = max(map(len, string_cols[2]))
    
    DTYPE = np.dtype([
        ("name", "S%d" % name_len), ("dist", "f8"), ("vr", "f8"),
        ("M_Ks", "f8"), ("M_Ks_group", "f8"), ("M_V", "f8"), ("B-V", "f8"),
        ("Ms", "f8"), ("r_cover", "f8"), ("source", "S%d" % source_len),
        ("ref", "S%d" % ref_len)
    ])

    out = np.zeros(len(float_cols[0]), dtype=DTYPE)
    out["name"] = string_cols[0]
    out["source"] = string_cols[1]
    out["ref"] = string_cols[2]
    out["dist"] = float_cols[0]
    out["vr"] = float_cols[1]
    out["M_Ks"] = float_cols[2]
    out["M_Ks_group"] = float_cols[3]
    out["M_V"] = float_cols[4]
    out["B-V"] = float_cols[5]
    out["Ms"] = float_cols[6]
    out["r_cover"] = float_cols[7]
    
    return out

def read_ELVES_sats():
    names = ["name", "host", "ra", "dec", "r_proj", "g", "g_err", "r", "r_err", "M_g", "M_V", "M_V_err", "Ms", "log_Ms_err", "mu_V", "mu_V_err", "r_eff", "r_eff_err", "early_type", "confirmed", "p_sat", "filter", "telescope", "bad_photo"]
    string_fields = ["name", "host", "filter", "telescope", "bad_photo"]
    bool_fields = ["early_type", "confirmed"]
    col_starts = np.array([1, 16, 24, 33, 42, 48, 54, 59, 65, 70, 77, 84, 89, 95, 100, 106, 111, 119, 125, 131, 137, 142, 147, 157], dtype=int)-1
    col_ends = np.array([14, 22, 31, 40, 46, 52, 57, 63, 68, 75, 82, 87, 93, 98, 104, 109, 117, 123, 129, 135, 140, 145, 155, 159], dtype=int)-1
    
    n = len(names)
    assert(len(col_ends) == n)
    assert(len(col_starts) == n)

    DTYPE = []
    for i in range(n):
        if names[i] in string_fields:
            DTYPE.append((names[i], "S%d" % (col_ends[i]-col_starts[i]+1)))
        elif names[i] in bool_fields:
            DTYPE.append((names[i], "?"))
        else:
            DTYPE.append((names[i], "f8"))
    DTYPE.append(("host_idx", "i8"))
            
    with open("tables/ELVES_sats.txt", "r") as f: text = f.read()
    rows = text.split("\n")
    rows = [row for row in rows[38:] if len(row.strip()) > 0]

    out = np.zeros(len(rows), dtype=DTYPE)

    host_idx = 0
    for ri in range(len(rows)):
        if ri != 0 and rows[ri][1] != rows[ri][1]:
            host_idx += 1
        for ci in range(n):
            token = rows[ri][col_starts[ci]: col_ends[ci]+1].strip()
            if names[ci] == "host_idx":
                out[ri][names[ci]] = host_idx
            elif names[ci] in string_fields:
                out[ri][names[ci]] = token
            elif names[ci] in bool_fields:
                out[ri][names[ci]] = token == "True"
            else:
                if len(token) == 0:
                    out[ri][names[ci]] = np.nan
                else:
                    out[ri][names[ci]] = float(token)

    out["Ms"] = 10**out["Ms"]
    return out
        
def read_SAGA():
    DTYPE = [
        ("host", "S40"), ("name", "S40"), ("photo_coverage", "S3"), ("z_source", "S40"),
        ("host_idx", "i8"),
        ("H_alpha", "?"),
        ("ra", "f8"), ("dec", "f8"), ("host_ra", "f8"), ("host_dec", "f8"), ("r_proj", "f8"), ("m_r", "f8"), ("M_r", "f8"), ("Ms", "f8"), ("dist", "f8"), ("vr_host", "f8"), ("d_vr", "f8"), ("mu_eff", "f8"), ("M_K_host", "f8"), ("spec_coverage", "f8"), ("g-r", "f8")
    ]

    # Taken almost verbatim from the SAGA website's example code.
    saga_hosts = astropy.table.Table.read("tables/saga_stage2_hosts.csv")
    saga_sats = astropy.table.Table.read("tables/saga_stage2_sats.csv")

    saga_joined = astropy.table.join(
        saga_sats, saga_hosts,
        keys="INTERNAL_HOSTID",
        join_type="left",
        uniq_col_name="{table_name}{col_name}",
        table_names=["", "HOST_"]
    )

    out = np.zeros(len(saga_joined), dtype=DTYPE)
    out["host"] = saga_joined["HOST_COMMON_NAME"]
    out["name"] = saga_joined["OBJID"]
    out["ra"] = saga_joined["RA"]
    out["dec"] = saga_joined["DEC"]
    out["r_proj"] = saga_joined["D_PROJ"]
    out["d_vr"] = saga_joined["DELTA_HRV"]
    out["m_r"] = saga_joined["R"]
    out["M_r"] = saga_joined["R_ABS"]
    out["g-r"] = saga_joined["GR"]
    out["mu_eff"] = saga_joined["MU_EFF"]
    out["Ms"] = 10**saga_joined["LOG_STELLAR_MASS"]
    out["H_alpha"] = [x == "Y" for x in saga_joined["H_ALPHA"]]
    out["z_source"] = saga_joined["Z_SOURCE"]
    out["host_ra"] = saga_joined["HOST_RA"]
    out["host_dec"] = saga_joined["HOST_DEC"]
    out["vr_host"] = saga_joined["HRV"]
    out["M_K_host"] = saga_joined["M_K"]
    out["dist"] = saga_joined["DIST"]
    out["photo_coverage"] = saga_joined["PHOT_COVERAGE"]
    out["spec_coverage"] = saga_joined["SPEC_COVERAGE"]
    
    host_idx = 0
    for i in range(len(out)):
        if i != 0 and out["host"][i] != out["host"][i-1]:
            host_idx += 1
        out["host_idx"][i] = host_idx
    
    return out
    
if __name__ == "__main__":
    hosts = read_ELVES_hosts()
    sats = read_ELVES_sats()
    saga = read_SAGA()
